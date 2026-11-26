HERE <- file.path(Sys.getenv("DTS_BUNDLE"),"pif")
source(file.path(HERE, "00_config.R")); source(file.path(HERE, "01_helpers.R"))
need(c("data.table", "Matrix", "statmod")); data.table::setDTthreads(1L)
save_session("PIF")
if (!file.exists(file.path(DIRS$prepared, "PREPARED_OK.rds"))) stop("Run 02_prepare.R successfully first.")
metric <- Sys.getenv("DTS_HEAT", "WBGT")
if (!metric %in% HEAT$metric) stop("Unknown DTS_HEAT.")
j <- match(metric, HEAT$metric); h <- HEAT[j, ]
m <- readRDS(file.path(DIRS$prepared, "GM3_compact_parameters.rds"))[[metric]]
if (!valid_diag(m$diagnostic)) stop("Invalid primary GM3.")
final_path <- file.path(DIRS$tables, paste0(h$stem, "_PIF_results.csv"))
if (file.exists(final_path)) stop("Completed PIF output already exists: ", final_path)
pred <- data.table::as.data.table(readRDS(file.path(DIRS$prepared, "PIF_common_full_trajectories.rds")))
keep <- unique(c("preg_uid", "mother_uid", COV_VARS, h$scaled, h$raw, h$monthly, h$iqr, ".g"))
pred <- pred[, ..keep]
pred <- align_to_model(pred, m)
data.table::setorder(pred, preg_uid, .g)
H <- length(m$horizon); N <- data.table::uniqueN(pred$preg_uid)
if (nrow(pred) != N * H) stop("Invalid full trajectory size.")
g <- matrix(pred$.g, ncol = H, byrow = TRUE)
if (!all(sweep(g, 2, m$horizon, "-") == 0)) stop("Trajectories have unexpected month order.")
ord <- rep(seq_len(H), N)
# Check the unchanged natural cumulative exposure reconstruction exactly.
natural <- pred[, {
 z <- get(h$monthly)[seq_len(3L)]; cs <- cumsum(z)/seq_len(3L)
 .(reconstructed = c(cs, rep(cs[3], .N - 3L)))
}, by = preg_uid]$reconstructed
if (max(abs(natural - pred[[h$raw]])) > NUMERIC_TOL) stop("Natural cumulative exposure does not match months 1-3; inspect upstream data.")
# Country cutpoints: one completed T1 mean per pregnancy, including reconstructed
# histories after early outcomes. Cutpoints are fixed for uncertainty estimation.
t1 <- pred[ord == 3L, .(country = as.character(country), T1 = get(h$raw), IQR = get(h$iqr))]
thresholds <- t1[, .(p90 = quantile(T1, .90, type = 7), p75 = quantile(T1, .75, type = 7),
 p50 = quantile(T1, .50, type = 7), IQR = IQR[1], n_distinct_IQR = data.table::uniqueN(IQR)), by = country]
floors <- pred[ord <= 3L, .(floor = min(get(h$monthly))), by = .(country = as.character(country))]
thresholds <- merge(thresholds, floors, by = "country")
if (any(thresholds$n_distinct_IQR != 1L)) stop("IQR denominators differ within country.")
write_csv(thresholds, file.path(DIRS$qc, paste0(h$stem, "_PIF_intervention_thresholds.csv")))
map <- match(as.character(pred$country), thresholds$country)
cf <- matrix(pred[[h$scaled]], nrow = nrow(pred), ncol = length(SCENARIOS),
 dimnames = list(NULL, SCENARIOS))
intervention_qc <- list()
for (scenario in SCENARIOS[-1]) {
  z <- pred[[h$monthly]]
  active <- ord <= 3L
  if (grepl("^cap", scenario)) {
    cutoff <- thresholds[[sub("^cap_", "", scenario)]][map]
    z[active] <- pmin(z[active], cutoff[active])
  } else {
    z[active] <- pmax(z[active] - pred[[h$iqr]][active], thresholds$floor[map[active]])
  }
  zmat <- matrix(z, ncol = H, byrow = TRUE)
  cum <- zmat
  for (k in 1:3) cum[, k] <- rowMeans(zmat[, 1:k, drop = FALSE])
  if (H > 3L) for (k in 4:H) cum[, k] <- cum[, 3]
  rawcf <- as.vector(t(cum))
  cf[, scenario] <- rawcf / pred[[h$iqr]]
  lo <- tapply(pred[[h$raw]], as.character(pred$country), min)
  out_of_range <- rawcf < lo[as.character(pred$country)] - NUMERIC_TOL
  intervention_qc[[scenario]] <- data.frame(metric = metric, scenario = scenario,
    changed_first_trimester_months = sum(abs(z[active] - pred[[h$monthly]][active]) > NUMERIC_TOL),
    pregnancies_changed = sum(rowSums(abs(cum - matrix(pred[[h$raw]], ncol = H, byrow = TRUE)) > NUMERIC_TOL) > 0),
    prediction_rows_below_country_observed_cumulative_min = sum(out_of_range),
    monthly_lower_bound_rule = "Country minimum T1 monthly exposure; no below-support monthly shift")
}
write_csv(data.table::rbindlist(intervention_qc), file.path(DIRS$qc, paste0(h$stem, "_PIF_intervention_QC.csv")))

# Standardized NET risk, integrating a common mother random intercept over its
# estimated N(0, sigma^2) distribution. Setting b=0 is NOT marginalization.
# R_i = E_b[1 - exp(-exp(b) * sum_t exp(X_it beta))].
# CIs use joint beta/log-SD covariance and a delta method on log(R_cf/R_nat).
# They condition on empirical covariates/cutpoints and are NOT bootstrap CIs.
gh_n <- as.integer(Sys.getenv("PIF_GH_N", "50"))
if (is.na(gh_n) || gh_n < 5L) stop("PIF_GH_N must be >=5.")
quad <- statmod::gauss.quad.prob(gh_n, dist = "normal")
quad_check <- statmod::gauss.quad.prob(gh_n + 10L, dist = "normal")
sigma <- exp(m$theta); p <- length(m$beta)
chunk_n <- as.integer(Sys.getenv("PIF_CHUNK_PREGNANCIES", "5000"))
if (is.na(chunk_n) || chunk_n < 1L) stop("Invalid chunk size.")
countries <- sort(unique(as.character(pred$country)))
groups <- c("Pooled", countries)
risk_sum <- matrix(0, length(groups), length(SCENARIOS), dimnames = list(groups, SCENARIOS))
grad_sum <- array(0, c(length(groups), length(SCENARIOS), p + 1L))
curve_sum <- matrix(0, H, length(SCENARIOS), dimnames = list(as.character(m$horizon), SCENARIOS))
curve_grad <- array(0, c(H, length(SCENARIOS), p + 1L))
n_group <- setNames(c(N, as.numeric(table(factor(t1$country, levels = countries)))), groups)
max_quad_difference <- 0
marginal <- function(A, quad) {
  R <- rep(0, length(A)); dA <- R; dt <- R
  for (l in seq_along(quad$nodes)) {
    mult <- exp(sigma * quad$nodes[l]); B <- mult * A
    surv <- exp(-B)
    R <- R + quad$weights[l] * (-expm1(-B))
    dA <- dA + quad$weights[l] * surv * mult
    dt <- dt + quad$weights[l] * surv * B * sigma * quad$nodes[l]
  }
  list(risk = R, dA = dA, dtheta = dt)
}
for (start in seq.int(1L, N, by = chunk_n)) {
  end <- min(N, start + chunk_n - 1L)
  ix <- seq.int((start - 1L) * H + 1L, end * H)
  d <- pred[ix]; n <- end - start + 1L
  X0 <- design_matrix(d, m)
  eta0 <- as.numeric(X0 %*% m$beta)
  ci <- match(as.character(d$country[seq.int(1L, nrow(d), by = H)]), countries) + 1L
  for (s in seq_along(SCENARIOS)) {
    delta <- cf[ix, s] - d[[h$scaled]]
    eta <- eta0 + m$beta[[h$scaled]] * delta
    if (any(!is.finite(eta)) || max(eta) > 700) stop("Unstable counterfactual linear predictor.")
    rate <- exp(eta)
    X <- X0
    X[, h$scaled] <- cf[ix, s]
    # One small matrix per chunk. No full N x N prediction covariance.
    WX <- as.matrix(X) * rate
    for (k in seq_len(H)) {
      idx <- which(rep(seq_len(H), n) <= k)
      ids <- rep(seq_len(n), each = H)[idx]
      A <- as.numeric(rowsum(rate[idx], ids, reorder = FALSE))
      dAb <- rowsum(WX[idx, , drop = FALSE], ids, reorder = FALSE)
      rr <- marginal(A, quad)
      dg <- cbind(dAb * rr$dA, rr$dtheta)
      curve_sum[k, s] <- curve_sum[k, s] + sum(rr$risk)
      curve_grad[k, s, ] <- curve_grad[k, s, ] + colSums(dg)
      if (k == H) {
        sums <- rowsum(cbind(rr$risk, dg), ci, reorder = FALSE)
        for (a in seq_len(nrow(sums))) {
          gi <- as.integer(rownames(sums)[a])
          risk_sum[gi, s] <- risk_sum[gi, s] + sums[a, 1]
          grad_sum[gi, s, ] <- grad_sum[gi, s, ] + sums[a, -1]
        }
        risk_sum[1, s] <- risk_sum[1, s] + sum(rr$risk)
        grad_sum[1, s, ] <- grad_sum[1, s, ] + colSums(dg)
        rr2 <- marginal(A, quad_check)$risk
        max_quad_difference <- max(max_quad_difference, max(abs(rr2 - rr$risk)))
      }
    }
  }
  if (start == 1L || end == N || start %% (10L * chunk_n) == 1L) message("PIF ", metric, ": ", end, "/", N, " pregnancies")
}
risks <- risk_sum / n_group
for (a in seq_along(groups)) grad_sum[a, , ] <- grad_sum[a, , ] / n_group[a]
curve_sum <- curve_sum/N; curve_grad <- curve_grad/N
if (max_quad_difference > 1e-5) stop("Quadrature check exceeded tolerance; inspect diagnostics and increase PIF_GH_N.")
tables <- list()
for (a in seq_along(groups)) for (s in 2:length(SCENARIOS)) {
  r0 <- risks[a, 1]; r1 <- risks[a, s]
  if (min(r0, r1) <= 0) stop("Nonpositive standardized risk.")
  gL <- grad_sum[a, s, ]/r1 - grad_sum[a, 1, ]/r0
  varL <- as.numeric(t(gL) %*% m$covariance %*% gL)
  if (!is.finite(varL) || varL < -1e-12) stop("Invalid PIF variance.")
  seL <- sqrt(max(varL, 0)); L <- log(r1/r0)
  tables[[length(tables) + 1L]] <- data.frame(metric = metric, population = groups[a], scenario = SCENARIOS[s],
    pregnancies = n_group[a], natural_net_risk = r0, counterfactual_net_risk = r1,
    risk_difference_per_1000 = 1000 * (r0 - r1), PIF_percent = 100 * (1 - exp(L)),
    lower_95 = 100 * (1 - exp(L + qnorm(.975)*seL)), upper_95 = 100 * (1 - exp(L - qnorm(.975)*seL)),
    SE_log_risk_ratio = seL, inference = "Joint beta/log-SD delta method; fixed empirical target and cutpoints",
    estimand = "Model-based net risk; binary-comparator target; competing events eliminated")
}
curves <- data.table::rbindlist(lapply(seq_along(SCENARIOS), function(s) {
  se <- vapply(seq_len(H), function(k) sqrt(max(0, as.numeric(t(curve_grad[k,s,]) %*% m$covariance %*% curve_grad[k,s,]))), numeric(1))
  data.frame(metric = metric, scenario = SCENARIOS[s], gestational_interval = m$horizon,
    net_risk = curve_sum[,s], lower_95 = pmax(0, curve_sum[,s]-qnorm(.975)*se),
    upper_95 = pmin(1, curve_sum[,s]+qnorm(.975)*se))
}))
write_csv(curves, file.path(DIRS$tables, paste0(h$stem, "_PIF_net_risk_curves.csv")))
write_csv(data.frame(metric = metric, GH_nodes = gh_n, check_nodes = gh_n + 10L,
 max_absolute_risk_difference = max_quad_difference, mother_SD = sigma, pregnancies = N,
 horizon_labels = paste(m$horizon, collapse = ";")), file.path(DIRS$qc, paste0(h$stem, "_PIF_numerical_diagnostics.csv")))
save_rds(list(model = m, risks = risks, gradients = grad_sum, curve_gradients = curve_grad,
 thresholds = thresholds, quadrature = quad), file.path(DIRS$models, paste0(h$stem, "_PIF_calculation_checkpoint.rds")))
write_csv(data.table::rbindlist(tables), final_path)
message("PIF completed: ", final_path)
