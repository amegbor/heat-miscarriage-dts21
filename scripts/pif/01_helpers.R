# This file defines functions only; entry scripts load config separately.
need <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Install required R packages first: ", paste(missing, collapse = ", "))
}
assert_cols <- function(d, cols) {
  missing <- setdiff(cols, names(d))
  if (length(missing)) stop("Missing required columns: ", paste(missing, collapse = ", "))
}
write_csv <- function(d, file) data.table::fwrite(d, file, na = "NA")
safe_num <- function(expr) tryCatch(as.numeric(expr), error = function(e) NA_real_)
save_session <- function(label) writeLines(capture.output(sessionInfo()), file.path(DIRS$qc, paste0(label, "_sessionInfo.txt")))
save_rds <- function(x, path) {
  tmp <- paste0(path, ".tmp_", Sys.getpid())
  saveRDS(x, tmp, compress = FALSE)
  if (!file.rename(tmp, path)) stop("Could not save checkpoint: ", path)
}
counts <- function(d, stage) data.frame(stage = stage, person_months = nrow(d),
 pregnancies = data.table::uniqueN(d$preg_uid), mothers = data.table::uniqueN(d$mother_uid),
 events = if ("miscarriage_event" %in% names(d)) sum(d$miscarriage_event) else NA_integer_)
complete_finite <- function(d, cols) {
  assert_cols(d, cols)
  ok <- complete.cases(d[, cols, with = FALSE])
  for (v in cols) if (is.numeric(d[[v]])) ok <- ok & is.finite(d[[v]])
  ok
}
fixed_formula <- function(expvars) reformulate(c(expvars, COV_TERMS), response = "miscarriage_event")
drop_factors <- function(d) {
  for (v in names(d)) if (is.factor(d[[v]])) data.table::set(d, j = v, value = droplevels(d[[v]]))
  d
}
check_risk <- function(d, full = FALSE, horizon = NULL) {
  assert_cols(d, c("preg_uid", "mother_uid", "country", "gest_month_f"))
  g <- suppressWarnings(as.integer(as.character(d$gest_month_f)))
  if (anyNA(g)) stop("gest_month_f must have integer month labels.")
  d[, .g := g]
  if (anyDuplicated(d[, .(preg_uid, .g)])) stop("Duplicate pregnancy-month records.")
  check <- d[, .(mothers = data.table::uniqueN(mother_uid), countries = data.table::uniqueN(country),
    gaps = any(diff(sort(.g)) != 1L)), by = preg_uid]
  if (any(check$mothers != 1 | check$countries != 1)) stop("Inconsistent pregnancy identifiers.")
  if ("miscarriage_event" %in% names(d)) {
    if (anyNA(d$miscarriage_event) || any(!d$miscarriage_event %in% 0:1)) stop("Invalid event indicator.")
    ck <- d[, .(n = sum(miscarriage_event), after = any(.g > max(.g[miscarriage_event == 1], -Inf))), by = preg_uid]
    if (any(ck$n > 1 | (ck$n == 1 & ck$after))) stop("Multiple events or records after an event.")
  }
  if (full) {
    ck <- d[, .(valid = identical(sort(.g), as.integer(horizon))), by = preg_uid]
    if (any(!ck$valid)) stop("Prediction dataset has incomplete/incorrect horizons; do not fill after loss by carrying exposures forward.")
  }
  invisible(check)
}
fit_diag <- function(fit, label) {
  b <- glmmTMB::fixef(fit)$cond
  V <- tryCatch(vcov(fit)$cond, error = function(e) matrix(NA_real_, length(b), length(b)))
  se <- suppressWarnings(sqrt(diag(V)))
  grad <- tryCatch(fit$sdr$gradient.fixed, error = function(e) NULL)
  if (is.null(grad) || !length(grad)) grad <- tryCatch(fit$obj$gr(fit$fit$par), error = function(e) NULL)
  data.frame(model = label, convergence_code = fit$fit$convergence,
    optimizer_message = paste(fit$fit$message, collapse = "; "), pdHess = isTRUE(fit$sdr$pdHess),
    all_fixed_SE_finite = length(se) == length(b) && all(is.finite(se)),
    all_fixed_estimates_finite = all(is.finite(b)),
    max_abs_gradient = if (length(grad) && all(is.finite(grad))) max(abs(grad)) else NA_real_,
    mother_random_intercept_SD = safe_num(attr(glmmTMB::VarCorr(fit)$cond$mother_uid, "stddev"))[1],
    logLik = safe_num(logLik(fit)), AIC = safe_num(AIC(fit)), BIC = safe_num(BIC(fit)),
    n_parameters = safe_num(attr(logLik(fit), "df")), N = nobs(fit),
    family = family(fit)$family, link = family(fit)$link)
}
valid_diag <- function(d) d$convergence_code == 0 && d$pdHess && d$all_fixed_SE_finite &&
 d$all_fixed_estimates_finite && all(is.finite(c(d$logLik, d$AIC))) && d$link == "cloglog"
save_fixed <- function(b, V, prefix) {
  se <- suppressWarnings(sqrt(diag(V))); z <- b / se
  write_csv(data.frame(term = names(b), estimate = b, SE = se, HR = exp(b),
   lower_95 = exp(b - qnorm(.975) * se), upper_95 = exp(b + qnorm(.975) * se),
   p_value = 2 * pnorm(abs(z), lower.tail = FALSE)), paste0(prefix, "_fixed_effects.csv"))
  save_rds(V, paste0(prefix, "_covariance.rds"))
}
compact_model <- function(fit, metric) {
  # Guard against silently using an older/different substantive model.
  ff <- lme4::nobars(formula(fit))
  expected <- fixed_formula(HEAT$scaled[match(metric, HEAT$metric)])
  if (!setequal(attr(terms(ff), "term.labels"), attr(terms(expected), "term.labels"))) stop("Unexpected GM3 fixed formula: ", metric)
  bars <- lme4::findbars(formula(fit))
  if (length(bars) != 1 || paste(deparse(bars[[1]]), collapse = "") != "1 | mother_uid") stop("Expected only a mother random intercept.")
  if (length(fit$modelInfo$reTrms$cond$cnms) != 1L) stop("Unexpected random effects.")
  b <- glmmTMB::fixef(fit)$cond
  par <- fit$fit$par
  ib <- which(names(par) == "beta"); it <- which(names(par) == "theta")
  if (length(ib) != length(b) || length(it) != 1L || length(par) != length(b) + 1L) stop("Unexpected top-level parameter layout.")
  # sdr covariance has the optimizer's top-level parameter order.
  VV <- fit$sdr$cov.fixed
  if (is.null(VV) || nrow(VV) != length(par) || any(!is.finite(VV))) stop("Invalid joint beta/theta covariance.")
  V <- VV[c(ib, it), c(ib, it), drop = FALSE]
  if (max(abs(V[seq_along(b), seq_along(b), drop = FALSE] - vcov(fit)$cond)) > 1e-8) stop("Joint covariance ordering differs from conditional covariance.")
  dimnames(V) <- list(c(names(b), "log_mother_SD"), c(names(b), "log_mother_SD"))
  if (min(eigen(V, symmetric = TRUE, only.values = TRUE)$values) <= 0) stop("Joint covariance is not positive definite.")
  xlevels <- lapply(fit$frame[vapply(fit$frame, is.factor, logical(1))], levels)
  # Preserve contrasts from the fitted conditional design matrix.
  X <- model.matrix(fit, component = "cond")
  ans <- list(metric = metric, beta = b, theta = unname(par[it]), covariance = V,
    formula = ff, xlevels = xlevels, contrasts = attr(X, "contrasts"),
    horizon = sort(as.integer(levels(fit$frame$gest_month_f))), diagnostic = fit_diag(fit, metric))
  rm(X)
  ans
}
align_to_model <- function(d, m) {
  for (v in intersect(names(m$xlevels), names(d))) {
    vals <- as.character(d[[v]])
    if (any(!vals %in% m$xlevels[[v]])) stop("Unseen model factor level in ", v)
    data.table::set(d, j = v, value = factor(vals, levels = m$xlevels[[v]]))
  }
  d
}
design_matrix <- function(d, m) {
  X <- Matrix::sparse.model.matrix(delete.response(terms(m$formula)), data = as.data.frame(d),
     contrasts.arg = m$contrasts)
  if (!setequal(colnames(X), names(m$beta))) stop("Prediction design columns do not match fitted model.")
  X[, names(m$beta), drop = FALSE]
}
rank_audit <- function(d, formula, prefix) {
  # Exact cross-product over ALL records, chunked to avoid a huge dense X.
  ff <- delete.response(terms(formula))
  G <- NULL
  for (start in seq.int(1L, nrow(d), by = 50000L)) {
    z <- as.data.frame(d[start:min(nrow(d), start + 49999L)])
    X <- Matrix::sparse.model.matrix(ff, data = z)
    if (is.null(G)) G <- matrix(0, ncol(X), ncol(X), dimnames = list(colnames(X), colnames(X)))
    G <- G + as.matrix(Matrix::crossprod(X))
  }
  den <- sqrt(diag(G)); den[den == 0] <- 1
  R <- G / outer(den, den)
  e <- eigen(R, symmetric = TRUE)
  tol <- 1e-10
  rank <- sum(e$values > tol * max(e$values))
  condition <- if (min(e$values) > 0) sqrt(max(e$values)/min(e$values)) else Inf
  write_csv(data.frame(columns = ncol(G), rank = rank, condition_number_scaled = condition,
    min_eigenvalue_scaled = min(e$values), tolerance = tol), paste0(prefix, "_design_rank.csv"))
  save_rds(list(crossproduct = G, scaled_crossproduct = R, eigenvalues = e$values), paste0(prefix, "_design_matrix_audit.rds"))
  if (rank < ncol(G)) {
    bad <- which(e$values <= tol * max(e$values))
    loads <- data.table::rbindlist(lapply(bad, function(j) data.frame(direction = j, term = rownames(G), loading = e$vectors[, j])))
    write_csv(loads, paste0(prefix, "_dependent_directions.csv"))
    stop("Rank-deficient design. See saved audit; no automatic dropping or changing covariates.")
  }
  invisible(condition)
}
figure_save <- function(p, name, width = 9, height = 6) {
  # Figures intentionally have no title, subtitle, caption, or figure number.
  custom_theme <- do.call(ggplot2::theme, p$theme)
  p <- p + ggplot2::labs(title = NULL, subtitle = NULL, caption = NULL) +
    ggplot2::theme_classic(base_size = 12) + custom_theme + ggplot2::theme(axis.text = ggplot2::element_text(colour = "black"))
  ggplot2::ggsave(file.path(DIRS$figures, paste0(name, ".pdf")), p, width = width, height = height)
  ggplot2::ggsave(file.path(DIRS$figures, paste0(name, ".tiff")), p, width = width, height = height,
    dpi = 600, compression = "lzw", device = "tiff")
  ggplot2::ggsave(file.path(DIRS$figures, paste0(name, ".png")), p, width = width, height = height, dpi = 200)
}
