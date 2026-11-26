HERE <- file.path(Sys.getenv("DTS_BUNDLE"),"pif")
source(file.path(HERE, "00_config.R")); source(file.path(HERE, "01_helpers.R"))
need(c("data.table", "qgcomp", "sandwich", "glmmTMB", "Matrix"))
data.table::setDTthreads(1L)
if (!file.exists(file.path(DIRS$prepared, "PREPARED_OK.rds"))) stop("Run 02_prepare.R first.")
mix <- Sys.getenv("DTS_MIXTURE", "Heat_only")
engine <- Sys.getenv("DTS_MIXTURE_ENGINE", "package")
if (!mix %in% names(MIXTURES)) stop("Unknown mixture: ", mix)
if (!engine %in% c("package", "mixed")) stop("Engine must be package or mixed.")
optimizer <- Sys.getenv("DTS_MIXTURE_OPTIMIZER", "nlminb")
if (!optimizer %in% c("nlminb", "BFGS")) stop("Choose nlminb or BFGS explicitly; there is no fallback loop.")
label <- paste(engine, mix, if (engine == "mixed") optimizer else "cloglog", sep = "_")
prefix <- file.path(DIRS$tables, label)
save_session(label)
model_file <- file.path(DIRS$models, paste0(label, ".rds"))
if (file.exists(model_file)) stop("Checkpoint exists: ", model_file, ". Inspect its diagnostics rather than overwriting/refitting automatically.")
dat <- data.table::as.data.table(readRDS(file.path(DIRS$prepared, "Mixture_common_quartiles.rds")))
vars <- MIXTURES[[mix]]
ff <- fixed_formula(vars)
rank_audit(dat, ff, file.path(DIRS$qc, label))
writeLines(capture.output(ff), paste0(prefix, "_fixed_formula.txt"))
warn <- character()
started <- Sys.time()
if (engine == "package") {
  # Native qgcomp noboot uses an additive GLM. It ignores id for robust SEs,
  # so the inferential covariance is explicitly replaced by DHS-cluster HC0.
  # Quantization was performed ONCE on the common sample in preparation.
  qfit <- withCallingHandlers(qgcomp::qgcomp.glm.noboot(f = ff, data = as.data.frame(dat),
    expnms = vars, q = NULL, family = binomial(link = "cloglog"), bayes = FALSE,
    control = glm.control(maxit = 100L, epsilon = 1e-8)),
    warning = function(w) {warn <<- c(warn, conditionMessage(w)); invokeRestart("muffleWarning")})
  fit <- qfit$fit
  save_rds(qfit, model_file)
  writeLines(if (length(warn)) unique(warn) else "No fit warnings captured.", paste0(prefix, "_warnings.txt"))
  if (family(fit)$link != "cloglog") stop("Package changed the link; do not report this fit as a hazard model.")
  b <- coef(fit)
  V <- sandwich::vcovCL(fit, cluster = dat$cluster_uid, type = "HC0", cadjust = TRUE)
  se <- sqrt(diag(V))
  dg <- data.frame(model = label, engine = engine, converged = fit$converged,
    pdHess = NA, all_fixed_SE_finite = all(is.finite(se)),
    all_fixed_estimates_finite = all(is.finite(b)), rank = fit$rank, columns = length(b),
    logLik = as.numeric(logLik(fit)), AIC = AIC(fit), BIC = BIC(fit),
    n_parameters = attr(logLik(fit), "df"), N = nobs(fit), iterations = fit$iter,
    family = family(fit)$family, link = family(fit)$link,
    dependence = "DHS-cluster HC0 sandwich; mother random intercept absent")
  package_psi <- as.numeric(qfit$psi)[1]
  write_csv(data.frame(package_psi = package_psi, coefficient_sum = sum(b[vars]),
    difference = package_psi - sum(b[vars]), package_model_based_varpsi = as.numeric(qfit$var.psi)[1]),
    paste0(prefix, "_package_contrast_check.csv"))
  write_csv(dg, paste0(prefix, "_diagnostics.csv"))
  if (!fit$converged || fit$rank != length(b) || !all(is.finite(c(b,V,se))) ||
     !all(is.finite(c(dg$logLik, dg$AIC))) || abs(package_psi - sum(b[vars])) > 1e-6) stop("Invalid qgcomp GLM or contrast. Checkpoint retained.")
  phat <- fitted(fit)
  save_rds(list(covariance = V, cluster_unit = "cluster_uid", type = "HC0", cadjust = TRUE),
    paste0(prefix, "_cluster_robust_covariance.rds"))
  rm(qfit)
} else {
  # Explicit mixed-model extension of additive quantile g-computation.
  # This is not a native qgcomp package GLMM implementation.
  f <- as.formula(paste(paste(deparse(ff), collapse = " "), "+ (1 | mother_uid)"))
  control <- if (optimizer == "nlminb") glmmTMB::glmmTMBControl(optCtrl = list(iter.max = 10000L, eval.max = 20000L), parallel = 1L) else
    glmmTMB::glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS"), optCtrl = list(maxit = 20000L), parallel = 1L)
  fit <- withCallingHandlers(glmmTMB::glmmTMB(f, data = as.data.frame(dat), family = binomial(link = "cloglog"), control = control),
    warning = function(w) {warn <<- c(warn, conditionMessage(w)); invokeRestart("muffleWarning")})
  save_rds(fit, model_file)
  writeLines(if (length(warn)) unique(warn) else "No fit warnings captured.", paste0(prefix, "_warnings.txt"))
  dg <- fit_diag(fit, label); dg$engine <- engine
  dg$dependence <- "Mother random intercept; GLMM model-based covariance"
  write_csv(dg, paste0(prefix, "_diagnostics.csv"))
  writeLines(capture.output(glmmTMB::VarCorr(fit)), paste0(prefix, "_VarCorr.txt"))
  diag_text <- tryCatch(capture.output(glmmTMB::diagnose(fit)), error = function(e) conditionMessage(e))
  writeLines(diag_text, paste0(prefix, "_diagnose.txt"))
  if (!valid_diag(dg)) stop("Mixed model failed diagnostic gate; checkpoint retained. No automatic optimizer retry.")
  b <- glmmTMB::fixef(fit)$cond; V <- vcov(fit)$cond
  phat <- predict(fit, type = "response")
}
writeLines(if (length(warn)) unique(warn) else "No fit warnings captured.", paste0(prefix, "_warnings.txt"))
save_fixed(b, V, prefix)
if (!all(vars %in% names(b))) stop("Missing mixture coefficients.")
psi <- sum(b[vars]); varpsi <- sum(V[vars, vars, drop = FALSE])
if (!is.finite(varpsi) || varpsi <= 0) stop("Invalid joint mixture variance.")
sepsi <- sqrt(varpsi)
result <- data.frame(model = label, mixture = mix, mixture_label = MIX_LABELS[[mix]], engine = engine,
 components = paste(vars, collapse = ";"), q = 4, psi = psi, SE = sepsi,
 HR = exp(psi), lower_95 = exp(psi - qnorm(.975) * sepsi), upper_95 = exp(psi + qnorm(.975) * sepsi),
 p_value = 2 * pnorm(abs(psi/sepsi), lower.tail = FALSE), person_months = nrow(dat),
 pregnancies = data.table::uniqueN(dat$preg_uid), mothers = data.table::uniqueN(dat$mother_uid),
 clusters = data.table::uniqueN(dat$cluster_uid), events = sum(dat$miscarriage_event),
 logLik = as.numeric(logLik(fit)), AIC = AIC(fit), BIC = BIC(fit),
 elapsed_minutes = as.numeric(difftime(Sys.time(), started, units = "mins")),
 estimand = "Conditional hazard ratio per simultaneous one-quartile increase; additive cloglog model",
 uncertainty = if (engine == "package") "DHS-cluster robust HC0 sandwich" else "GLMM model-based delta method")
weights <- data.frame(model = label, mixture = mix, engine = engine, component = vars, coefficient = unname(b[vars]),
 direction = ifelse(b[vars] > 0, "Positive", ifelse(b[vars] < 0, "Negative", "Zero")), weight = 0)
for (dr in c("Positive", "Negative")) {
 ix <- weights$direction == dr
 if (any(ix)) weights$weight[ix] <- abs(weights$coefficient[ix]) / sum(abs(weights$coefficient[ix]))
}
write_csv(weights, paste0(prefix, "_weights.csv"))
write_csv(result,paste0(prefix,"_results.csv"))
message('Joint heat-air mixture complete: ',label)
