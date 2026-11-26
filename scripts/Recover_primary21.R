# Run in a fresh Slurm R process using the existing Torch container.
suppressPackageStartupMessages(library(glmmTMB))
root <- Sys.getenv('DTS_ROOT', file.path(getwd(),'output'))
bundle <- Sys.getenv('DTS_BUNDLE', file.path(getwd(),'scripts'))
primary <- file.path(bundle, 'Primary_linear.R')
stopifnot(file.exists(primary))

# Diagnostic extraction corrections are already incorporated in Primary_linear.R.
md <- file.path(root, 'models', 'Linear_Nonlinear')
qd <- file.path(root, 'qc', 'Primary_recovery21')
dir.create(qd, recursive = TRUE, showWarnings = FALSE)

check_fit <- function(fit) {
  se <- sqrt(diag(vcov(fit)$cond))
  g <- fit$sdr$gradient.fixed
  V <- fit$sdr$cov.fixed
  scaled <- NA_real_
  if (length(g) && is.matrix(V) && nrow(V) == length(g) &&
      all(is.finite(g)) && all(is.finite(V)) && all(diag(V) > 0)) {
    scaled <- max(abs(drop(V %*% g) / sqrt(diag(V))))
  }
  data.frame(
    convergence_code = if(length(fit$fit$convergence)) fit$fit$convergence[1] else NA_integer_,
    optimizer_message = if(length(fit$fit$message)) fit$fit$message[1] else NA_character_,
    pdHess = isTRUE(fit$sdr$pdHess),
    all_fixed_SE_finite = length(se) > 0 && all(is.finite(se)),
    max_abs_gradient = if(length(g) && all(is.finite(g))) max(abs(g)) else NA_real_,
    approximate_change_SE_units = scaled,
    logLik = as.numeric(logLik(fit)), AIC = AIC(fit), N = nobs(fit)
  )
}
valid <- function(d) {
  isTRUE(d$convergence_code == 0) && isTRUE(d$pdHess) &&
    isTRUE(d$all_fixed_SE_finite) && is.finite(d$logLik)
}
export_metric <- function(metric) {
  Sys.setenv(DTS_HEAT = metric, DTS_ROOT = root)
  message('Regenerating primary tables from checkpoints: ', metric)
  # Original runner refuses invalid checkpoints; valid saved models are reused.
  env <- new.env(parent = globalenv())
  sys.source(primary, envir = env)
  rm(env)
  invisible(gc(full = TRUE))
}

# HI: inspect checkpoint and export without fitting any new model.
hi_id <- 'GM3_Heat_Index_DTS_IQR_CountryFE_MotherRE_BFGS'
hi <- readRDS(file.path(md, paste0(hi_id, '.rds')))
dhi <- check_fit(hi)
write.csv(dhi, file.path(qd, 'HI_GM3_checkpoint_diagnostics.csv'), row.names = FALSE)
print(dhi)
if (!valid(dhi)) stop('HI checkpoint failed direct checks; no refit or replacement attempted.')
rm(hi)
invisible(gc(full = TRUE))
export_metric('Heat Index')

# MRT: warm-start only GM3, preserving its original parameterisation and rows.
mrt_id <- 'GM3_MRT_DTS_IQR_CountryFE_MotherRE_NoCountryMonthInt'
canonical <- file.path(md, paste0(mrt_id, '.rds'))
recovery_path <- file.path(md, paste0(mrt_id, '_BFGS_confirmation.rds'))
original_path <- file.path(md, paste0(mrt_id, '_before_recovery21.rds'))
old <- readRDS(if(file.exists(original_path)) original_path else canonical)
before <- check_fit(old)
before_beta <- fixef(old)$cond
before_se <- sqrt(diag(vcov(old)$cond))
write.csv(before, file.path(qd, 'MRT_GM3_before.csv'), row.names = FALSE)
if (!file.exists(recovery_path)) {
  message('Starting MRT GM3 BFGS confirmation: ', Sys.time())
  recovered <- update(old, data = old$frame,
    start = list(beta = unname(fixef(old)$cond),
                 theta = unname(glmmTMB::getME(old, 'theta'))),
    control = glmmTMBControl(optimizer = optim,
      optArgs = list(method = 'BFGS'),
      optCtrl = list(maxit = 20000L, reltol = 1e-10)))
  saveRDS(recovered, recovery_path, compress = FALSE)
} else {
  message('Reusing saved BFGS confirmation checkpoint.')
  rm(old)
  invisible(gc(full = TRUE))
  recovered <- readRDS(recovery_path)
}
rm(list = intersect('old', ls()))
invisible(gc(full = TRUE))
after <- check_fit(recovered)
write.csv(after, file.path(qd, 'MRT_GM3_after.csv'), row.names = FALSE)
b <- fixef(recovered)$cond
s <- sqrt(diag(vcov(recovered)$cond))
comparison <- data.frame(term = names(b), before = before_beta[names(b)],
  after = b, before_SE = before_se[names(b)], after_SE = s,
  change_in_original_SE_units = (b - before_beta[names(b)]) / before_se[names(b)])
write.csv(comparison, file.path(qd, 'MRT_GM3_coefficient_comparison.csv'), row.names = FALSE)
print(after)
print(comparison[comparison$term == 'mrt_T1_country_iqr', ])
# Conservative operational guard, not a universal convergence threshold.
if (!valid(after) || !is.finite(after$approximate_change_SE_units) ||
    after$approximate_change_SE_units > 0.05 ||
    after$N != before$N || after$logLik < before$logLik - 0.01) {
  stop('Confirmation requires review. Candidate and diagnostics saved; original unchanged.')
}
rm(recovered)
invisible(gc(full = TRUE))
# Preserve original bytes before publishing candidate under downstream's expected name.
if (!file.exists(original_path)) stopifnot(file.rename(canonical, original_path))
if (!file.copy(recovery_path, canonical, overwrite = TRUE)) {
  stop('Could not install candidate; original remains at ', original_path)
}
export_metric('MRT')
writeLines(paste('Recovery and exports completed:', Sys.time()),
           file.path(qd, 'RECOVERY_OK.txt'))
message('RECOVERY COMPLETE. See ', qd)
