HERE <- file.path(Sys.getenv("DTS_BUNDLE"),"pif")
source(file.path(HERE, "00_config.R")); source(file.path(HERE, "01_helpers.R"))
need(c("data.table", "glmmTMB", "lme4", "Matrix", "qgcomp", "sandwich", "statmod", "ggplot2"))
data.table::setDTthreads(1L)
save_session("prepare")
marker <- file.path(DIRS$prepared, "PREPARED_OK.rds")
if (file.exists(marker)) stop("Preparation already exists. Use a new DTS_ROOT/output directory for a different analysis.")

# Validate and extract each large primary checkpoint sequentially.
audit <- list(); compact <- list()
for (metric in HEAT$metric) {
  f <- MODEL_FILES[[metric]]
  if (!file.exists(f)) stop("Required validated GM3 checkpoint missing: ", f,
    "\nComplete and validate the primary models before preparation.")
  message("Checking GM3: ", metric)
  fit <- readRDS(f)
  if (!inherits(fit, "glmmTMB")) stop("Checkpoint is not a glmmTMB object: ", f)
  dg <- fit_diag(fit, metric)
  dg$source_file <- f; dg$size_bytes <- file.info(f)$size
  dg$source_mtime <- as.character(file.info(f)$mtime)
  audit[[metric]] <- dg
  write_csv(data.table::rbindlist(audit), file.path(DIRS$qc, "Primary_GM3_validation.csv"))
  writeLines(capture.output(glmmTMB::VarCorr(fit)), file.path(DIRS$qc, paste0(HEAT$stem[match(metric, HEAT$metric)], "_primary_VarCorr.txt")))
  if (!valid_diag(dg)) stop("GM3 failed validation: ", metric, ". Invalid results will not be used.")
  compact[[metric]] <- compact_model(fit, metric)
  save_fixed(compact[[metric]]$beta, vcov(fit)$cond,
    file.path(DIRS$tables, paste0(HEAT$stem[match(metric, HEAT$metric)], "_primary_GM3")))
  rm(fit); invisible(gc(full = TRUE))
}
horizon <- compact[[1]]$horizon
if (!all(vapply(compact, function(m) identical(m$horizon, horizon), logical(1)))) stop("GM3 models have different risk horizons.")
if (!(identical(horizon, 1:6) || identical(horizon, 0:6))) stop("Unexpected month indexing; inspect data preparation.")
save_rds(compact, file.path(DIRS$prepared, "GM3_compact_parameters.rds"))

message("Loading observed risk set")
obs <- data.table::as.data.table(readRDS(OBS_FILE))
assert_cols(obs, c("miscarriage_event", "preg_uid", "mother_uid", "cluster_uid", COV_VARS, HEAT$scaled, EXPOSURES))
if (anyNA(obs$preg_uid) || anyNA(obs$mother_uid) || anyNA(obs$cluster_uid)) stop("Missing linkage IDs.")
check_risk(obs)
primary_counts <- data.table::rbindlist(lapply(seq_len(nrow(HEAT)), function(j) {
  v <- c(COV_VARS, HEAT$scaled[j], "mother_uid", "miscarriage_event")
  okj <- complete_finite(obs, v)
  z <- counts(obs[okj], HEAT$metric[j])
  z$checkpoint_person_months <- compact[[HEAT$metric[j]]]$diagnostic$N
  z
}))
write_csv(primary_counts, file.path(DIRS$qc, "Primary_GM3_exact_sample_counts.csv"))
if (any(primary_counts$person_months != primary_counts$checkpoint_person_months)) stop("Observed data complete-case counts differ from the fitted GM3 models.")
flow <- list(counts(obs, "Observed primary risk set"))
needed <- c("miscarriage_event", "preg_uid", "mother_uid", "cluster_uid", COV_VARS, EXPOSURES)
write_csv(data.frame(variable = needed, missing = vapply(obs[, ..needed], function(x) sum(is.na(x)), integer(1)),
  percent_missing = 100 * vapply(obs[, ..needed], function(x) mean(is.na(x)), numeric(1))), file.path(DIRS$qc, "Mixture_missingness.csv"))
ok <- complete_finite(obs, needed)
qdat <- drop_factors(data.table::copy(obs[ok, ..needed]))
factor_audit <- data.table::rbindlist(lapply(COV_VARS, function(v) {
  if (!is.factor(obs[[v]])) return(NULL)
  data.frame(variable = v, stored_levels = nlevels(obs[[v]]),
    observed_levels_before_complete_cases = data.table::uniqueN(obs[[v]][!is.na(obs[[v]])]),
    retained_levels = nlevels(qdat[[v]]), removed_levels = nlevels(obs[[v]]) - nlevels(qdat[[v]]),
    retained_reference_level = levels(qdat[[v]])[1])
}))
write_csv(factor_audit, file.path(DIRS$qc, "Mixture_factor_level_audit.csv"))
flow[[2]] <- counts(qdat, "Common complete-case sample for ALL mixture specifications")
qc <- check_risk(qdat)
write_csv(data.frame(pregnancies_with_internal_missing_intervals = sum(qc$gaps)), file.path(DIRS$qc, "Mixture_followup_gaps.csv"))
write_csv(qdat[, .(person_months = .N, pregnancies = data.table::uniqueN(preg_uid),
 mothers = data.table::uniqueN(mother_uid), clusters = data.table::uniqueN(cluster_uid), events = sum(miscarriage_event)), by = country],
 file.path(DIRS$qc, "Mixture_country_counts.csv"))
# Cluster robust package reference requires a mother to belong to only one DHS cluster.
nest <- qdat[, .(clusters = data.table::uniqueN(cluster_uid)), by = mother_uid]
if (any(nest$clusters != 1L)) stop("Mothers cross DHS clusters; single-cluster sandwich needs reconsideration.")

# Cutpoints are defined on common complete-case PERSON-MONTH exposures, pooled
# across countries, exactly as qgcomp's default quantization. No silent jitter.
breaks0 <- lapply(EXPOSURES, function(v) quantile(qdat[[v]], probs = seq(0, 1, .25), type = 7, names = FALSE))
names(breaks0) <- EXPOSURES
tiecheck <- vapply(breaks0, function(z) anyDuplicated(z) > 0, logical(1))
write_csv(data.table::rbindlist(lapply(EXPOSURES, function(v) data.frame(variable = v,
 probability = seq(0, 1, .25), value = breaks0[[v]], tied_cutpoints = tiecheck[[v]]))), file.path(DIRS$qc, "Mixture_quartile_cutpoints.csv"))
if (any(tiecheck)) stop("Tied quartile cutpoints. Inspect saved table; do not randomly break ties.")
quant <- qgcomp::quantize(as.data.frame(qdat), expnms = EXPOSURES, q = 4)
qraw <- data.table::copy(qdat[, ..EXPOSURES])
qdat <- data.table::as.data.table(quant$data)
save_rds(quant$breaks, file.path(DIRS$prepared, "Mixture_quantize_breaks.rds"))
rm(quant)
cor0 <- cor(qraw, method = "spearman")
corq <- cor(qdat[, ..EXPOSURES], method = "spearman")
write_csv(data.frame(component = rownames(cor0), cor0, check.names = FALSE), file.path(DIRS$qc, "Mixture_Spearman_raw.csv"))
write_csv(data.frame(component = rownames(corq), corq, check.names = FALSE), file.path(DIRS$qc, "Mixture_Spearman_quantized.csv"))
rm(qraw)
write_csv(data.table::rbindlist(lapply(EXPOSURES, function(v) qdat[, .(component = v, category = get(v), events = miscarriage_event)][,
 .(person_months = .N, events = sum(events)), by = .(component, category)])), file.path(DIRS$qc, "Mixture_category_counts.csv"))
for (mix in names(MIXTURES)) {
  # Every model is audited on the same sample before submission.
  message("Rank audit: ", mix)
  rank_audit(qdat, fixed_formula(MIXTURES[[mix]]), file.path(DIRS$qc, mix))
  # Positivity/support summaries for the hypothetical equal-quartile contrast.
  vars <- MIXTURES[[mix]]
  mat <- as.matrix(qdat[, ..vars])
  write_csv(data.frame(joint_category = 0:3,
   n_all_components_equal = vapply(0:3, function(k) sum(rowSums(mat == k) == length(vars)), numeric(1)),
   n_can_increase_all_by_one = sum(rowSums(mat < 3) == length(vars))),
   file.path(DIRS$qc, paste0(mix, "_joint_support.csv")))
}
invisible(drop_factors(qdat))
save_rds(qdat, file.path(DIRS$prepared, "Mixture_common_quartiles.rds"))
save_rds(flow, file.path(DIRS$prepared, "Mixture_sample_flow.rds"))
write_csv(data.table::rbindlist(flow), file.path(DIRS$qc, "Mixture_sample_flow.csv"))

# PIF target must also be eligible for the primary single-exposure models.
primary_ok <- complete_finite(obs, c(COV_VARS, HEAT$scaled, "mother_uid", "miscarriage_event"))
primary_ids <- unique(as.character(obs$preg_uid[primary_ok]))
obs_check <- obs[primary_ok, c("preg_uid", "gest_month_f", HEAT$scaled), with = FALSE]
obs_check[, gest_month_f := as.character(gest_month_f)]
rm(obs, qdat); invisible(gc(full = TRUE))
message("Loading existing full prediction histories")
pred <- data.table::as.data.table(readRDS(PRED_FILE))
required <- unique(c("preg_uid", "mother_uid", COV_VARS, HEAT$scaled, HEAT$raw, HEAT$monthly, HEAT$iqr))
assert_cols(pred, required)
pred <- pred[, ..required]
pred[, gest_month_f := as.character(gest_month_f)]
check_risk(pred, full = TRUE, horizon = horizon)
pflow <- list(counts(pred, "Full binary-comparator prediction cohort"))
pred <- pred[as.character(preg_uid) %in% primary_ids]
pflow[[2]] <- counts(pred, "Eligible in primary GM3 observed risk set")
ok <- complete_finite(pred, required)
# Exclude whole trajectories if ANY horizon exposure/covariate is missing.
bad_ids <- unique(pred$preg_uid[!ok])
pred <- pred[!preg_uid %in% bad_ids]
pflow[[3]] <- counts(pred, "Complete common full prediction trajectories (all heat metrics)")
if (!nrow(pred)) stop("No complete PIF trajectories.")
check_risk(pred, full = TRUE, horizon = horizon)
data.table::setorder(pred, preg_uid, .g)
for (j in seq_len(nrow(HEAT))) {
  s <- HEAT$scaled[j]; r <- HEAT$raw[j]; iq <- HEAT$iqr[j]
  if (any(pred[[iq]] <= 0)) stop("Invalid country IQR.")
  if (max(abs(pred[[s]] - pred[[r]] / pred[[iq]])) > NUMERIC_TOL) stop("Raw/IQR exposure mismatch: ", s)
}
linked <- merge(obs_check, pred[, c("preg_uid", "gest_month_f", HEAT$scaled), with = FALSE],
 by = c("preg_uid", "gest_month_f"), suffixes = c("_observed", "_prediction"))
comparison <- data.table::rbindlist(lapply(HEAT$scaled, function(v) data.frame(variable = v, matched_rows = nrow(linked),
 max_absolute_difference = max(abs(linked[[paste0(v, "_observed")]] - linked[[paste0(v, "_prediction")]])))))
write_csv(comparison, file.path(DIRS$qc, "PIF_observed_prediction_exposure_agreement.csv"))
if (!nrow(linked) || any(comparison$max_absolute_difference > NUMERIC_TOL)) stop("Prediction exposure histories disagree with observed model inputs.")
save_rds(pred, file.path(DIRS$prepared, "PIF_common_full_trajectories.rds"))
write_csv(data.table::rbindlist(pflow), file.path(DIRS$qc, "PIF_sample_flow.csv"))
write_csv(pred[, .(pregnancies = data.table::uniqueN(preg_uid), mothers = data.table::uniqueN(mother_uid)), by = country], file.path(DIRS$qc, "PIF_country_counts.csv"))
save_rds(list(horizon = horizon, models = MODEL_FILES, prepared_at = Sys.time(),
 data_files = c(observed = OBS_FILE, prediction = PRED_FILE),
 input_info = file.info(c(OBS_FILE, PRED_FILE)),
 target = "Binary comparator; model-based net risk; competing events eliminated",
 common_mixture_sample = TRUE), marker)
message("Preparation passed. Horizon labels: ", paste(horizon, collapse = ", "))
