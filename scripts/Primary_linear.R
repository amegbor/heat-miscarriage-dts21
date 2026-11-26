# ============================================================================
# 01A_DTS_Linear_Batch_REVISED_v2.R
# Revised batch-safe DTS linear models
#
# Key changes after diagnostic WBGT GM1:
#   * country = fixed effect
#   * cluster_uid random intercept removed (variance ~0)
#   * country random intercept removed (variance ~0)
#   * mother_uid random intercept retained
#   * GM2/GM3 retain pregnancy-start year/month adjustment
#   * country x pregnancy-start-month interaction removed
#   * no automatic BFGS retry
#   * invalid fits stop the batch job immediately
# ============================================================================

# Fresh-process runner; preserve caller variables.
invisible(gc(reset = TRUE))

suppressPackageStartupMessages({
  library(glmmTMB)
  library(dplyr)
  library(tibble)
  library(broom.mixed)
  library(readr)
})

Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1"
)

BASE_DIR <- getwd()
DTS_DIR <- Sys.getenv("DTS_ROOT",file.path(BASE_DIR,"output"))

DATA_FILE <- file.path(DTS_DIR, "data", "DTS_ObservedRiskSet.rds")
MODEL_DIR <- file.path(DTS_DIR, "models", "Linear_Nonlinear")
TABLE_DIR <- file.path(DTS_DIR, "tables", "Linear_Nonlinear")
LOG_DIR   <- file.path(DTS_DIR, "logs")

for (d in c(MODEL_DIR, TABLE_DIR, LOG_DIR)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

if (!file.exists(DATA_FILE)) {
  stop("Observed DTS risk-set file not found: ", DATA_FILE)
}

HEAT_TO_RUN <- Sys.getenv("DTS_HEAT", unset = "")
allowed_heat <- c("WBGT", "UTCI", "MRT", "Heat Index")

if (!nzchar(HEAT_TO_RUN)) {
  stop("DTS_HEAT was not supplied.")
}
if (!HEAT_TO_RUN %in% allowed_heat) {
  stop(
    "Invalid DTS_HEAT='", HEAT_TO_RUN,
    "'. Expected one of: ", paste(allowed_heat, collapse = ", ")
  )
}

message("Loading: ", DATA_FILE)
DT <- readRDS(DATA_FILE)

message("Loaded person-month records: ", format(nrow(DT), big.mark = ","))
message("Pregnancies: ", format(dplyr::n_distinct(DT$preg_uid), big.mark = ","))
message("Events: ", format(sum(DT$miscarriage_event, na.rm = TRUE), big.mark = ","))

heat_lookup <- tibble::tribble(
  ~heat_label,   ~scaled_var,                ~raw_var,
  "WBGT",        "WBGT_T1_country_iqr",      "WBGT_T1_cummean",
  "UTCI",        "utci_T1_country_iqr",      "utci_T1_cummean",
  "MRT",         "mrt_T1_country_iqr",       "mrt_T1_cummean",
  "Heat Index",  "HI_T1_country_iqr",        "HI_T1_cummean"
)

heat_row <- heat_lookup %>%
  dplyr::filter(heat_label == HEAT_TO_RUN)

if (nrow(heat_row) != 1L) {
  stop("Could not uniquely identify heat metric: ", HEAT_TO_RUN)
}

heat_label <- heat_row$heat_label[[1]]
heat_var   <- heat_row$scaled_var[[1]]
metric_name_file <- gsub(" ", "_", heat_label)

# Final manuscript covariates; cooking fuel and marital status intentionally excluded.
model3_covariates <- c(
  "age_at_preg_start",
  "I(age_at_preg_start^2)",
  "rc_empl",
  "v106",
  "v190",
  "rc_media_none",
  "rc_hins_any",
  "v025"
)

base_hazard_term <- "gest_month_f"
temporal_terms <- c("preg_start_year_f", "preg_start_month_f")
country_fixed_term <- "country"
mother_random_term <- "(1 | mother_uid)"

make_formula <- function(heat_term, model_number = 1L) {
  model_number <- as.integer(model_number)

  if (!model_number %in% 1:3) {
    stop("model_number must be 1, 2, or 3.")
  }

  fixed_terms <- c(
    base_hazard_term,
    heat_term,
    country_fixed_term
  )

  if (model_number >= 2L) {
    fixed_terms <- c(
      fixed_terms,
      temporal_terms
    )
  }

  if (model_number >= 3L) {
    fixed_terms <- c(
      fixed_terms,
      model3_covariates
    )
  }

  stats::as.formula(
    paste(
      "miscarriage_event ~",
      paste(c(fixed_terms, mother_random_term), collapse = " + ")
    )
  )
}

get_model_data <- function(formula, data) {
  formula_vars <- unique(all.vars(formula))
  missing_vars <- setdiff(formula_vars, names(data))

  if (length(missing_vars) > 0L) {
    stop("Missing model variables: ", paste(missing_vars, collapse = ", "))
  }

  keep_vars <- unique(c(formula_vars, "preg_uid"))

  complete_rows <- stats::complete.cases(
    data[, formula_vars, drop = FALSE]
  )

  out <- data[
    complete_rows,
    keep_vars,
    drop = FALSE
  ]

  rm(complete_rows)
  droplevels(out)
}

model_diagnostics <- function(fit) {
  fixed_effects <- tryCatch(
    glmmTMB::fixef(fit)$cond,
    error = function(e) numeric(0)
  )

  vc <- tryCatch(
    stats::vcov(fit)$cond,
    error = function(e) NULL
  )

  se <- if (
    is.null(vc) ||
    length(vc) == 0L ||
    any(!is.finite(vc))
  ) {
    rep(NA_real_, length(fixed_effects))
  } else {
    suppressWarnings(sqrt(diag(vc)))
  }

  gradient <- tryCatch(
    if(length(fit$sdr$gradient.fixed)) fit$sdr$gradient.fixed else fit$fit$gradient,
    error = function(e) NULL
  )

  max_grad <- if (
    is.null(gradient) ||
    length(gradient) == 0L ||
    all(!is.finite(gradient))
  ) {
    NA_real_
  } else {
    max(abs(gradient), na.rm = TRUE)
  }

  tibble::tibble(
    convergence_code = tryCatch(as.integer(fit$fit$convergence),
                                error = function(e) NA_integer_),
    optimizer_message = tryCatch(if(length(fit$fit$message)) as.character(fit$fit$message)[1L] else NA_character_,
                                 error = function(e) NA_character_),
    pdHess = tryCatch(isTRUE(fit$sdr$pdHess),
                      error = function(e) FALSE),
    all_fixed_SE_finite = (
      length(se) == length(fixed_effects) &&
      length(se) > 0L &&
      all(is.finite(se))
    ),
    max_abs_gradient = max_grad,
    logLik = tryCatch(as.numeric(stats::logLik(fit)),
                      error = function(e) NA_real_),
    AIC = tryCatch(as.numeric(stats::AIC(fit)),
                   error = function(e) NA_real_),
    N = tryCatch(stats::nobs(fit),
                 error = function(e) NA_integer_)
  )
}

fit_is_inferentially_valid <- function(diagnostics) {
  isTRUE(diagnostics$convergence_code[[1]] == 0L) &&
    isTRUE(diagnostics$pdHess[[1]]) &&
    isTRUE(diagnostics$all_fixed_SE_finite[[1]])
}

fit_glmmtmb_safe <- function(formula, data, checkpoint_file) {
  diagnostic_file <- sub(
    "\\.rds$",
    "_optimizer_diagnostics.csv",
    checkpoint_file
  )

  if (file.exists(checkpoint_file)) {
    message("Loading existing checkpoint: ", basename(checkpoint_file))
    fit_existing <- readRDS(checkpoint_file)
    dg_existing <- model_diagnostics(fit_existing)
    print(dg_existing)

    if (!fit_is_inferentially_valid(dg_existing)) {
      stop(
        "Existing checkpoint is not inferentially valid: ",
        checkpoint_file
      )
    }

    return(fit_existing)
  }

  model_data <- get_model_data(
    formula = formula,
    data = data
  )

  message(
    "Fitting N person-months = ",
    format(nrow(model_data), big.mark = ","),
    "; pregnancies = ",
    format(dplyr::n_distinct(model_data$preg_uid), big.mark = ","),
    "; events = ",
    format(sum(model_data$miscarriage_event, na.rm = TRUE), big.mark = ",")
  )

  message(if(heat_label == "Heat Index" && m == 3L) "Fitting with BFGS optimizer..." else "Fitting with nlminb optimizer...")

  fit <- suppressWarnings(
    glmmTMB::glmmTMB(
      formula = formula,
      data = model_data,
      family = stats::binomial(link = "cloglog"),
      control = if(heat_label == "Heat Index" && m == 3L) glmmTMB::glmmTMBControl(optimizer=optim,optArgs=list(method="BFGS"),optCtrl=list(maxit=20000L)) else glmmTMB::glmmTMBControl(optCtrl=list(iter.max=10000L,eval.max=20000L))
    )
  )

  rm(model_data)
  invisible(gc(full = TRUE))

  dg <- model_diagnostics(fit)

  message("")
  message("Model diagnostics:")
  print(dg)

  readr::write_csv(dg, diagnostic_file)

  saveRDS(
    fit,
    checkpoint_file,
    compress = FALSE
  )

  message("Model checkpoint saved: ", checkpoint_file)
  message("Optimizer diagnostics saved: ", diagnostic_file)

  if (!fit_is_inferentially_valid(dg)) {

    stop(
      paste0(
        "\nModel completed but failed inferential diagnostics.\n",
        "Required: convergence_code=0, pdHess=TRUE, ",
        "all_fixed_SE_finite=TRUE.\n",
        "No subsequent model will be fitted in this batch job.\n",
        "Checkpoint retained: ", checkpoint_file
      )
    )
  }

  fit
}

extract_linear_heat <- function(
    fit,
    heat_var,
    heat_label,
    model_label,
    scale_label
) {
  tt <- broom.mixed::tidy(
    fit,
    effects = "fixed"
  ) %>%
    dplyr::filter(term == heat_var)

  if (nrow(tt) != 1L) {
    stop("Could not uniquely extract heat coefficient: ", heat_var)
  }

  if (
    !is.finite(tt$estimate[[1]]) ||
    !is.finite(tt$std.error[[1]])
  ) {
    stop("Heat coefficient or SE is non-finite for: ", heat_var)
  }

  tibble::tibble(
    Model = model_label,
    Heat_metric = heat_label,
    Exposure_scale = scale_label,
    estimate_log_HR = tt$estimate,
    SE = tt$std.error,
    HR = exp(tt$estimate),
    lower_95 = exp(tt$estimate - 1.96 * tt$std.error),
    upper_95 = exp(tt$estimate + 1.96 * tt$std.error),
    p_value = tt$p.value,
    N_person_months = stats::nobs(fit)
  )
}

message("")
message("============================================================")
message("DTS LINEAR BATCH ANALYSIS - REVISED v2")
message("Heat metric: ", heat_label)
message("Exposure variable: ", heat_var)
message("Country: fixed effect")
message("Cluster random intercept: removed")
message("Mother random intercept: retained")
message("Country x pregnancy-start-month interaction: excluded")
message("Start time: ", Sys.time())
message("============================================================")

metric_results <- list()
metric_diagnostics <- list()

for (m in 1:3) {
  if (m == 1L) {
    # Reuse the valid GM1 checkpoint from the previous run.
    model_name <- paste0(
      "GM1_", metric_name_file,
      "_DTS_IQR_CountryFE_MotherRE"
    )
  } else {
    # New names prevent the failed prior GM2 from being loaded.
    model_name <- paste0(
      "GM", m, "_", metric_name_file,
      "_DTS_IQR_CountryFE_MotherRE_NoCountryMonthInt"
    )
  }

  if(heat_label == "Heat Index" && m == 3L) model_name <- "GM3_Heat_Index_DTS_IQR_CountryFE_MotherRE_BFGS"
  checkpoint_file <- file.path(
    MODEL_DIR,
    paste0(model_name, ".rds")
  )

  message("")
  message("============================================================")
  message("RUNNING: ", model_name)
  message("MODEL ", m, " OF 3 FOR ", heat_label)
  message("============================================================")

  fml <- make_formula(
    heat_term = heat_var,
    model_number = m
  )

  message("Formula:")
  print(fml)

  fit <- fit_glmmtmb_safe(
    formula = fml,
    data = DT,
    checkpoint_file = checkpoint_file
  )

  result_one <- extract_linear_heat(
    fit = fit,
    heat_var = heat_var,
    heat_label = heat_label,
    model_label = paste0("Model ", m),
    scale_label = "Country-specific IQR"
  )

  diagnostic_one <- model_diagnostics(fit) %>%
    dplyr::mutate(
      Model = paste0("Model ", m),
      Heat_metric = heat_label,
      Exposure_scale = "Country-specific IQR",
      Model_name = model_name,
      Country_effect = "Fixed",
      Cluster_random_effect = "Excluded",
      Mother_random_effect = "Included",
      Country_month_interaction = "Excluded"
    )

  print(result_one)
  print(diagnostic_one)

  result_file <- file.path(
    TABLE_DIR,
    paste0(model_name, "_HR.csv")
  )

  diagnostic_result_file <- file.path(
    TABLE_DIR,
    paste0(model_name, "_diagnostics.csv")
  )

  readr::write_csv(result_one, result_file)
  readr::write_csv(diagnostic_one, diagnostic_result_file)

  metric_results[[model_name]] <- result_one
  metric_diagnostics[[model_name]] <- diagnostic_one

  message("Completed and saved: ", model_name)

  rm(
    fit,
    result_one,
    diagnostic_one,
    fml
  )

  invisible(gc(full = TRUE))

  message("Memory cleanup completed before next model.")
}

metric_results_df <- dplyr::bind_rows(metric_results)
metric_diagnostics_df <- dplyr::bind_rows(metric_diagnostics)

combined_results_file <- file.path(
  TABLE_DIR,
  paste0(
    "DTS_Linear_",
    metric_name_file,
    "_GM1_GM2_GM3_CountryFE_MotherRE_v2.csv"
  )
)

combined_diagnostics_file <- file.path(
  TABLE_DIR,
  paste0(
    "DTS_Linear_",
    metric_name_file,
    "_Diagnostics_CountryFE_MotherRE_v2.csv"
  )
)

readr::write_csv(metric_results_df, combined_results_file)
readr::write_csv(metric_diagnostics_df, combined_diagnostics_file)

session_file <- file.path(
  TABLE_DIR,
  paste0(
    "SessionInfo_01A_",
    metric_name_file,
    "_CountryFE_MotherRE_v2.txt"
  )
)

writeLines(
  capture.output(sessionInfo()),
  session_file
)

message("")
message("============================================================")
message(heat_label, " LINEAR DTS ANALYSIS COMPLETED")
message("End time: ", Sys.time())
message("Combined HR results: ", combined_results_file)
message("Combined diagnostics: ", combined_diagnostics_file)
message("Model checkpoints: ", MODEL_DIR)
message("============================================================")

print(metric_results_df)
print(metric_diagnostics_df)

invisible(gc(full = TRUE))
