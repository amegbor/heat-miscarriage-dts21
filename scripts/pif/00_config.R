# Shared configuration. Copy ALL bundle files into DTS_Scripts.
# Environment overrides are optional; defaults match the current Torch workflow.
ROOT <- Sys.getenv("DTS_ROOT", file.path(getwd(),'output'))
OBS_FILE <- file.path(ROOT, "data", "DTS_ObservedRiskSet.rds")
PRED_FILE <- Sys.getenv("DTS_PIF_PRED_FILE", file.path(ROOT, "data", "DTS_FullPredictionRiskSet_BinaryComparator.rds"))
OUT <- file.path(ROOT, "PIF_QGComp_v3")
DIRS <- setNames(as.list(file.path(OUT, c("prepared", "models", "tables", "qc", "figures"))),
                 c("prepared", "models", "tables", "qc", "figures"))
for (d in DIRS) dir.create(d, recursive = TRUE, showWarnings = FALSE)
HEAT <- data.frame(metric = c("WBGT", "UTCI", "MRT", "Heat Index"),
  stem = c("WBGT", "utci", "mrt", "HI"),
  scaled = c("WBGT_T1_country_iqr", "utci_T1_country_iqr", "mrt_T1_country_iqr", "HI_T1_country_iqr"),
  stringsAsFactors = FALSE)
HEAT$raw <- paste0(HEAT$stem, "_T1_cummean")
HEAT$monthly <- paste0(HEAT$stem, "_month")
HEAT$iqr <- paste0(HEAT$stem, "_IQR")
MODEL_FILES <- setNames(file.path(ROOT, "models", "Linear_Nonlinear", c(
  "GM3_WBGT_DTS_IQR_CountryFE_MotherRE_NoCountryMonthInt.rds",
  "GM3_UTCI_DTS_IQR_CountryFE_MotherRE_NoCountryMonthInt.rds",
  "GM3_MRT_DTS_IQR_CountryFE_MotherRE_NoCountryMonthInt.rds",
  "GM3_Heat_Index_DTS_IQR_CountryFE_MotherRE_BFGS.rds")), HEAT$metric)
# No fallback to the invalid HI nlminb checkpoint.
AIR <- paste0(c("CO", "DUST", "PM", "SO", "SU"), "_T1_cummean")
EXPOSURES <- c(HEAT$raw, AIR)
COV_TERMS <- c("gest_month_f", "country", "preg_start_year_f", "preg_start_month_f",
 "age_at_preg_start", "I(age_at_preg_start^2)", "rc_empl", "v106", "v190",
 "rc_media_none", "rc_hins_any", "v025")
COV_VARS <- all.vars(reformulate(COV_TERMS))
MIXTURES <- list(Air_WBGT=c(HEAT$raw[1],AIR),Air_UTCI=c(HEAT$raw[2],AIR),Air_MRT=c(HEAT$raw[3],AIR),Air_HI=c(HEAT$raw[4],AIR))
MIX_LABELS <- c(Air_WBGT='Air + WBGT',Air_UTCI='Air + UTCI',Air_MRT='Air + MRT',Air_HI='Air + Heat Index')
SCENARIOS <- c("natural", "cap_p90", "cap_p75", "cap_p50", "shift_down_1iqr")
SCENARIO_LABELS <- c(natural = "Observed exposure", cap_p90 = "Cap at country 90th percentile",
 cap_p75 = "Cap at country 75th percentile", cap_p50 = "Cap at country median",
 shift_down_1iqr = "Reduce by one country IQR")
NUMERIC_TOL <- 1e-7
