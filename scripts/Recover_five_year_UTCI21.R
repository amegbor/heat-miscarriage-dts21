# Targeted recovery: reuse completed nlminb models, fit raw GM2/3 with BFGS.
bundle <- Sys.getenv('DTS_BUNDLE',file.path(getwd(),'scripts'))
Sys.setenv(DTS_HEAT='UTCI', DTS_OPTIMIZER='nlminb')
source(file.path(bundle,'supp','00_settings.R'))
source(file.path(bundle,'supp','01_helpers.R'))
stopifnot(HEAT_RUN=='UTCI')
# These source nlminb checkpoints must exist and pass their original checks.
reuse <- c(paste0('Five_year_reported_UTCI_IQR_GM',1:3),
           'Five_year_reported_UTCI_RAW_GM1')
stopifnot(all(file.exists(file.path(OUT,'models',paste0(reuse,'.rds')))))
fit_saved_original <- fit_saved
fit_saved <- function(f,d,id) {
  targets <- paste0('Five_year_reported_UTCI_RAW_GM',2:3)
  if(!id %in% targets) return(fit_saved_original(f,d,id))
  previous <- OPTIMIZER
  OPTIMIZER <<- 'BFGS'
  on.exit({OPTIMIZER <<- previous},add=TRUE)
  rescue_id <- paste0(id,'_BFGS_recovery21')
  message('Targeted recovery: ',rescue_id)
  # Separate cache IDs preserve the original failed fit and its metadata.
  # Existing BFGS recovery checkpoints are validated/reused on a repeat run.
  result <- fit_saved_original(f,d,rescue_id)
  if(id==targets[1]) {
    old <- fread(file.path(OUT,'tables',paste0(id,'_fixed_effects.csv')))
    b <- fixef(result)$cond
    se <- sqrt(diag(vcov(result)$cond))
    comparison <- data.frame(term=names(b),estimate_before=old$estimate[match(names(b),old$term)],
      estimate_BFGS=unname(b),SE_before=old$SE[match(names(b),old$term)],SE_BFGS=unname(se))
    comparison$change_in_original_SE_units <-
      (comparison$estimate_BFGS-comparison$estimate_before)/comparison$SE_before
    write_tab(comparison,'Five_year_UTCI_RAW_GM2_optimizer_comparison','qc')
  }
  result
}
# Same stored five-year flag, sample rules, IQR denominators and pooled years.
# Combined six-model output tables keep their expected filenames.
source(file.path(bundle,'supp','05_five_year.R'))
write_tab(data.frame(
  analysis=paste0('Five_year_reported_UTCI_RAW_GM',2:3),
  checkpoint=paste0('Five_year_reported_UTCI_RAW_GM',2:3,'_BFGS_recovery21.rds'),
  optimizer='BFGS',note='Targeted recovery; original checkpoints preserved'),
  'Five_year_UTCI_recovery21_provenance','qc')
writeLines(paste('UTCI five-year recovery and combined tables completed:',Sys.time()),
  file.path(OUT,'qc','FIVE_YEAR_UTCI_RECOVERY_OK.txt'))
message('UTCI FIVE-YEAR RECOVERY COMPLETE')
