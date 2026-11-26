# Replicate the binary sensitivity selection EXACTLY: outcome_last5yr == 1.
# Select pregnancies, then retain their complete available observed DTS follow-up.
# No additional pregnancy-initiation/calendar-interval restriction is applied.
h <- metric_row();d <- load_obs();preg <- as.data.table(readRDS(COHORT))
need(preg,c('preg_uid','outcome_last5yr'))
if(anyDuplicated(preg$preg_uid)||anyNA(preg$preg_uid)) stop('Pregnancy cohort key invalid')
flag <- suppressWarnings(as.numeric(as.character(preg$outcome_last5yr)))
if(any(!is.na(preg$outcome_last5yr)&is.na(flag)) || any(!is.na(flag)&!flag %in% c(0,1))) stop('outcome_last5yr must be 0, 1 or missing')
keep <- !is.na(flag)&flag==1
write_tab(data.frame(selection='Stored outcome_last5yr == 1; same rule as binary sensitivity',
 cohort_pregnancies=nrow(preg),eligible_pregnancies=sum(keep),excluded_flag0=sum(flag==0,na.rm=TRUE),
 missing_flag=sum(is.na(flag))),'Five_year_restriction','qc')
# Audit timing if dates are present; never replace the stored-flag selection.
if(all(c('v008','p3') %in% names(preg))) {
 lag <- preg$v008-preg$p3
 expected <- !is.na(lag)&lag>=0&lag<=60
 write_tab(data.frame(selected_by_stored_flag=sum(keep),selected_by_outcome_CMC_check=sum(expected),
  selection_disagreements=sum(keep!=expected),
  note='Audit only. Analysis selection uses the stored flag without recomputing its definition.'),
  'Five_year_flag_timing_audit','qc')
}
d <- d[as.character(preg_uid) %in% as.character(preg$preg_uid[keep])]
yy <- as.character(d$preg_start_year_f)
yy[yy %in% c('2015','2016')] <- '2015-2016'
d[,preg_start_year_f:=factor(yy)]
write_tab(d[,.(person_months=.N,events=sum(miscarriage_event)),by=preg_start_year_f],paste0('Five_year_',h$key,'_pooled_year_events'),'qc')

write_tab(counts(d),'Five_year_reported_observed_counts','qc')
rm(preg);gc()
# Preserve original country-IQR denominators. Mean-centering changes intercept,
# not the HR per country IQR. For DTS, centre over eligible person-month values.
original_var <- h$iqr
centre <- mean(d[[original_var]][is.finite(d[[original_var]])])
if(!is.finite(centre)) stop('No finite IQR exposures in five-year subset')
model_var <- paste0(original_var,'_centered')
d[, (model_var) := get(original_var)-centre]
write_tab(data.frame(Heat_metric=h$metric,original_variable=original_var,model_variable=model_var,
 centering_mean=centre,centering_unit='Restricted DTS person-months; finite original scaled values',
 IQR_denominators='Original saved country IQRs preserved'),paste0('Five_year_',h$key,'_IQR_centering'),'qc')
h_centered <- h;h_centered$iqr <- model_var
run_linear(d,h_centered,'IQR','Five_year_reported')
run_linear(d,h,'RAW','Five_year_reported')
rm(d);gc()
