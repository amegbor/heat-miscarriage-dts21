h <- metric_row();d <- load_obs()
countries <- sort(unique(as.character(d$country)))
if(COUNTRY_RUN!='ALL') {if(!COUNTRY_RUN %in% countries) stop('Unknown country');countries <- COUNTRY_RUN}
rows <- omitted <- list()
for(cc in countries) {
 id <- paste('Country',cc,h$key,'GM3',sep='_')
 rows[[cc]] <- tryCatch({
  z <- d[as.character(country)==cc]
  ff <- form(h$iqr,3,country=FALSE)
  dd <- analysis_data(z,ff)
  # Constant within-country terms are non-estimable; record every omission.
  covs <- c('gest_month_f',TEMP,'age_c','rc_empl','v106','v190','rc_media_none','rc_hins_any','v025')
  constants <- covs[vapply(covs,function(v) uniqueN(dd[[v]])<2,logical(1))]
  if(uniqueN(dd[[h$iqr]])<2) stop('No within-country heat variation')
  ff <- form(h$iqr,3,country=FALSE)
  for(v in constants) {
   ff <- update(ff,as.formula(paste('. ~ . -',v)))
   if(v=='age_c') ff <- update(ff,.~.-I(age_c^2))
  }
  omitted[[cc]] <- data.frame(country=cc,omitted_terms=paste(constants,collapse=';'))
  fit <- fit_saved(ff,dd,id)
  ans <- heat_result(fit,h$iqr,h$metric,3,'Country-specific IQR',id)
  ans$country <- cc;ans$status <- 'passed';rm(fit,dd,z);gc();ans
 },error=function(e) data.frame(country=cc,Heat_metric=h$metric,status='no_stable_estimate',error=conditionMessage(e)))
 write_tab(rbindlist(rows,fill=TRUE),paste0('Country_',h$key,'_GM3_results'))
 if(length(omitted)) write_tab(rbindlist(omitted),paste0('Country_',h$key,'_omitted_terms'),'qc')
}
rm(d);gc()
