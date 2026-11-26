# Saved-object audit of all 12 primary IQR fits, sequentially. No refitting.
d <- load_obs();rows <- list()
for(j in 1:nrow(HEAT)) for(m in 1:3) {
 h <- HEAT[j,];path <- primary_path(h$metric,m);id <- paste(h$key,paste0('GM',m),sep='_')
 ans <- tryCatch({
  fit <- readRDS(path)
  if(family(fit)$link!='cloglog') stop('Wrong link')
  actual <- paste(deparse(formula(fit)),collapse=' ')
  if(!grepl('country',actual)||!grepl('mother_uid',actual)||grepl('cook_fuel|marital',actual)) stop('Unexpected primary specification')
  ff <- formula(fit);dd <- analysis_data(d,ff)
  if(nrow(dd)!=nobs(fit)) stop('Reconstructed complete-case N differs from checkpoint')
  dg <- cbind(Model=paste('Model',m),Heat_metric=h$metric,source=path,diag_fit(fit),counts(dd)[,-1,drop=FALSE])
  save_gradient_details(fit,paste0('Primary_',id))
  dg$status <- if(valid(dg)) 'passed_stored_fit_checks' else 'requires_review'
  writeLines(capture.output(formula(fit)),file.path(OUT,'qc',paste0('Primary_',id,'_formula.txt')))
  writeLines(capture.output(VarCorr(fit)),file.path(OUT,'qc',paste0('Primary_',id,'_VarCorr.txt')))
  write_tab(data.frame(term=names(fixef(fit)$cond),estimate=fixef(fit)$cond,
   SE=sqrt(diag(vcov(fit)$cond))),paste0('Primary_',id,'_fixed_effects'))
  rm(fit,dd);gc();dg
 },error=function(e) data.frame(Model=paste('Model',m),Heat_metric=h$metric,source=path,status='requires_review',error=conditionMessage(e)))
 rows[[id]] <- ans;write_tab(rbindlist(rows,fill=TRUE),'Supplementary_Table_S3_primary_IQR_diagnostics')
}
rm(d);gc()

if(any(vapply(rows,function(x) x$status != "passed_stored_fit_checks",logical(1))))
 stop("S3 found a primary model requiring review; inspect Supplementary_Table_S3_primary_IQR_diagnostics.csv. Dependent jobs must wait.")
