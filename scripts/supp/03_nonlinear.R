# Fits GM1-GM3 on BOTH exposure scales; exports curves and pregnancy histograms.
h <- metric_row();d <- load_obs()
for(scale in c('IQR','RAW')) {
 var <- h[[tolower(scale)]]
 # Common plotting distribution: last OBSERVED interval in the GM3-complete sample.
 dd3 <- analysis_data(d,form(var,3),extra=c(var,'gest_month_f'))
 terminal <- dd3[order(as.integer(as.character(gest_month_f))),.SD[.N],by=preg_uid]
 distribution <- data.frame(pregnancy_id=as.character(terminal$preg_uid),Heat_metric=h$metric,exposure=terminal[[var]])
 write_tab(distribution,paste0('Nonlinear_',h$key,'_',scale,'_distribution'))
 lim <- range(dd3[[var]])
 basis0 <- splines::ns(dd3[[var]],df=3)
 spec <- list(knots=attr(basis0,'knots'),Boundary.knots=attr(basis0,'Boundary.knots'),intercept=FALSE,
  definition='3-df natural cubic spline; knots determined on GM3-complete person-month sample',
  reference=median(terminal[[var]]),grid=seq(quantile(terminal[[var]],.01),quantile(terminal[[var]],.99),length.out=200),
  distribution='Pregnancies: final observed eligible risk interval in GM3-complete sample')
 saveRDS(spec,file.path(OUT,'qc',paste0('Nonlinear_',h$key,'_',scale,'_spline_spec.rds')))
 basis_names <- paste0('heat_s',1:3)
 curves <- comparisons <- list()
 for(m in 1:3) {
  dd <- analysis_data(d,form(var,m),extra=var)
  B <- splines::ns(dd[[var]],knots=spec$knots,Boundary.knots=spec$Boundary.knots,intercept=FALSE)
  for(j in 1:3) set(dd,j=basis_names[j],value=B[,j])
  ff <- form(paste(basis_names,collapse=' + '),m)
  id <- paste('Nonlinear',h$key,scale,paste0('GM',m),sep='_')
  fit <- fit_saved(ff,dd,id)
  gridB <- splines::ns(spec$grid,knots=spec$knots,Boundary.knots=spec$Boundary.knots,intercept=FALSE)
  refB <- splines::ns(spec$reference,knots=spec$knots,Boundary.knots=spec$Boundary.knots,intercept=FALSE)
  contrast <- sweep(gridB,2,as.numeric(refB),'-')
  beta <- fixef(fit)$cond[basis_names];V <- vcov(fit)$cond[basis_names,basis_names]
  eta <- drop(contrast%*%beta);se <- sqrt(rowSums((contrast%*%V)*contrast))
  curves[[m]] <- data.frame(Model=paste('Model',m),Heat_metric=h$metric,exposure=spec$grid,
   estimate=exp(eta),lower_95=exp(eta-1.96*se),upper_95=exp(eta+1.96*se),reference_exposure=spec$reference,
   prediction_type='Conditional HR relative to common GM3-sample median; pointwise Wald CI')
  # Linear comparator refitted on IDENTICAL rows; no cross-sample AIC/LRT.
  if(scale=='IQR') {
   lin <- readRDS(primary_path(h$metric,m))
   if(!valid(diag_fit(lin)) || nobs(lin)!=nrow(dd)) stop('Primary linear comparator differs in validity/N')
   # Verify stored rows directly; row-count equality alone is insufficient.
   mf <- lin$frame
   comparable <- intersect(names(dd),names(mf))
   comparable <- setdiff(comparable,c(basis_names,'age_c'))
   for(v in comparable) if(!isTRUE(all.equal(as.character(dd[[v]]),as.character(mf[[v]])))) stop('Primary comparator row/value mismatch: ',v)
   if(m==3 && max(abs(dd$age_c+attr(dd,'age_centre')-mf$age_at_preg_start))>1e-8) stop('Primary comparator age mismatch')
   } else {
   rawpath <- file.path(OUT,'models',paste0('Raw_linear_',h$key,'_RAW_GM',m,'.rds'))
   if(file.exists(rawpath)) {
    lin <- readRDS(rawpath)
    if(!valid(diag_fit(lin)) || nobs(lin)!=nrow(dd)) stop('Raw linear comparator differs in validity/N')
    mf <- lin$frame;comparable <- setdiff(intersect(names(dd),names(mf)),basis_names)
    for(v in comparable) if(!isTRUE(all.equal(as.character(dd[[v]]),as.character(mf[[v]])))) stop('Raw comparator row/value mismatch: ',v)
   } else lin <- fit_saved(form(var,m),dd,paste0(id,'_matched_linear'))
  }
  ll1 <- logLik(lin);ll2 <- logLik(fit);stat <- 2*(as.numeric(ll2)-as.numeric(ll1));df0 <- attr(ll2,'df')-attr(ll1,'df')
  comparisons[[m]] <- data.frame(Model=paste('Model',m),Heat_metric=h$metric,scale=scale,
   N=nobs(fit),linear_AIC=AIC(lin),spline_AIC=AIC(fit),delta_AIC=AIC(fit)-AIC(lin),
   LRT=stat,df=df0,p_nonlinearity=if(stat>=0&&df0>0) pchisq(stat,df0,lower.tail=FALSE) else NA_real_)
  write_tab(rbindlist(curves),paste0('Nonlinear_',h$key,'_',scale,'_predictions'))
  write_tab(rbindlist(comparisons),paste0('Nonlinear_',h$key,'_',scale,'_comparisons'))
  rm(fit,lin,dd);gc()
 }
 write_tab(data.frame(Heat_metric=h$metric,reference_exposure=spec$reference),paste0('Nonlinear_',h$key,'_',scale,'_reference'))
 rm(dd3,terminal,basis0);gc()
}
rm(d);gc()
