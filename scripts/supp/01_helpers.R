suppressPackageStartupMessages({library(data.table);library(glmmTMB);library(lme4);library(digest)})
setDTthreads(1L)
Sys.setenv(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
write_tab <- function(x,name,folder='tables') fwrite(x,file.path(OUT,folder,paste0(name,'.csv')))
need <- function(d,v) {if(length(setdiff(v,names(d)))) stop('Missing variables: ',paste(setdiff(v,names(d)),collapse=', '))}
finite_rows <- function(d,v) {
 need(d,v);ok <- complete.cases(d[,..v])
 for (z in v) if(is.numeric(d[[z]])) ok <- ok & is.finite(d[[z]])
 ok
}
metric_row <- function() {j <- match(HEAT_RUN,HEAT$metric);if(is.na(j)) stop('Invalid HEAT_RUN');HEAT[j,]}
form <- function(heat,m,country=TRUE,re='(1 | mother_uid)') {
 terms <- c('gest_month_f',heat,if(country) 'country',if(m>=2) TEMP,if(m>=3) SOC,re)
 as.formula(paste('miscarriage_event ~',paste(terms,collapse=' + ')))
}
load_obs <- function() {
 d <- as.data.table(readRDS(OBS));need(d,c('preg_uid','mother_uid','cluster_uid','miscarriage_event','gest_month_f','country','age_at_preg_start'))
 d[,age_c:=age_at_preg_start] # actual centring is performed on each analysis sample
 for(v in intersect(c('country','mother_uid','cluster_uid','gest_month_f',TEMP,'rc_empl','v106','v190','rc_media_none','rc_hins_any','v025'),names(d))) set(d,j=v,value=factor(d[[v]]))
 d
}
analysis_data <- function(d,f,extra=character()) {
 v <- unique(c(all.vars(f),extra,'preg_uid','cluster_uid','mother_uid'))
 out <- copy(d[finite_rows(d,v),..v])
 centre <- if('age_c' %in% v) mean(out$age_c) else NA_real_
 if('age_c' %in% v) out[,age_c:=age_c-centre]
 for(vv in names(out)) if(is.factor(out[[vv]])) set(out,j=vv,value=droplevels(out[[vv]]))
 setattr(out,'age_centre',centre)
 if(!nrow(out) || sum(out$miscarriage_event)==0) stop('No usable data/events.')
 out
}
counts <- function(d) data.frame(N=nrow(d),pregnancies=uniqueN(d$preg_uid),mothers=uniqueN(d$mother_uid),clusters=uniqueN(d$cluster_uid),events=if("miscarriage_event" %in% names(d)) sum(d$miscarriage_event) else NA_real_)
diag_fit <- function(fit) {
 b <- fixef(fit)$cond;V <- vcov(fit)$cond;se <- suppressWarnings(sqrt(diag(V)))
 g <- fit$sdr$gradient.fixed;C <- fit$sdr$cov.fixed
 scaled <- NA_real_
 if(length(g) && is.matrix(C) && nrow(C)==length(g) && all(is.finite(c(g,C)))) {
  sd0 <- suppressWarnings(sqrt(diag(C)))
  if(all(is.finite(sd0)&sd0>0)) scaled <- max(abs(as.vector(C%*%g)/sd0))
 }
 ll <- logLik(fit)
 data.frame(convergence_code=fit$fit$convergence,pdHess=isTRUE(fit$sdr$pdHess),
 all_fixed_SE_finite=all(is.finite(c(b,V,se)))&&all(se>0),
 max_abs_raw_gradient=if(length(g)) max(abs(g)) else NA_real_,
 max_abs_approximate_change_in_SE_units=scaled,logLik=as.numeric(ll),AIC=AIC(fit),BIC=BIC(fit),
 parameters=attr(ll,'df'),N=nobs(fit),optimizer_message=paste(fit$fit$message,collapse=';'))
}
valid <- function(dg) isTRUE(dg$convergence_code==0)&&isTRUE(dg$pdHess)&&isTRUE(dg$all_fixed_SE_finite)&&all(is.finite(c(dg$logLik,dg$AIC)))
save_gradient_details <- function(fit,id) {
 p <- fit$fit$par;g <- fit$sdr$gradient.fixed;V <- fit$sdr$cov.fixed
 if(length(g)!=length(p) || !is.matrix(V) || nrow(V)!=length(p) || any(!is.finite(c(g,V)))) return(invisible(NULL))
 se <- suppressWarnings(sqrt(diag(V)));step <- -drop(V%*%g)
 labels <- names(p);ib <- which(labels=='beta')
 if(length(ib)==length(fixef(fit)$cond)) labels[ib] <- names(fixef(fit)$cond)
 write_tab(data.frame(parameter=labels,estimate=unname(p),raw_gradient=g,parameter_SE=se,
  approximate_Newton_step=step,approximate_change_in_SE_units=step/se),paste0(id,'_gradient_details'),'qc')
}
# Cache signature checks formula, full analytic model data, and optimizer.
fit_saved <- function(f,d,id) {
 path <- file.path(OUT,'models',paste0(id,'.rds'));meta <- paste0(path,'.meta.rds')
 sig <- digest(list(formula=paste(deparse(f),collapse=' '),data=as.data.frame(d),optimizer=OPTIMIZER),algo='xxhash64')
 if(file.exists(path)) {
  if(!file.exists(meta) || !identical(readRDS(meta)$signature,sig)) stop('Checkpoint signature differs: ',path)
  fit <- readRDS(path)
 } else {
  ctl <- switch(OPTIMIZER,nlminb=glmmTMBControl(optCtrl=list(iter.max=10000L,eval.max=20000L)),
   BFGS=glmmTMBControl(optimizer=optim,optArgs=list(method='BFGS'),optCtrl=list(maxit=20000L)),stop('Invalid OPTIMIZER'))
  warnings0 <- character();started <- Sys.time()
  message('Fitting ',id,'; N=',nrow(d),'; ',Sys.time())
  fit <- withCallingHandlers(glmmTMB(f,data=as.data.frame(d),family=binomial('cloglog'),control=ctl),
   warning=function(w) {warnings0 <<- c(warnings0,conditionMessage(w));invokeRestart('muffleWarning')})
  saveRDS(fit,path,compress=FALSE)
  saveRDS(list(signature=sig,formula=f,age_centre=attr(d,'age_centre'),counts=counts(d),optimizer=OPTIMIZER,
   elapsed_minutes=as.numeric(difftime(Sys.time(),started,units='mins')),created=Sys.time()),meta)
  writeLines(unique(warnings0),file.path(OUT,'qc',paste0(id,'_warnings.txt')))
 }
 dg <- cbind(model=id,diag_fit(fit),counts(d)[,-1,drop=FALSE])
 write_tab(dg,paste0(id,'_diagnostics'),'qc')
 save_gradient_details(fit,id)
 writeLines(capture.output(VarCorr(fit)),file.path(OUT,'qc',paste0(id,'_VarCorr.txt')))
 writeLines(capture.output(f),file.path(OUT,'qc',paste0(id,'_formula.txt')))
 b <- fixef(fit)$cond;V <- vcov(fit)$cond
 write_tab(data.frame(term=names(b),estimate=b,SE=suppressWarnings(sqrt(diag(V)))),paste0(id,'_fixed_effects'))
 saveRDS(V,file.path(OUT,'qc',paste0(id,'_fixed_covariance.rds')))
 writeLines(capture.output(sessionInfo()),file.path(OUT,'qc',paste0(id,'_sessionInfo.txt')))
 # No diagnose() Richardson numerical Hessian: retained stored gradient/covariance only.
 if(!valid(dg)) stop('Fit failed inferential checks; checkpoint/diagnostics retained: ',id)
 fit
}
heat_result <- function(fit,var,metric,m,scale,id) {
 b <- fixef(fit)$cond[var];se <- sqrt(vcov(fit)$cond[var,var])
 if(length(b)!=1L || !is.finite(b+se) || se<=0) stop('Heat coefficient not estimable')
 data.frame(analysis=id,Model=paste('Model',m),Heat_metric=metric,Exposure_scale=scale,
 estimate_log_HR=unname(b),SE=unname(se),HR=exp(b),lower_95=exp(b-1.96*se),upper_95=exp(b+1.96*se),
 p_value=2*pnorm(abs(b/se),lower.tail=FALSE),N_person_months=nobs(fit),row.names=NULL)
}
run_linear <- function(d,h,scale,prefix) {
 var <- h[[tolower(scale)]]; rows <- list()
 for(m in 1:3) {
  f <- form(var,m);dd <- analysis_data(d,f);id <- paste(prefix,h$key,scale,paste0('GM',m),sep='_')
  fit <- fit_saved(f,dd,id);ans <- heat_result(fit,var,h$metric,m,if(scale=='IQR') 'Country-specific IQR' else 'Raw (per degree Celsius)',id)
  if(prefix=='Five_year_reported') ans$calendar_year_adjustment <- if(m>=2)'2015-2016 pooled; remaining initiation years categorical' else 'Not included in Model 1'
  rows[[m]] <- ans;write_tab(ans,paste0(id,'_HR'));rm(fit,dd);gc()
 }
 write_tab(rbindlist(rows),paste(prefix,h$key,scale,'GM1_GM2_GM3',sep='_'))
}
