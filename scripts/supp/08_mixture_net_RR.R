# Four single-heat + five-air mixtures only. No all-heat mixture is fitted.
# Model-standardized cumulative NET risk over full trajectories, not interval RR.
# Source settings/helpers. Choose package (existing GLM) or mixed (new mother GLMM).
suppressPackageStartupMessages({library(sandwich);library(statmod)})
if(is.na(GH_N)||GH_N<5L) stop('DTS_GH_N must be at least 5')
h <- metric_row();vars <- c(h$raw,AIR);key <- paste0('Air_',h$key)
prepared <- file.path(ROOT,'PIF_QGComp_v3','prepared')
qdat <- as.data.table(readRDS(file.path(prepared,'Mixture_common_quartiles.rds')))
need(qdat,c(vars,'preg_uid','mother_uid','cluster_uid','miscarriage_event'))
if(!MIXTURE_ENGINE %in% c('package','mixed')) stop('Invalid MIXTURE_ENGINE')
if(MIXTURE_ENGINE=='package') {
 source_path <- file.path(ROOT,'PIF_QGComp_v3','models',paste0('package_',key,'_cloglog.rds'))
 qfit <- readRDS(source_path);fit <- qfit$fit
 if(!inherits(fit,'glm') || family(fit)$link!='cloglog' || !isTRUE(fit$converged)) stop('Invalid existing cloglog GLM')
 if(nobs(fit)!=nrow(qdat)) stop('Existing GLM N differs from prepared sample')
 b <- coef(fit);V <- vcovCL(fit,cluster=qdat$cluster_uid,type='HC0',cadjust=TRUE)
 if(any(!is.finite(c(b,V)))) stop('Invalid GLM coefficients/covariance')
 ff <- formula(fit);sigma <- 0;centre <- NA_real_
 write_tab(data.frame(model=key,engine='package',N=nobs(fit),converged=fit$converged,rank=fit$rank,
  columns=length(b),logLik=as.numeric(logLik(fit)),AIC=AIC(fit),BIC=BIC(fit),
  uncertainty='DHS-cluster HC0 robust delta method'),paste0(key,'_package_net_RR_fit_diagnostics'),'qc')
} else {
 if(!'age_c' %in% names(qdat)) qdat[,age_c:=age_at_preg_start]
 ff <- form(paste(vars,collapse=' + '),3)
 dd <- analysis_data(qdat,ff,extra=vars)
 if(nrow(dd)!=nrow(qdat)) stop('Mixed comparison has a different analytic sample')
 centre <- attr(dd,'age_centre')
 fit <- fit_saved(ff,dd,paste0('Mixture_mixed_',key))
 b <- fixef(fit)$cond
 p <- fit$fit$par;ib <- which(names(p)=='beta');it <- which(names(p)=='theta')
 if(length(ib)!=length(b)||length(it)!=1||length(p)!=length(b)+1) stop('Expected beta plus one mother log-SD')
 if(max(abs(p[ib]-b))>1e-9) stop('Stored fixed parameter mismatch')
 V <- fit$sdr$cov.fixed[c(ib,it),c(ib,it),drop=FALSE]
 if(!all(is.finite(V))) stop('Invalid joint beta/log-SD covariance')
 sigma <- exp(p[it]);rm(dd)
}
id <- paste0(key,'_',MIXTURE_ENGINE,'_net_RR')
ff <- lme4::nobars(ff)
vc_hr <- if(MIXTURE_ENGINE=='mixed') vcov(fit)$cond else V
psi <- sum(b[vars]);sepsi <- sqrt(sum(vc_hr[vars,vars,drop=FALSE]))
write_tab(data.frame(mixture=key,engine=MIXTURE_ENGINE,HR=exp(psi),lower_95=exp(psi-1.96*sepsi),upper_95=exp(psi+1.96*sepsi),
 estimand='Conditional hazard ratio per simultaneous one-quartile increase'),paste0(id,'_conditional_HR'))
w <- data.frame(mixture=key,engine=MIXTURE_ENGINE,component=vars,coefficient=b[vars],direction=ifelse(b[vars]>0,'Positive',ifelse(b[vars]<0,'Negative','Zero')),weight=0)
for(dr in c('Positive','Negative')) {ix <- w$direction==dr;if(any(ix)) w$weight[ix] <- abs(w$coefficient[ix])/sum(abs(w$coefficient[ix]))}
write_tab(w,paste0(id,'_directional_weights'))
# Save full covariance used for beta/log-SD delta propagation.
saveRDS(list(beta=b,sigma=sigma,covariance=V,age_centre=centre,formula=ff),file.path(OUT,'qc',paste0(id,'_parameters.rds')))
# Full trajectories are required. Observed follow-up ending at an event cannot
# simply be multiplied to estimate a common-horizon standardized net risk.
target <- as.data.table(readRDS(FULL))
if(MIXTURE_ENGINE=='mixed') target[,age_c:=age_at_preg_start-centre]
required <- unique(c(all.vars(delete.response(terms(ff))),vars,'preg_uid','mother_uid','cluster_uid','gest_month_f'))
need(target,required)
target <- target[as.character(preg_uid) %in% as.character(unique(qdat$preg_uid)),..required]
# Set factor levels/reference contrasts to the fitted analysis data.
for(v in intersect(names(target),names(qdat))) if(is.factor(qdat[[v]])) set(target,j=v,value=factor(target[[v]],levels=levels(qdat[[v]])))
bad <- unique(target$preg_uid[!finite_rows(target,required)])
target <- target[!preg_uid %in% bad]
setorder(target,preg_uid,gest_month_f)
horizon <- sort(unique(as.integer(as.character(qdat$gest_month_f))))
ck <- target[,.(n=.N,labels=paste(sort(as.integer(as.character(gest_month_f))),collapse=',')),by=preg_uid]
if(!nrow(target)||any(ck$n!=length(horizon))||any(ck$labels!=paste(horizon,collapse=','))) stop('Prediction trajectories do not match the full model horizon')
# Quantize observed target exposure with EXISTING fit cutpoints; tails extend to Inf.
cutpoints <- fread(file.path(ROOT,'PIF_QGComp_v3','qc','Mixture_quartile_cutpoints.csv'))
for(v in vars) {
 br <- cutpoints[variable==v][order(probability)]$value
 if(length(br)!=5L||any(diff(br)<=0)) stop('Invalid saved cutpoints for ',v)
 set(target,j=v,value=as.integer(cut(target[[v]],breaks=c(-Inf,br[2:4],Inf),labels=FALSE,right=TRUE))-1L)
}
write_tab(cbind(mixture=key,engine=MIXTURE_ENGINE,horizon=paste(horizon,collapse=','),counts(target),excluded_incomplete_pregnancies=length(bad)),paste0(id,'_target_counts'),'qc')
# Report empirical joint-category support: interventions all0..all3 may be sparse.
write_tab(data.frame(category=0:3,n_observed_rows_all_components_equal=vapply(0:3,function(k)
 sum(rowSums(as.matrix(qdat[,..vars])==k)==length(vars)),numeric(1))),paste0(id,'_joint_support'),'qc')

# Compute risk and analytic derivative wrt beta and log-SD in pregnancy chunks.
# Conditional cloglog survival over a path: S=exp(-sum(exp(eta_t))*exp(b_m)).
# Integrate ONE shared mother intercept over the whole pregnancy path.
# Derivative includes log-SD uncertainty for mixed fits.
standard_risk <- function(k,nodes_n) {
 gh <- if(sigma==0) list(nodes=0,weights=1) else statmod::gauss.quad.prob(nodes_n,dist='normal')
 ids <- unique(target$preg_uid);acc <- 0;grad <- numeric(nrow(V));n <- 0
 for(start in seq.int(1L,length(ids),by=CHUNK_PREG)) {
  end <- min(start+CHUNK_PREG-1L,length(ids))
  ix <- seq.int((start-1L)*length(horizon)+1L,end*length(horizon))
  z <- copy(target[ix])
  if(!is.na(k)) for(v in vars) set(z,j=v,value=rep(k,nrow(z)))
  X <- model.matrix(delete.response(terms(ff)),z)
  if(!setequal(colnames(X),names(b))) stop('Prediction matrix differs from fitted coefficients')
  X <- X[,names(b),drop=FALSE]
  expeta <- exp(drop(X%*%b));if(any(!is.finite(expeta))) stop('Overflow in cumulative hazard')
  group <- as.character(z$preg_uid)
  L <- drop(rowsum(matrix(expeta,ncol=1),group,reorder=FALSE))
  dL <- rowsum(X*expeta,group,reorder=FALSE)
  pr <- derivL <- derivTheta <- numeric(length(L))
  for(j in seq_along(gh$nodes)) {
   e <- exp(sigma*gh$nodes[j]);u <- L*e;s <- exp(-u)
   pr <- pr+gh$weights[j]*(-expm1(-u))
   derivL <- derivL+gh$weights[j]*s*e
   if(sigma>0) derivTheta <- derivTheta+gh$weights[j]*s*u*sigma*gh$nodes[j]
  }
  acc <- acc+sum(pr);grad[seq_along(b)] <- grad[seq_along(b)]+colSums(dL*derivL)
  if(MIXTURE_ENGINE=='mixed') grad[length(grad)] <- grad[length(grad)]+sum(derivTheta)
  n <- n+length(L)
 }
 list(risk=acc/n,gradient=grad/n)
}
risks <- lapply(c(NA_integer_,0:3),function(k) standard_risk(k,GH_N))
names(risks) <- c('observed',paste0('joint_',0:3))
# Quadrature check includes the gradient used for uncertainty, not only risk.
checks <- if(sigma==0) risks else lapply(c(NA_integer_,0:3),function(k) standard_risk(k,GH_N+10L))
qcheck <- data.frame(scenario=names(risks),risk_absolute_difference=vapply(seq_along(risks),function(j) abs(risks[[j]]$risk-checks[[j]]$risk),numeric(1)),
 max_gradient_absolute_difference=vapply(seq_along(risks),function(j) max(abs(risks[[j]]$gradient-checks[[j]]$gradient)),numeric(1)))
write_tab(qcheck,paste0(id,'_quadrature_check'),'qc')
if(any(qcheck$risk_absolute_difference>1e-5)||any(qcheck$max_gradient_absolute_difference>1e-5)) stop('Quadrature check failed: increase DTS_GH_N and rerun prediction stage')
ref <- risks$joint_0
results <- rbindlist(lapply(seq_along(risks),function(j) {
 r <- risks[[j]];g <- r$gradient/r$risk-ref$gradient/ref$risk
 variance <- drop(t(g)%*%V%*%g)
 if(!is.finite(variance)||variance < -1e-10) stop('Invalid risk-ratio variance')
 se <- sqrt(max(0,variance));rr <- r$risk/ref$risk
 data.frame(mixture=key,Heat_metric=h$metric,engine=MIXTURE_ENGINE,scenario=names(risks)[j],
 standardized_net_risk=r$risk,reference_net_risk=ref$risk,RR=rr,lower_95=exp(log(rr)-1.96*se),upper_95=exp(log(rr)+1.96*se),SE_log_RR=se,
 horizon=paste(horizon,collapse=','),target='Complete binary-comparator prediction trajectories; equally weighted pregnancies',
 uncertainty=if(MIXTURE_ENGINE=='package') 'Cluster HC0 delta; empirical target and cutpoints treated fixed' else 'Joint beta/log-SD model-based delta; empirical target and cutpoints treated fixed',
 estimand='Population-standardized cumulative NET risk ratio versus all mixture components at category zero; competing outcomes eliminated')
}))
write_tab(results,paste0(id,'_results'))
writeLines(capture.output(sessionInfo()),file.path(OUT,'qc',paste0(id,'_sessionInfo.txt')))
rm(qdat,target,fit);gc()
