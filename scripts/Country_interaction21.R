# Overall country-by-heat tests for the expanded 21-country DTS analysis.
# Read-only access to canonical primary models. New results use a separate folder.
suppressPackageStartupMessages({library(glmmTMB); library(data.table)})
Sys.setenv(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
root <- Sys.getenv("DTS_ROOT", file.path(getwd(),"output"))
action <- Sys.getenv("DTS_COUNTRY_TEST_ACTION", "fit")
metrics <- c("WBGT", "UTCI", "MRT", "Heat Index")
out <- file.path(root, "Supplementary_DTS_v1", "Country_interaction21")
for (s in c("models", "tables", "qc")) dir.create(file.path(out,s),recursive=TRUE,showWarnings=FALSE)
files <- setNames(file.path(out,"tables",paste0(c("WBGT","UTCI","MRT","HI"),"_country_interaction_LRT.csv")),metrics)
if (action=="collect") {
  inventory <- data.table(Heat_metric=metrics,file=unname(files),available=file.exists(files))
  fwrite(inventory,file.path(out,"qc","Collection_inventory.csv"))
  if (!all(inventory$available)) stop("Not all four validated tests are available; inspect Collection_inventory.csv and array logs.")
  z <- rbindlist(lapply(files,fread),use.names=TRUE)
  stopifnot(nrow(z)==4L,setequal(z$Heat_metric,metrics),all(z$status=="passed"),all(z$df==20L),all(is.finite(z$p_value)))
  z <- z[match(c("UTCI","MRT","WBGT","Heat Index"),Heat_metric)]
  z[,p_Holm:=p.adjust(p_value,method="holm",n=4L)]
  fwrite(z,file.path(out,"tables","Country_interaction_GM3_IQR_combined.csv"))
  display <- z[,.(Heat_metric,Person_months=N,Countries=countries,LR_statistic=round(LRT,3),df,
                  P_value=ifelse(p_value<0.001,"<0.001",sprintf("%.3f",p_value)),
                  P_Holm=ifelse(p_Holm<0.001,"<0.001",sprintf("%.3f",p_Holm)))]
  fwrite(display,file.path(out,"tables","Country_interaction_GM3_IQR_display.csv"))
  writeLines(paste("All four tests collected:",Sys.time()),file.path(out,"qc","COUNTRY_INTERACTION_OK.txt"))
  print(display); quit(save="no",status=0L)
}
if (action!="fit") stop("DTS_COUNTRY_TEST_ACTION must be fit or collect.")
metric <- Sys.getenv("DTS_HEAT", "")
if (!metric %in% metrics) stop("DTS_HEAT must identify one of the four heat metrics.")
key <- c(WBGT="WBGT",UTCI="UTCI",MRT="MRT","Heat Index"="HI")[[metric]]
heat <- c(WBGT="WBGT_T1_country_iqr",UTCI="utci_T1_country_iqr",MRT="mrt_T1_country_iqr","Heat Index"="HI_T1_country_iqr")[[metric]]
qc <- file.path(out,"qc",key); dir.create(qc,recursive=TRUE,showWarnings=FALSE)
writeLines(capture.output(sessionInfo()),file.path(qc,"sessionInfo.txt"))
# Prevent a failed rerun leaving an earlier validated-looking result row.
if (file.exists(files[[metric]])) stop("A result already exists for ",metric,". No refit needed; collect it or archive that result explicitly before rerunning.")
source_table <- file.path(root,"Supplementary_DTS_v1","tables","Supplementary_Table_S3_primary_IQR_diagnostics.csv")
canonical <- file.path(root,"models","Linear_Nonlinear",paste0("GM3_",gsub(" ","_",metric),"_DTS_IQR_CountryFE_MotherRE_",if(metric=="Heat Index")"BFGS" else "NoCountryMonthInt",".rds"))
if (file.exists(source_table)) {
  s <- fread(source_table); pick <- s[Heat_metric==metric & Model=="Model 3"]
  if (nrow(pick)!=1L || !"source" %in% names(pick)) stop("S3 does not uniquely identify the primary GM3 checkpoint.")
  recorded <- pick$source[[1L]]
  if (file.exists(recorded)) {
    canonical <- recorded
  } else if (file.exists(file.path(root,"models","Linear_Nonlinear",basename(recorded)))) {
    canonical <- file.path(root,"models","Linear_Nonlinear",basename(recorded))
  } else stop("S3 checkpoint is unavailable; update its path or supply the canonical model.")
}
if (!file.exists(canonical)) stop("Primary GM3 checkpoint not found: ",canonical)
message("Loading primary model: ",canonical,"; ",Sys.time())
base <- readRDS(canonical)
if (!inherits(base,"glmmTMB")) stop("Primary checkpoint is not a glmmTMB fit.")
fam <- family(base)
if (fam$family!="binomial" || fam$link!="cloglog") stop("Expected binomial cloglog primary model.")
if (isTRUE(base$modelInfo$REML)) stop("This test requires ML, not REML.")
f0 <- formula(base,component="cond")
fixed0 <- formula(base,fixed.only=TRUE,component="cond")
expected <- c("gest_month_f",heat,"country","preg_start_year_f","preg_start_month_f",
              "age_at_preg_start","I(age_at_preg_start^2)","rc_empl","v106","v190","rc_media_none","rc_hins_any","v025")
normalise <- function(x) gsub("[[:space:]]","",x)
if (!setequal(normalise(attr(terms(fixed0),"term.labels")),normalise(expected))) stop("Primary fixed effects differ from the expected manuscript GM3 specification. Review before fitting.")
vc0 <- VarCorr(base)$cond
if (length(vc0)!=1L || !identical(names(vc0),"mother_uid") || nrow(vc0[[1]])!=1L) stop("Expected only a mother random intercept.")
if (length(fixef(base)$zi)>0L) stop("Zero inflation was not expected.")
dat <- base$frame
if (is.null(dat) || !is.data.frame(dat)) stop("The stored primary model has no usable model frame.")
if ("(weights)" %in% names(dat) || "(offset)" %in% names(dat)) stop("Unexpected weights or offset: review the baseline specification before fitting.")
if (!is.factor(dat$country) || !is.factor(dat$mother_uid)) stop("Country and mother must retain their fitted factor encodings.")
counts <- as.data.table(table(dat$country));setnames(counts,c("country","person_months"))
fwrite(counts,file.path(qc,"Country_sample_counts.csv"))
if (sum(counts$person_months>0L)!=21L || any(counts$person_months==0L)) stop("Expected exactly 21 nonempty country levels.")
expected_N <- if (metric %in% c("WBGT","Heat Index"))4823692L else 4823874L
expected_events <- if (metric %in% c("WBGT","Heat Index"))41000L else 41004L
if (nrow(dat)!=expected_N || sum(dat$miscarriage_event)!=expected_events) stop("Primary sample counts differ from the reviewed 21-country analysis.")
diagnose_fit <- function(fit,label) {
  beta <- fixef(fit)$cond; v <- as.matrix(vcov(fit)$cond)
  se <- suppressWarnings(sqrt(diag(v)))
  g <- fit$sdr$gradient.fixed; V <- fit$sdr$cov.fixed
  scaled <- NA_real_
  if (length(g)>0L && is.matrix(V) && nrow(V)==length(g) && all(is.finite(g)) &&
      all(is.finite(V)) && all(diag(V)>0)) {
    scaled <- max(abs(as.numeric(V %*% g)/sqrt(diag(V))))
  }
  ll <- tryCatch(logLik(fit),error=function(e)NA_real_)
  k <- attr(ll,"df"); if(is.null(k)||length(k)!=1L)k<-NA_integer_
  msg <- fit$fit$message; if(is.null(msg)||length(msg)==0L)msg<-NA_character_
  data.table(model=label,convergence_code=as.integer(fit$fit$convergence),
    optimizer_message=as.character(msg)[1L],pdHess=isTRUE(fit$sdr$pdHess),
    all_fixed_estimates_finite=length(beta)>0L && all(is.finite(beta)),
    all_fixed_SE_finite=length(se)==length(beta) && all(is.finite(se)) && all(se>0),
    max_abs_raw_gradient=if(length(g)&&all(is.finite(g)))max(abs(g)) else NA_real_,
    max_abs_approximate_change_in_SE_units=scaled,logLik=as.numeric(ll),
    K=as.integer(k),AIC=tryCatch(AIC(fit),error=function(e)NA_real_),
    BIC=tryCatch(BIC(fit),error=function(e)NA_real_),N=as.integer(nobs(fit)))
}
valid <- function(x) isTRUE(x$convergence_code==0L) && isTRUE(x$pdHess) &&
  isTRUE(x$all_fixed_estimates_finite) && isTRUE(x$all_fixed_SE_finite) &&
  is.finite(x$logLik) && is.finite(x$max_abs_approximate_change_in_SE_units) &&
  x$max_abs_approximate_change_in_SE_units<=0.05
d0 <- diagnose_fit(base,"Primary_GM3_common_slope")
fwrite(d0,file.path(qc,"Null_diagnostics.csv"));print(d0)
if (!valid(d0)) stop("Primary checkpoint requires diagnostic review; it has not been modified.")
f1 <- update.formula(f0,paste(". ~ . + country:",heat))
fixed1 <- update.formula(fixed0,paste(". ~ . + country:",heat))
environment(f1)<-baseenv(); environment(fixed1)<-baseenv()
writeLines(c("PRIMARY",deparse(f0),"COUNTRY INTERACTION",deparse(f1)),file.path(qc,"Model_formulas.txt"))
# Use one row with all original factor levels retained to obtain coefficient
# names without constructing an extra 4.8-million-row design matrix.
contr <- base$modelInfo$contrasts
x0 <- model.matrix(fixed0,data=dat[1L,,drop=FALSE],contrasts.arg=contr)
x1 <- model.matrix(fixed1,data=dat[1L,,drop=FALSE],contrasts.arg=contr)
b0 <- fixef(base)$cond
if (!setequal(names(b0),colnames(x0)) || ncol(x1)-ncol(x0)!=20L) stop("Unexpected contrasts or interaction dimension.")
start_beta <- setNames(rep(0,ncol(x1)),colnames(x1));start_beta[names(b0)]<-b0
initial <- list(beta=unname(start_beta),theta=unname(getME(base,"theta")),b=unname(getME(base,"b")))
same_frame <- function(fit) {
  ff<-fit$frame
  is.data.frame(ff) && nrow(ff)==nrow(dat) && all(names(dat)%in%names(ff)) &&
    all(vapply(names(dat),function(v)identical(dat[[v]],ff[[v]]),logical(1L)))
}
validate_comparison <- function(fit,dg) {
  valid(dg) && same_frame(fit) && setequal(names(fixef(fit)$cond),colnames(x1)) &&
    isTRUE(dg$K-d0$K==20L) && isTRUE(dg$logLik>=d0$logLik-1e-4)
}
chosen <- NULL; d1 <- NULL; chosen_optimizer <- NA_character_
# A bounded recovery attempt keeps the model, data and diagnostic thresholds fixed.
for (optimizer in c("nlminb","BFGS")) {
  checkpoint<-file.path(out,"models",paste0(key,"_GM3_country_interaction_",optimizer,".rds"))
  warnings_seen<-character()
  fit<-tryCatch({
    if(file.exists(checkpoint)) {
      message("Reading interaction checkpoint: ",checkpoint)
      readRDS(checkpoint)
    } else {
      message("Fitting ",metric," country interactions with ",optimizer,"; ",Sys.time())
      ctl<-if(optimizer=="BFGS")glmmTMBControl(optimizer=optim,optArgs=list(method="BFGS"),
        optCtrl=list(maxit=20000L,reltol=1e-10),rank_check="stop",parallel=1L) else
        glmmTMBControl(optCtrl=list(iter.max=10000L,eval.max=20000L),rank_check="stop",parallel=1L)
      f<-withCallingHandlers(glmmTMB(formula=f1,data=dat,family=binomial(link="cloglog"),
        REML=FALSE,na.action=na.fail,contrasts=contr,start=initial,control=ctl),
        warning=function(w){warnings_seen<<-c(warnings_seen,conditionMessage(w));invokeRestart("muffleWarning")})
      saveRDS(f,checkpoint,compress=FALSE);f
    }
  },error=function(e){writeLines(conditionMessage(e),file.path(qc,paste0(optimizer,"_error.txt")));NULL})
  writeLines(warnings_seen,file.path(qc,paste0(optimizer,"_warnings.txt")))
  if(!is.null(fit)) {
    dg<-tryCatch(diagnose_fit(fit,paste0("Country_interaction_",optimizer)),error=function(e){
      writeLines(conditionMessage(e),file.path(qc,paste0(optimizer,"_diagnostic_error.txt")))
      NULL
    })
    if(is.null(dg)){rm(fit);invisible(gc());next}
    dg[,`:=`(identical_primary_frame=same_frame(fit),expected_fixed_terms=setequal(names(fixef(fit)$cond),colnames(x1)),
              df_added=K-d0$K,logLik_change=logLik-d0$logLik)]
    fwrite(dg,file.path(qc,paste0(optimizer,"_diagnostics.csv")));print(dg)
    if(validate_comparison(fit,dg)){chosen<-fit;d1<-dg;chosen_optimizer<-optimizer;break}
    if(optimizer=="nlminb" && same_frame(fit) && setequal(names(fixef(fit)$cond),colnames(x1))) {
      attempt<-tryCatch(list(beta=unname(fixef(fit)$cond[colnames(x1)]),
        theta=unname(getME(fit,"theta")),b=unname(getME(fit,"b"))),error=function(e)NULL)
      if(!is.null(attempt)&&all(vapply(attempt,function(v)all(is.finite(v)),logical(1L))))initial<-attempt
    }
    rm(fit);invisible(gc())
  }
}
if(is.null(chosen))stop("Interaction fit failed validation after the bounded recovery attempt. Review retained checkpoints and QC; no P value exported.")
lr<-2*(d1$logLik-d0$logLik)
if(lr<0){writeLines(sprintf("Tiny negative LR before truncation: %.12g",lr),file.path(qc,"Likelihood_roundoff.txt"));lr<-0}
result<-data.table(Heat_metric=metric,N=nrow(dat),events=sum(dat$miscarriage_event),countries=21L,
  null_K=d0$K,interaction_K=d1$K,null_logLik=d0$logLik,interaction_logLik=d1$logLik,
  null_AIC=d0$AIC,interaction_AIC=d1$AIC,delta_AIC=d1$AIC-d0$AIC,
  LRT=lr,df=d1$K-d0$K,p_value=pchisq(lr,df=d1$K-d0$K,lower.tail=FALSE),
  optimizer=chosen_optimizer,primary_checkpoint=canonical,identical_primary_frame=TRUE,status="passed")
fwrite(result,files[[metric]])
co<-summary(chosen)$coefficients$cond
fwrite(data.table(term=rownames(co),as.data.frame(co)),file.path(out,"tables",paste0(key,"_interaction_fixed_effects.csv")))
writeLines(paste("Validated interaction test completed:",Sys.time()),file.path(qc,"TEST_OK.txt"))
print(result)
