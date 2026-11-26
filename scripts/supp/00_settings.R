# Run stages separately in OOD RStudio; no automatic submissions.
ROOT <- Sys.getenv('DTS_ROOT',file.path(getwd(),'output'))
OUT <- Sys.getenv('DTS_SUPP_OUT',file.path(ROOT,'Supplementary_DTS_v1'))
OBS <- file.path(ROOT,'data','DTS_ObservedRiskSet.rds')
COHORT <- file.path(ROOT,'data','DTS_PregnancyCohort_PrimaryCauseSpecific.rds')
FULL <- file.path(ROOT,'data','DTS_FullPredictionRiskSet_BinaryComparator.rds')
for (d in c('models','tables','qc','figures','prepared')) dir.create(file.path(OUT,d),recursive=TRUE,showWarnings=FALSE)
HEAT <- data.frame(metric=c('WBGT','UTCI','MRT','Heat Index'),key=c('WBGT','UTCI','MRT','HI'),
 raw=c('WBGT_T1_cummean','utci_T1_cummean','mrt_T1_cummean','HI_T1_cummean'),
 iqr=c('WBGT_T1_country_iqr','utci_T1_country_iqr','mrt_T1_country_iqr','HI_T1_country_iqr'))
# Set in the RStudio console AFTER sourcing settings, or as environment variables.
HEAT_RUN <- Sys.getenv('DTS_HEAT','WBGT')
OPTIMIZER <- Sys.getenv('DTS_OPTIMIZER','nlminb') # explicit BFGS option; no automatic retry
# Five-year stage always uses the stored outcome_last5yr == 1 flag.
COUNTRY_RUN <- Sys.getenv('DTS_COUNTRY','ALL')
MIXTURE_ENGINE <- Sys.getenv('DTS_MIXTURE_ENGINE','package') # package or mixed
GH_N <- as.integer(Sys.getenv('DTS_GH_N','50'))
CHUNK_PREG <- 3000L
SOC <- c('age_c','I(age_c^2)','rc_empl','v106','v190','rc_media_none','rc_hins_any','v025')
TEMP <- c('preg_start_year_f','preg_start_month_f')
AIR <- paste0(c('CO','DUST','PM','SO','SU'),'_T1_cummean')
primary_path <- function(metric,m) {
 label <- if(metric=='Heat Index') 'Heat_Index' else metric
 suffix <- if(m==1) '' else '_NoCountryMonthInt'
 if(metric=='Heat Index' && m==3) suffix <- '_BFGS'
 file.path(ROOT,'models','Linear_Nonlinear',paste0('GM',m,'_',label,'_DTS_IQR_CountryFE_MotherRE',suffix,'.rds'))
}
