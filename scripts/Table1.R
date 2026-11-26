# Revised manuscript Table 1. Run in OOD RStudio; no refitting.
suppressPackageStartupMessages({library(data.table);library(glmmTMB)})
setDTthreads(1L)
plain <- function(x) {
 if(inherits(x,c('haven_labelled','haven_labelled_spss','labelled'))) attributes(x)<-NULL
 x
}
chr <- function(x) as.character(plain(x))
root <- Sys.getenv('DTS_ROOT',file.path(getwd(),'output'))
out <- file.path(root,'Supplementary_DTS_v1','tables','DTS_Descriptive_v2')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
need <- function(d,v) if(length(setdiff(v,names(d)))) stop('Missing: ',paste(setdiff(v,names(d)),collapse=', '))
cohort <- as.data.table(readRDS(file.path(root,'data','DTS_PregnancyCohort_PrimaryCauseSpecific.rds')))
d <- as.data.table(readRDS(file.path(root,'data','DTS_ObservedRiskSet.rds')))
model_path <- file.path(root,'models','Linear_Nonlinear','GM3_UTCI_DTS_IQR_CountryFE_MotherRE_NoCountryMonthInt.rds')
fit <- readRDS(model_path)
if(!isTRUE(fit$fit$convergence==0)||!isTRUE(fit$sdr$pdHess)) stop('Primary model failed checks')
vars<-all.vars(formula(fit));need(d,vars)
ok<-complete.cases(d[,..vars])
for(v in vars) if(is.numeric(d[[v]])) ok<-ok & is.finite(d[[v]])
analytic<-d[ok]
if(nrow(analytic)!=nobs(fit)) stop('Analytic row count differs from saved UTCI GM3 model')
for(v in intersect(vars,names(fit$frame))) if(!identical(chr(analytic[[v]]),chr(fit$frame[[v]]))) stop('Model-frame mismatch: ',v)
ids<-chr(unique(analytic$preg_uid))
need(cohort,c('preg_uid','outcome_type','age_at_preg_start','rc_empl','v106','v190','rc_media_none','rc_hins_any','v025'))
if(anyDuplicated(cohort$preg_uid)) stop('Cohort has duplicate pregnancy IDs')
p<-copy(cohort[chr(preg_uid) %chin% ids])
if(nrow(p)!=length(ids)) stop('Pregnancy cohort linkage failed')
if(anyNA(p$outcome_type)||any(!chr(p$outcome_type) %chin% c('live_birth','miscarriage','stillbirth','abortion'))) stop('Unexpected outcome categories')

# Each pregnancy contributes one exposure value, at its FINAL observed interval.
# Early losses therefore use exposure accumulated only up to that interval.
# Completed first-trimester means remain frozen after gestational month 3.
exposures<-c('WBGT_T1_cummean','utci_T1_cummean','mrt_T1_cummean','HI_T1_cummean',
 'CO_T1_cummean','DUST_T1_cummean','PM_T1_cummean','SO_T1_cummean','SU_T1_cummean')
need(d,c('preg_uid','gest_month',exposures))
terminal<-d[chr(preg_uid) %chin% ids,c('preg_uid','gest_month',exposures),with=FALSE]
terminal[,gest_month:=as.integer(chr(gest_month))]
if(anyNA(terminal$gest_month)||anyDuplicated(terminal[,.(preg_uid,gest_month)])) stop('Invalid interval keys')
setorder(terminal,preg_uid,gest_month)
terminal<-terminal[, .SD[.N],by=preg_uid]
idx<-match(chr(p$preg_uid),chr(terminal$preg_uid))
if(anyNA(idx)) stop('Missing terminal exposure row')
for(v in exposures) set(p,j=v,value=plain(terminal[[v]][idx]))
rm(d,analytic,fit,terminal);gc()

groups<-c(Overall='ALL','Live birth'='live_birth',Miscarriage='miscarriage',Stillbirth='stillbirth',Abortion='abortion')
indices<-lapply(groups,function(g) if(g=='ALL') seq_len(nrow(p)) else which(chr(p$outcome_type)==g))
sizes<-vapply(indices,length,integer(1))
headers<-paste0(names(groups),' (N = ',format(sizes,big.mark=',',trim=TRUE),')')
table_rows<-list();availability<-list()
add <- function(label,values) {
 row<-as.list(c(label,values));names(row)<-c('Characteristic',headers)
 table_rows[[length(table_rows)+1L]]<<-as.data.table(row)
}
section <- function(label) add(label,rep('',length(groups)))
continuous <- function(v,label) {
 cells<-vapply(seq_along(indices),function(i) {
  x<-plain(p[[v]][indices[[i]]]);good<-is.finite(x)
  availability[[length(availability)+1L]]<<-data.table(variable=v,group=names(groups)[i],pregnancies=sizes[i],available=sum(good),missing=sum(!good))
  if(!any(good)) return('Not available')
  q<-quantile(x[good],c(.5,.25,.75),names=FALSE)
  cell<-sprintf('%.2f (%.2f, %.2f)',q[1],q[2],q[3])
  if(any(!good)) cell<-paste0(cell,' [n = ',format(sum(good),big.mark=',',trim=TRUE),']')
  cell
 },character(1))
 add(label,cells)
}
categorical <- function(v,title,codes,labels) {
 section(title)
 x<-chr(p[[v]])
 if(anyNA(x)||any(!x %chin% codes)) stop('Unexpected/missing category for ',v,'; review coding before creating table')
 for(k in seq_along(codes)) {
  cells<-vapply(seq_along(indices),function(i) {
   n<-sum(x[indices[[i]]]==codes[k]);sprintf('%s (%.1f%%)',format(n,big.mark=',',trim=TRUE),100*n/sizes[i])
  },character(1))
  add(labels[k],cells)
 }
}
section('Cumulative first-trimester heat exposure at the final observed interval')
for(i in 1:4) continuous(exposures[i],c('Wet-bulb globe temperature (WBGT), degrees Celsius',
 'Universal Thermal Climate Index (UTCI), degrees Celsius','Mean radiant temperature (MRT), degrees Celsius',
 'Heat Index (HI), degrees Celsius')[i])
section('Cumulative first-trimester ambient air-pollution exposure at the final observed interval')
for(i in 5:9) continuous(exposures[i],c('Carbon monoxide (CO)','Mineral dust','Fine particulate matter (PM2.5)',
 'Sulphur dioxide (SO2)','Sulphate (SO4)')[i-4])
section('Individual and household characteristics')
continuous('age_at_preg_start','Age at pregnancy initiation, years')
# Employment recode meanings are not supplied; keep the original numeric codes.
categorical('rc_empl','Employment status',c('0','1','2'),c('Code 0','Code 1','Code 2'))
categorical('v106','Educational attainment',as.character(0:3),c('No education','Primary','Secondary','Higher'))
categorical('v190','Household wealth quintile',as.character(1:5),c('Poorest','Poorer','Middle','Richer','Richest'))
categorical('rc_media_none','No media access',c('0','1'),c('No','Yes'))
categorical('rc_hins_any','Health insurance coverage',c('0','1'),c('No','Yes'))
categorical('v025','Place of residence',c('1','2'),c('Urban','Rural'))
z<-rbindlist(table_rows)
fwrite(z,file.path(out,'Table1_DTS_revised.csv'))
fwrite(rbindlist(availability),file.path(out,'Table1_DTS_exposure_and_age_available_counts.csv'))
title<-'Table 1. First-trimester heat and ambient air-pollution exposures and individual and household characteristics of pregnancies included in the discrete-time survival analyses, overall and by pregnancy outcome.'
notes<-c('Values are median (25th, 75th percentiles) or number (column percentage). All descriptions are unweighted; each pregnancy contributes once.',
 'The reference descriptive sample is the primary UTCI/MRT Model 3 sample. The WBGT/HI sample contains 717,281 pregnancies rather than 717,314; use S1 for metric-specific denominators.',
 'Exposures are cumulative means at each pregnancy\'s final observed interval, up to gestational month 3. The completed first-trimester mean is held fixed thereafter; pregnancies ending earlier use exposure accrued up to their observed ending interval.',
 'Bracketed n values identify non-missing exposure denominators. Age and household characteristics use the pregnancy denominator stated in the column heading.',
 'Outcome columns describe observed pregnancy outcomes. Miscarriage is the event of interest; other outcomes contribute follow-up until censoring or the end of the six-month risk window.',
 'Pregnancy-level summaries in this table differ from person-month exposure summaries and should not be labelled as person-month-weighted estimates.',
 'Employment labels remain codes pending confirmation of the study-specific recoding. Confirm pollution source units before adding units to their row labels; units have not been inferred.',
 'Mineral dust is distinct from total PM2.5. Cooking fuel is not included.')
md<-c(title,'',paste0('| ',paste(names(z),collapse=' | '),' |'),
 paste0('| ',paste(rep('---',ncol(z)),collapse=' | '),' |'),
 apply(as.matrix(z),1,function(r)paste0('| ',paste(r,collapse=' | '),' |')),'',paste0('Note. ',notes))
writeLines(md,file.path(out,'Table1_DTS_revised.md'))
writeLines(c(title,notes),file.path(out,'Table1_DTS_caption_and_notes.txt'))
print(z)
message('Saved revised Table 1 to: ',out)
