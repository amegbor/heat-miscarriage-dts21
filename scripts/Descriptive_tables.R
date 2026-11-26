# Run in OOD RStudio. No model fitting; original tables remain unchanged.
suppressPackageStartupMessages({library(data.table);library(glmmTMB)})
setDTthreads(1L)
# Convert labelled DHS vectors by their underlying values, without changing factors.
plain_values <- function(x) {
 if(inherits(x,c('haven_labelled','haven_labelled_spss','labelled'))) attributes(x)<-NULL
 x
}
safe_chr <- function(x) as.character(plain_values(x))
root <- Sys.getenv('DTS_ROOT',file.path(getwd(),'output'))
out <- file.path(root,'Supplementary_DTS_v1','tables','DTS_Descriptive_v2')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
put <- function(z,name) fwrite(z,file.path(out,paste0(name,'.csv')))
need <- function(d,v) if(length(setdiff(v,names(d)))) stop('Missing variables: ',paste(setdiff(v,names(d)),collapse=', '))
finite <- function(d,v) {
 need(d,v);ok<-complete.cases(d[,..v])
 for(x in v) if(is.numeric(d[[x]])) ok<-ok & is.finite(d[[x]])
 ok
}
country_names <- c(AO='Angola',BD='Bangladesh',BF='Burkina Faso',CD='Democratic Republic of the Congo',
 CI="Cote d'Ivoire",ET='Ethiopia',GH='Ghana',JO='Jordan',KE='Kenya',KH='Cambodia',LS='Lesotho',
 ML='Mali',MW='Malawi',MZ='Mozambique',NG='Nigeria',NP='Nepal',PH='Philippines',SN='Senegal',
 TJ='Tajikistan',TZ='Tanzania',ZM='Zambia')
label_country <- function(z) {
 code<-toupper(sub('[0-9]+$','',safe_chr(z)))
 ans<-unname(country_names[code]);ans[is.na(ans)]<-safe_chr(z)[is.na(ans)];ans
}
cohort <- as.data.table(readRDS(file.path(root,'data','DTS_PregnancyCohort_PrimaryCauseSpecific.rds')))
d <- as.data.table(readRDS(file.path(root,'data','DTS_ObservedRiskSet.rds')))
need(cohort,c('preg_uid','country','outcome_type'))
need(d,c('preg_uid','mother_uid','cluster_uid','country','gest_month_f','miscarriage_event'))
if(anyDuplicated(cohort$preg_uid)) stop('Pregnancy cohort must contain one row per pregnancy')
if(anyNA(cohort$preg_uid)||anyNA(d$preg_uid)) stop('Missing pregnancy IDs')
if(anyDuplicated(d[,.(preg_uid,gest_month_f)])) stop('Duplicate observed pregnancy/interval keys')
if(length(setdiff(unique(d$preg_uid),cohort$preg_uid))) stop('Observed pregnancies not found in cohort')

rows<-flow<-baseline<-exposure_stats<-list(); manifest<-list()
summarize_exposures <- function(z,vars,sample) rbindlist(lapply(vars,function(v) {
 x<-plain_values(z[[v]]);good<-is.finite(x);x<-x[good]
 data.table(sample,variable=v,unit_of_description='Retained person-months',
  person_months=nrow(z),available=sum(good),missing=sum(!good),
  mean=if(length(x)) mean(x) else NA_real_,SD=if(length(x)>1) sd(x) else NA_real_,
  median=if(length(x)) median(x) else NA_real_,
  Q1=if(length(x)) unname(quantile(x,.25)) else NA_real_,
  Q3=if(length(x)) unname(quantile(x,.75)) else NA_real_,
  minimum=if(length(x)) min(x) else NA_real_,maximum=if(length(x)) max(x) else NA_real_)
}))
catvars<-c('rc_empl','v106','v190','rc_media_none','rc_hins_any','v025')
baseline_vars<-c('age_at_preg_start',catvars)
baseline_vars<-intersect(baseline_vars,names(cohort))
# Optional old FINAL analytic pregnancy dataset, not the unfiltered NR_Data file.
oldpath<-Sys.getenv('DTS_OLD_ANALYTIC_COHORT','')
old<-if(nzchar(oldpath)) as.data.table(readRDS(oldpath)) else NULL
if(!is.null(old)) {
 need(old,c('preg_uid','country','outcome_type',baseline_vars))
 if(anyNA(old$preg_uid)||anyDuplicated(old$preg_uid)) stop('Old cohort must have one identified row per analytic pregnancy')
}
same_baseline <- function(p) {
 if(is.null(old)) return(NA)
 if(nrow(p)!=nrow(old)||!setequal(safe_chr(p$preg_uid),safe_chr(old$preg_uid))) return(FALSE)
 a<-p[order(safe_chr(preg_uid))];b<-old[order(safe_chr(preg_uid))]
 all(vapply(c('country','outcome_type',baseline_vars),function(v)
 identical(safe_chr(a[[v]]),safe_chr(b[[v]])),logical(1)))
}
country_table <- function(p,z,sample) {
 # Pregnancy outcomes describe actual outcomes even if their event interval was excluded.
 a<-p[,.(pregnancies=.N,
  live_birth=sum(outcome_type=='live_birth'),miscarriage=sum(outcome_type=='miscarriage'),
  stillbirth=sum(outcome_type=='stillbirth'),abortion=sum(outcome_type=='abortion'),
  other_outcome=sum(!outcome_type %chin% c('live_birth','miscarriage','stillbirth','abortion'))),by=country]
 b<-z[,.(person_months=.N,mothers=uniqueN(mother_uid),clusters=uniqueN(cluster_uid),
  retained_miscarriage_events=sum(miscarriage_event)),by=country]
 a<-merge(a,b,by='country',all.x=TRUE)
 a[,`:=`(sample=sample,country_name=label_country(country),
  pregnancy_share_percent=100*pregnancies/sum(pregnancies),
  miscarriage_outcome_percent=100*miscarriage/pregnancies)]
 for(v in c('live_birth','miscarriage','stillbirth','abortion')) a[,(paste0(v,'_percent')):=100*get(v)/pregnancies]
 a
}
summarize_baseline <- function(p,sample) {
 numeric<-intersect('age_at_preg_start',baseline_vars)
 a<-rbindlist(lapply(numeric,function(v) p[,{
  x<-plain_values(get(v));good<-is.finite(x)
  list(variable=v,category=NA_character_,denominator=.N,n=sum(good),missing=sum(!good),
   mean=if(any(good)) mean(x[good]) else NA_real_,SD=if(sum(good)>1) sd(x[good]) else NA_real_,
   median=if(any(good)) median(x[good]) else NA_real_,
   Q1=if(any(good)) unname(quantile(x[good],.25)) else NA_real_,
   Q3=if(any(good)) unname(quantile(x[good],.75)) else NA_real_)
 },by=outcome_type]),fill=TRUE)
 b<-rbindlist(lapply(intersect(catvars,baseline_vars),function(v) {
  ans<-p[,.(n=.N),by=.(outcome_type,category=safe_chr(get(v)))]
  ans[,`:=`(variable=v,denominator=sum(n),percent=100*n/sum(n)),by=outcome_type];ans
 }),fill=TRUE)
 ans<-rbindlist(list(a,b),fill=TRUE);ans[,sample:=sample];ans
}
add_sample <- function(z,sample,source,check_old=FALSE) {
 p<-cohort[safe_chr(preg_uid) %chin% safe_chr(unique(z$preg_uid))]
 rows[[sample]]<<-country_table(p,z,sample)
 per_preg<-z[,.(intervals=.N,events=sum(miscarriage_event)),by=preg_uid]
 flow[[sample]]<<-data.table(sample,source,person_months=nrow(z),pregnancies=nrow(p),
  mothers=uniqueN(z$mother_uid),clusters=uniqueN(z$cluster_uid),events=sum(z$miscarriage_event),
  miscarriage_outcome_pregnancies=sum(p$outcome_type=='miscarriage'),
  retained_event_pregnancies=sum(per_preg$events>0),median_retained_intervals=median(per_preg$intervals))
 if(check_old) {
  equal<-same_baseline(p)
  manifest[[sample]]<<-data.table(sample,unchanged_from_old=equal,
   recommendation=if(isTRUE(equal)) 'Omit repeated pregnancy baseline descriptions' else
    if(is.na(equal)) 'Candidate only: old final analytic cohort not supplied; compare before including' else 'Include updated baseline descriptions; analytic population/values differ')
  if(!isTRUE(equal)) baseline[[sample]]<<-summarize_baseline(p,sample)
 }
}
add_sample(d,'Observed DTS risk set before model-specific complete cases',file.path(root,'data','DTS_ObservedRiskSet.rds'))
# Cohort outcome counts are separate from the analytic risk-set counts.
cohort_counts<-cohort[,.(pregnancies=.N),by=.(country,outcome_type)]
cohort_counts[,country_name:=label_country(country)]
put(cohort_counts,'S1_eligible_cohort_country_outcomes_REFERENCE')

model_files<-c(WBGT='GM3_WBGT_DTS_IQR_CountryFE_MotherRE_NoCountryMonthInt.rds',
 UTCI='GM3_UTCI_DTS_IQR_CountryFE_MotherRE_NoCountryMonthInt.rds',
 MRT='GM3_MRT_DTS_IQR_CountryFE_MotherRE_NoCountryMonthInt.rds',
 HI='GM3_Heat_Index_DTS_IQR_CountryFE_MotherRE_BFGS.rds')
for(metric in names(model_files)) {
 path<-file.path(root,'models','Linear_Nonlinear',model_files[[metric]])
 if(!file.exists(path)) stop('Missing primary GM3 model: ',path)
 fit<-readRDS(path)
 if(!isTRUE(fit$fit$convergence==0)||!isTRUE(fit$sdr$pdHess)) stop('Invalid primary fit: ',metric)
 v<-all.vars(formula(fit));need(d,v)
 z<-d[finite(d,v)]
 if(nrow(z)!=nobs(fit)) stop('Reconstructed GM3 N differs: ',metric)
 mf<-fit$frame
 if(is.null(mf)) stop('Missing saved model frame: ',metric)
 for(x in intersect(v,names(mf))) if(!identical(safe_chr(z[[x]]),safe_chr(mf[[x]])))
  stop('GM3 model-frame row/value mismatch: ',metric,' / ',x)
 sample_label<-paste(metric,'primary GM3 complete-case')
 add_sample(z,sample_label,path,TRUE)
 stem<-switch(metric,UTCI='utci',MRT='mrt',metric)
 description_vars<-intersect(c(paste0(stem,'_T1_cummean'),paste0(stem,'_T1_country_iqr')),names(z))
 exposure_stats[[metric]]<-summarize_exposures(z,description_vars,sample_label)
 rm(fit,z,mf);gc()
}

# Reconstruct the raw common mixture sample; the saved prepared data are quantized.
heat<-c('WBGT_T1_cummean','utci_T1_cummean','mrt_T1_cummean','HI_T1_cummean')
air<-paste0(c('CO','DUST','PM','SO','SU'),'_T1_cummean')
exposures<-c(heat,air)
cov<-c('gest_month_f','country','preg_start_year_f','preg_start_month_f',
 'age_at_preg_start',catvars,'miscarriage_event','preg_uid','mother_uid','cluster_uid')
z<-d[finite(d,c(cov,exposures))]
prepared<-file.path(root,'PIF_QGComp_v3','prepared','Mixture_common_quartiles.rds')
if(!file.exists(prepared)) stop('Missing saved mixture sample')
q<-as.data.table(readRDS(prepared));need(q,c('preg_uid','gest_month_f'))
if(nrow(z)!=nrow(q) || !identical(safe_chr(z$preg_uid),safe_chr(q$preg_uid)) ||
 !identical(safe_chr(z$gest_month_f),safe_chr(q$gest_month_f))) stop('Mixture sample row/key mismatch')
add_sample(z,'Common complete-case mixture sample',prepared,TRUE)
exposure_stats[['mixture']]<-summarize_exposures(z,exposures,'Common complete-case mixture sample')
put(rbindlist(exposure_stats,fill=TRUE),'DTS_UPDATED_person_month_exposure_descriptives')

# S2: Spearman correlations of raw cumulative exposures on retained person-months.
correlation<-cor(as.matrix(z[,..exposures]),method='spearman')
labels<-c('WBGT','UTCI','MRT','Heat Index','Carbon monoxide','Mineral dust','PM2.5','Sulphur dioxide','Sulphate')
dimnames(correlation)<-list(labels,labels)
oldcorpath<-Sys.getenv('DTS_OLD_CORRELATION_CSV','')
same_cor<-NA
if(nzchar(oldcorpath)) {
 oc<-fread(oldcorpath,check.names=FALSE)
 need(oc,c('component',labels))
 if(!setequal(oc$component,labels)) stop('Old correlation CSV needs matching component labels')
 om<-as.matrix(oc[match(labels,component),..labels])
 same_cor<-isTRUE(all.equal(unname(correlation),unname(om),tolerance=1e-12,check.attributes=FALSE))
}
if(!isTRUE(same_cor)) put(data.table(component=labels,as.data.table(correlation)),
 if(is.na(same_cor)) 'S2_DTS_Spearman_correlations_COMPARE_WITH_OLD' else 'S2_DTS_Spearman_correlations_UPDATED')
put(data.table(sample='Common complete-case mixture sample',person_months=nrow(z),
 pregnancies=uniqueN(z$preg_uid),method='Spearman; raw cumulative first-trimester exposures; unweighted retained person-months',
 unchanged_from_old=same_cor,recommendation=if(isTRUE(same_cor)) 'Reuse old correlation table; update denominator/caption if needed' else
 if(is.na(same_cor)) 'Compare against old correlation table before final inclusion' else 'Use updated S2 matrix'),'S2_population_and_comparison')
put(rbindlist(rows,fill=TRUE),'S1_DTS_country_outcomes_and_analytic_denominators')
put(rbindlist(flow,fill=TRUE),'DTS_analytic_sample_flow')
put(rbindlist(manifest,fill=TRUE),'Baseline_comparison_and_inclusion_decisions')
if(length(baseline)) {
 b<-rbindlist(baseline,fill=TRUE)
 m<-rbindlist(manifest,fill=TRUE)
 changed<-m[unchanged_from_old==FALSE,sample]
 unknown<-m[is.na(unchanged_from_old),sample]
 if(length(changed)) put(b[sample %chin% changed],'DTS_UPDATED_baseline_descriptives_by_outcome')
 if(length(unknown)) put(b[sample %chin% unknown],'DTS_baseline_descriptives_CANDIDATE_compare_before_including')
}
writeLines(c('Unweighted descriptions. No model refitting.',
 'S1 counts actual pregnancy outcomes separately from retained miscarriage event intervals.',
 'Stillbirth, abortion and live birth are censoring outcomes in the cause-specific DTS workflow.',
 'The window ends after interval 6; later outcomes are descriptive outcomes, not events within that window.',
 'Pregnancy baseline summaries use each retained pregnancy once; person-month denominators describe retained intervals.',
 'S2 uses raw cumulative exposure values from the common complete-case person-month sample, not quantized categories.',
 'No independence-based correlation p-values are reported for repeated pregnancy intervals.',
 'The eligible-cohort reference table is for reconciliation; do not repeat it if the earlier cohort table is unchanged.',
 'Pending five-year sensitivity and full prediction/PIF cohorts are not combined with the primary analytic sample.',
 'Candidate baseline summaries are not automatically selected for the manuscript.'),file.path(out,'Table_population_notes.txt'))
print(rbindlist(flow,fill=TRUE));print(rbindlist(manifest,fill=TRUE))
message('Saved descriptive outputs: ',out)
