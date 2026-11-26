#!/usr/bin/env Rscript
# Run one explicit stage in a fresh process; never submits or fits other stages.
a <- commandArgs(trailingOnly=TRUE)
if (!length(a)) stop("Usage: Rscript scripts/run.R STAGE [WBGT|UTCI|MRT|HI]")
script_arg <- grep("^--file=", commandArgs(FALSE), value=TRUE)
if (!length(script_arg)) stop("Use Rscript scripts/run.R rather than source() for this runner")
bundle <- dirname(normalizePath(sub("^--file=", "", script_arg[1]), mustWork=TRUE))
repo <- dirname(bundle)
if (file.exists(file.path(repo,"config.R"))) source(file.path(repo,"config.R"))
root <- Sys.getenv("DTS_ROOT")
if (!nzchar(root)) stop("Set DTS_ROOT or copy config.example.R to config.R and edit it")
root <- normalizePath(root,mustWork=TRUE)
Sys.setenv(DTS_BUNDLE=bundle,DTS_ROOT=root,DTS_SUPP_OUT=file.path(root,"Supplementary_DTS_v1"))
stage <- a[1]
metric <- if (length(a)>1) a[2] else "WBGT"
if (metric=="HI") metric <- "Heat Index"
if (!metric %in% c("WBGT","UTCI","MRT","Heat Index")) stop("Unknown heat metric")
Sys.setenv(DTS_HEAT=metric,DTS_COUNTRY="ALL")
if (!nzchar(Sys.getenv("PIF_GH_N"))) Sys.setenv(PIF_GH_N="50")
if (!nzchar(Sys.getenv("DTS_GH_N"))) Sys.setenv(DTS_GH_N="50")
# Initial supplementary fits used the metric-specific optimizer choice.
if (!nzchar(Sys.getenv("DTS_OPTIMIZER"))) Sys.setenv(DTS_OPTIMIZER=if(metric=="Heat Index")"BFGS" else "nlminb")
s <- function(x) source(file.path(bundle,x))
settings <- function() {s("supp/00_settings.R");s("supp/01_helpers.R")}
allowed <- c("primary","primary_recovery","raw","nonlinear","country","five_year",
 "five_year_utci_recovery","s3","prepare","pif","mixture_package","mixture_RR_package",
 "mixture_RR_mixed","country_interaction","collect_interaction","descriptive","collect_s4","figures","map_flow")
if (!stage %in% allowed) stop("Unknown stage: ",stage)
if(stage=="primary") {s("Primary_linear.R")} else if(stage=="primary_recovery") {s("Recover_primary21.R")} else if(stage=="five_year_utci_recovery") {s("Recover_five_year_UTCI21.R")} else if(stage %in% c("raw","nonlinear","country","five_year","s3","mixture_RR_package","mixture_RR_mixed")) {
 settings()
 if(stage %in% c("mixture_RR_package","mixture_RR_mixed")) {
  MIXTURE_ENGINE <- if(stage=="mixture_RR_package")"package" else "mixed"
  FULL <- file.path(root,"data","DTS_FullPredictionRiskSet_Mixture.rds")
  s("supp/08_mixture_net_RR.R")
 } else s(file.path("supp",switch(stage,raw="02_raw_linear.R",nonlinear="03_nonlinear.R",
   country="06_country_specific.R",five_year="05_five_year.R",s3="04_primary_S3.R")))
} else if(stage=="prepare") {
 if (!file.exists(file.path(root,"data","DTS_FullPredictionRiskSet_Mixture.rds")))
  stop("Supply the preconstructed mixture prediction risk set; raw linkage is outside this release")
 s("pif/02_prepare.R")
} else if(stage=="pif") {s("pif/03_PIF.R")} else if(stage=="mixture_package") {
 key <- c(WBGT="WBGT",UTCI="UTCI",MRT="MRT","Heat Index"="HI")[[metric]]
 Sys.setenv(DTS_MIXTURE=paste0("Air_",key),DTS_MIXTURE_ENGINE="package")
 s("pif/04_QGComp.R")
} else if(stage %in% c("country_interaction","collect_interaction")) {
 Sys.setenv(DTS_COUNTRY_TEST_ACTION=if(stage=="collect_interaction")"collect" else "fit")
 s("Country_interaction21.R")
} else if(stage=="descriptive") {s("Descriptive_tables.R");s("Table1.R")} else if(stage=="collect_s4") {s("Collect_S4.R")} else if(stage=="figures") {s("Manuscript_figures.R")} else if(stage=="map_flow") {s("Study_map_and_flow.R")}
message("Completed: ",stage," / ",metric)
