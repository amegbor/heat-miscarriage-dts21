#!/usr/bin/env Rscript
# Checks syntax and dependencies without loading study data or fitting models.
args <- grep("^--file=",commandArgs(FALSE),value=TRUE)
repo <- dirname(dirname(normalizePath(sub("^--file=","",args[1]))))
files <- list.files(file.path(repo,"scripts"),pattern="\\.R$",recursive=TRUE,full.names=TRUE)
for(f in files) {parse(file=f);cat("Parsed:",basename(f),"\n")}
required <- c("data.table","glmmTMB","dplyr","tibble","broom.mixed","readr","lme4",
 "digest","Matrix","qgcomp","sandwich","statmod","ggplot2","patchwork","maps")
installed <- vapply(required,requireNamespace,logical(1),quietly=TRUE)
print(data.frame(package=required,available=installed))
if(any(!installed)) stop("Missing dependencies: ",paste(required[!installed],collapse=", "))
cat("Syntax and dependency checks passed. No models were fitted.\n")
