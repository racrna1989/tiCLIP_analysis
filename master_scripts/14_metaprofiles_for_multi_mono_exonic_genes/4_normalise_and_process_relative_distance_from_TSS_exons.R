#!/usr/bin/env Rscript

#get current wd
WD<-getwd()

#source scripts
source(paste0(WD,"/master_scripts/14_metaprofiles_for_multi_mono_exonic_genes/process_dist_to_TSS_1.R"))
source(paste0(WD,"/master_scripts/14_metaprofiles_for_multi_mono_exonic_genes/scale_and_normalise_2_3.R"))

#load libraries
suppressMessages(library(scales))
suppressMessages(library(dplyr))
suppressMessages(library(tibble))
suppressMessages(library(tidyr))
suppressMessages(library(ggplot2))

#define arguments
args = commandArgs(trailingOnly=TRUE)

DATA_FRAME_FILEPATH<-args[1]
ANNOTATION_BED_FILEPATH<-args[2]
rRNA_FACTOR_FILEPATH<-args[3]
RPM_FACTOR_FILEPATH<-args[4]
OUTFILE_PATH<-args[5]

#run scripts on dataframes

DF_1<-process_distToTSS_1(DATA_FRAME_FILEPATH,ANNOTATION_BED_FILEPATH)
  
DF_2<-ScaledNormalise_2_3(DF_1, 200, 10000000, ANNOTATION_BED_FILEPATH, rRNA_FACTOR_FILEPATH,RPM_FACTOR_FILEPATH)

#write table
write.table(DF_2, OUTFILE_PATH, quote = F, sep = "\t", row.names = F)
