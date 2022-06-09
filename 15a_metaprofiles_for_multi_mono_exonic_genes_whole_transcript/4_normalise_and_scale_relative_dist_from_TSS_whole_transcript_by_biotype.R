#!/usr/bin/env Rscript

#get current wd
WD<-getwd()

#source scripts
source(paste0(WD,"/master_scripts/15a_metaprofiles_for_multi_mono_exonic_genes_whole_transcript/process_dist_to_TSS_1_whole_transcript.R"))
source(paste0(WD,"/master_scripts/15a_metaprofiles_for_multi_mono_exonic_genes_whole_transcript/scale_and_normalise_2_3_whole_transcript_by_biotype.R"))

#load libraries
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))

#define arguments
args = commandArgs(trailingOnly=TRUE)

DATA_FRAME_FILEPATH<-args[1]
ANNOTATION_BED_FILEPATH<-args[2]
rRNA_FACTOR_FILEPATH<-args[3]
RPM_FACTOR_FILEPATH<-args[4]
EXPRESSION_VECTOR<-args[6]
OUTFILE_PATH<-args[5]

#run scripts on dataframes

DF_1<-process_distToTSS_1(DATA_FRAME_FILEPATH,ANNOTATION_BED_FILEPATH)
  
DF_2<-ScaledNormalise_2_3(DF_1, 40, 10000000, ANNOTATION_BED_FILEPATH, rRNA_FACTOR_FILEPATH,RPM_FACTOR_FILEPATH, EXPRESSION_VECTOR)

#write table
write.table(DF_2, OUTFILE_PATH, quote = F, sep = "\t", row.names = F)
