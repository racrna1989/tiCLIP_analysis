#!/usr/bin/env Rscript

suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(fuzzyjoin))
suppressMessages(library(stringr))
args = commandArgs(trailingOnly=TRUE)

ANNOTATION_BED_FILEPATH<-args[1]
EXPRESSION_VECTOR_FILEPATH<-args[2]
IN_FILE_PATH<-args[3]
HOST_GENE_INDEX_FILEPATH<-args[4]
OUTFILE_PATH<-args[5]

#Function to process count files -------------------------
wrangle_bed_counts <- function(dataframe) {
  dataframe %>%
    #the gene ID is complicated and has different number of columns for some snoRNAs, best to label by grepl in new column.
    mutate(snoRNALocation = case_when(
      grepl(":::intronic:::", geneID) ~ "intronic", 
      grepl(":::nonintronic:::", geneID) ~ "non_intronic"
    )) %>% 
    mutate(region = case_when(
      grepl("intron3ss", geneID) ~ "intron3ss",
      grepl("intron5ss", geneID) ~ "intron5ss",
      grepl("mature3end", geneID) ~ "mature3end",
      grepl("mature5end", geneID) ~ "mature5end",
    )) %>%
      mutate(snoRNAType = case_when(
      grepl("SNORD|snoU|U3|U8|snoMe28S-Am2634|snoMBII|snoZ|snosnR66", geneID) ~ "cdBox", 
      grepl("SNORA|ACA|RNU105C|RNU105B", geneID) ~ "HACA",
      grepl("SCARNA", geneID) ~ "scaRNA",
      TRUE ~ "other"
    )) %>%
    mutate(geneID = as.character(geneID)) %>%
    #fuzzy_left join allows string matching rather than exact matching. Noticed more output afterwards, but unique fixed this. 
    #fuzzy_left_join(host_gene_index, by =c("geneID" = "geneID"), match_fun = str_detect) %>%
    #separate(hostGeneID, into=c("hostGeneID","transcriptBiotype","intronNumber","distFromTSS"),sep = ":::") %>% 
    #unique() %>% 
    #left join the expression vector using specific predefined columns. 
    #left_join(expression_vector, by = c("hostGeneID" = "geneID")) %>%
    #mutate(norm_count = count/ctrl_RNAseq_expr) %>%
    #select(-count, -ctrl_RNAseq_expr) %>%
    group_by(Sample,region,snoRNAType, snoRNALocation, DistToLandmark) %>% 
    summarise(sum_count = sum(count)) %>% 
    left_join(number_of_intron_annotations) %>% 
    mutate(sum_count_norm_annotation_number = sum_count/intron_count)
}




#Load expression vector -----------------------------------
print("load expression vector")

load(EXPRESSION_VECTOR_FILEPATH)

expression_vector<-left_join( 
  (as.data.frame(ctrl_RNAseq_expr) %>% 
     tibble::rownames_to_column(var = "geneID")),
  (as.data.frame(ctrl_TTseq_expr) %>% 
     tibble::rownames_to_column(var = "geneID"))
) %>% 
  select(-ctrl_TTseq_expr) %>% 
  mutate(ctrl_RNAseq_expr = case_when(
    ctrl_RNAseq_expr ==0 ~ min(ctrl_RNAseq_expr[ctrl_RNAseq_expr > 0]), 
    TRUE ~  ctrl_RNAseq_expr
  ))


#load annobed -----------------------
print("load annotation bed")

annoBed<-read.table(ANNOTATION_BED_FILEPATH, sep = "\t", header = F)
print(ncol(annoBed))
colnames(annoBed)<-c("chr", "start", "end", "geneID", "score", "strand")

number_of_intron_annotations<-
  annoBed %>% 
  mutate(snoRNALocation = case_when(
    grepl(":::intronic", geneID) ~ "intronic", 
    grepl(":::nonintronic", geneID) ~ "non_intronic"
  )) %>% 
  mutate(snoRNAType = case_when(
      grepl("SNORD|snoU|U3|U8|snoMe28S-Am2634|snoMBII|snoZ|snosnR66", geneID) ~ "cdBox", 
      grepl("SNORA|ACA|RNU105C|RNU105B", geneID) ~ "HACA",
      grepl("SCARNA", geneID) ~ "scaRNA",
      TRUE ~ "other"
    )) %>% 
  group_by(snoRNALocation, snoRNAType) %>%
  summarise(intron_count =n())

#load index file ------------------
#print("load index file")

#host_gene_index<-read.table(HOST_GENE_INDEX_FILEPATH, sep = "\t", header = F)
#colnames(host_gene_index)<-c("geneID", "hostGeneID")
#host_gene_index$geneID<-as.character(host_gene_index$geneID)


#load count file -----------------
print("load count file")

int_exon_junction_coverage<-read.table(IN_FILE_PATH, sep = "\t", header = F)
colnames(int_exon_junction_coverage)<-c("Sample","chr", "start", "end", "geneID", "DistToLandmark", "strand", "count")

int_exon_junction_coverage<-
  wrangle_bed_counts(int_exon_junction_coverage) 

#rbind the files together and write to output filename

OUTFILE<-int_exon_junction_coverage

write.table(OUTFILE, OUTFILE_PATH, quote = F, sep = "\t", row.names = F)



