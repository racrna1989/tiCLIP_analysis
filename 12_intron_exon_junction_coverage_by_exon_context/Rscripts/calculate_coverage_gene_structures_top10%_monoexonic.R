#!/usr/bin/env Rscript

suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
args = commandArgs(trailingOnly=TRUE)

ANNOTATION_BED_FILEPATH<-args[1]
EXPRESSION_VECTOR_FILEPATH<-args[2]
IN_FILE_PATH_3_END<-args[3]
IN_FILE_PATH_5_END<-args[4]
OUTFILE_PATH<-args[5]



#Load expression vector -----------------------------------

load(EXPRESSION_VECTOR_FILEPATH)

expression_vector<-left_join( 
  (as.data.frame(ctrl_RNAseq_expr) %>% 
     add_rownames(var = "geneID")),
  (as.data.frame(ctrl_TTseq_expr) %>% 
     add_rownames(var = "geneID"))
) %>% 
  select(-ctrl_TTseq_expr) %>% 
  mutate(ctrl_RNAseq_expr = case_when(
    ctrl_RNAseq_expr ==0 ~ min(ctrl_RNAseq_expr[ctrl_RNAseq_expr > 0]), 
    TRUE ~  ctrl_RNAseq_expr
  ))


#load annobed -----------------------

annoBed<-read.table(ANNOTATION_BED_FILEPATH, sep = "\t", header = F)
colnames(annoBed)<-c("chr", "start", "end", "geneID", "score", "strand")


#calculate top 10 of expressed genes ------------

top_10pct_expressed<-
  annoBed %>% 
  separate(geneID, into = c("geneID", "Biotype", "ExonNumber", "TotalNumberOfExon", "ExonSize", "ExonicDisance", "ExonDistFromTSS","ExonStature", "GeneStructure"), sep=":::") %>% 
  select(Biotype, geneID, GeneStructure) %>% 
  filter(grepl("single", GeneStructure)) %>% 
  select(geneID) %>%
  unique() %>% 
  left_join(expression_vector) %>%
  arrange(desc(ctrl_RNAseq_expr)) %>% 
  top_frac(0.1)

#calculate number of annotations (must have same filter as wrangle bed counts) -------

number_of_intron_annotations<-
  annoBed %>% 
  separate(geneID, into = c("geneID", "Biotype", "IntronNumber", "TotalNumberOfIntrons", "ExonSize", "ExonicDisance", "ExonDistFromTSS","ExonStature", "GeneStructure"), sep=":::") %>%
  filter(geneID %in% top_10pct_expressed$geneID) %>%
  group_by(GeneStructure) %>% 
  summarise(intron_count =n())


#Function to process count files -------------------------

wrangle_bed_counts <- function(dataframe) {
  dataframe %>%
    separate(geneID, into = c("geneID", "Biotype", "IntronNumber", "TotalNumberOfIntrons", "ExonSize", "ExonicDisance", "ExonDistFromTSS","ExonStature", "GeneStructure"), sep=":::") %>%
    filter(geneID %in% top_10pct_expressed$geneID) %>%
    left_join(expression_vector) %>% 
    #this replaces NAs introduced by no value present in expression_vector, and replaces them with min value in expression_vector
    mutate_at(vars(ctrl_RNAseq_expr), ~replace(., is.na(.), min(expression_vector$ctrl_RNAseq_expr))) %>% 
    mutate(norm_count = count/ctrl_RNAseq_expr) %>%
    select(-count, -ctrl_RNAseq_expr) %>%
    group_by(Sample,GeneStructure, DistToLandmark) %>% 
    summarise(sum_RNAseq_norm_count_norm_annotation_number = sum(norm_count)) %>% 
    left_join(number_of_intron_annotations) %>% 
    mutate(sum_RNAseq_norm_count_norm_annotation_number = sum_RNAseq_norm_count_norm_annotation_number/intron_count)
}




#load 5end count file -----------------

int_exon_junction_coverage_five_end<-read.table(IN_FILE_PATH_5_END, sep = "\t", header = F)
colnames(int_exon_junction_coverage_five_end)<-c("Sample","chr", "start", "end", "geneID", "DistToLandmark", "strand", "count")

int_exon_junction_coverage_five_end<-
  int_exon_junction_coverage_five_end %>% 
  wrangle_bed_counts() %>%
  mutate(Position = "five_end")

#load 3end count file -----------------

int_exon_junction_coverage_three_end<-read.table(IN_FILE_PATH_3_END, sep = "\t", header = F)
colnames(int_exon_junction_coverage_three_end)<-c("Sample","chr", "start", "end", "geneID", "DistToLandmark", "strand", "count")

int_exon_junction_coverage_three_end <-
  int_exon_junction_coverage_three_end %>% 
  wrangle_bed_counts() %>%
  mutate(Position = "three_end")

#rbind the files together and write to output filename

OUTFILE<-rbind(int_exon_junction_coverage_three_end,int_exon_junction_coverage_five_end)

write.table(OUTFILE, OUTFILE_PATH, quote = F, sep = "\t", row.names = F)


