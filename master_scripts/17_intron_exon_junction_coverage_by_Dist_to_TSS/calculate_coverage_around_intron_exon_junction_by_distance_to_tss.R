#!/usr/bin/env Rscript

suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
args = commandArgs(trailingOnly=TRUE)

ANNOTATION_BED_FILEPATH<-args[1]
EXPRESSION_VECTOR_FILEPATH<-args[2]
IN_FILE_PATH_3_END<-args[3]
IN_FILE_PATH_5_END<-args[4]
OUTFILE_PATH<-args[5]


#Function to process count files -------------------------

#3end
wrangle_bed_counts_for_3end <- function(dataframe) {
  dataframe %>%
    separate(geneID, into = c("geneID", "Biotype", "ExonNumber", "TotalNumberOfExons", "ExonSize", "ExonicDisanceToTSS", "GenomicDistToTSS","ExonStature", "GeneStructure"), sep=":::") %>%
    left_join(expression_vector, by = "geneID") %>% 
    #this replaces NAs introduced by no value present in expression_vector, and replaces them with min value in expression_vector
  	mutate_at(vars(ctrl_RNAseq_expr), ~replace(., is.na(.), min(expression_vector$ctrl_RNAseq_expr))) %>% 
    mutate(norm_count = count/ctrl_RNAseq_expr) %>%
    select(-count, -ctrl_RNAseq_expr) %>%
    inner_join(DistGroup_3end, by = c("geneID", "ExonNumber")) %>% 
    group_by(Sample, DistToLandmark, DistGroup) %>% 
    summarise(sum_RNAseq_norm_count_norm_annotation_number = sum(norm_count)) %>% 
    left_join(number_of_intron_annotations_3end, by = "DistGroup") %>% 
    mutate(sum_RNAseq_norm_count_norm_annotation_number = sum_RNAseq_norm_count_norm_annotation_number/annotation_count)
}


#5end
wrangle_bed_counts_for_5end <- function(dataframe) {
  dataframe %>%
    separate(geneID, into = c("geneID", "Biotype", "ExonNumber", "TotalNumberOfExons", "ExonSize", "ExonicDisanceToTSS", "GenomicDistToTSS","ExonStature", "GeneStructure"), sep=":::") %>%
    left_join(expression_vector, by = "geneID") %>% 
    #this replaces NAs introduced by no value present in expression_vector, and replaces them with min value in expression_vector
  	mutate_at(vars(ctrl_RNAseq_expr), ~replace(., is.na(.), min(expression_vector$ctrl_RNAseq_expr))) %>% 
    mutate(norm_count = count/ctrl_RNAseq_expr) %>%
    select(-count, -ctrl_RNAseq_expr) %>%
    inner_join(DistGroup_3end, by = c("geneID", "ExonNumber")) %>% 
    group_by(Sample, DistToLandmark, DistGroup) %>% 
    summarise(sum_RNAseq_norm_count_norm_annotation_number = sum(norm_count)) %>% 
    left_join(number_of_intron_annotations_5end, by = "DistGroup") %>% 
    mutate(sum_RNAseq_norm_count_norm_annotation_number = sum_RNAseq_norm_count_norm_annotation_number/annotation_count)
}


#Load expression vector -----------------------------------

load(EXPRESSION_VECTOR_FILEPATH)

expression_vector<-left_join( 
  (as.data.frame(ctrl_RNAseq_expr) %>% 
     tibble::rownames_to_column(var = "geneID")),
  (as.data.frame(ctrl_TTseq_expr) %>% 
     tibble::rownames_to_column(var = "geneID")), by = "geneID"
) %>% 
  select(-ctrl_TTseq_expr) %>% 
  mutate(ctrl_RNAseq_expr = case_when(
    ctrl_RNAseq_expr ==0 ~ min(ctrl_RNAseq_expr[ctrl_RNAseq_expr > 0]), 
    TRUE ~  ctrl_RNAseq_expr
  ))


#load annobed -----------------------

annoBed<-read.table(ANNOTATION_BED_FILEPATH, sep = "\t", header = F)
colnames(annoBed)<-c("chr", "start", "end", "geneID", "score", "strand")

#make groupings for 3 end of exon
DistGroup_3end<-
  annoBed %>% 
  separate(geneID, into = c("geneID", "Biotype", "ExonNumber", "TotalNumberOfExons", "ExonSize", "ExonicDisanceToTSS", "GenomicDistToTSS","ExonStature", "GeneStructure"), sep=":::") %>%
  mutate(three_end = as.numeric(GenomicDistToTSS) + as.numeric(ExonSize)) %>%
    mutate(DistGroup = case_when(
      three_end < 10000 ~ "<10k",
      three_end %in% c(10000:20000) ~ "10k-20k",
      three_end %in% c(20001:30000) ~ "20k-30k",
      three_end %in% c(30001:40000) ~ "30k-40k",
      three_end %in% c(40001:50000) ~ "40k-50k",
      three_end %in% c(50001:60000) ~ "50k-60k",
      three_end %in% c(60001:70000) ~ "60k-70k",
      three_end %in% c(70001:80000) ~ "70k-80k",
      three_end %in% c(80001:90000) ~ "80k-90k",
      three_end %in% c(90001:100000) ~ "90k-100k",
      three_end %in% c(100001:110000) ~ "100k-110k",
      three_end %in% c(110001:120000) ~ "110k-120k",
      three_end %in% c(120001:130000) ~ "120k-130k",
      three_end %in% c(130001:140000) ~ "130k-140k",
      three_end %in% c(140001:150000) ~ "140k-150k",
      three_end %in% c(150001:160000) ~ "150k-160k",
      three_end %in% c(160001:170000) ~ "160k-170k",
      three_end %in% c(170001:180000) ~ "170k-180k",
      three_end %in% c(180001:190000) ~ "180k-190k",
      three_end %in% c(190001:200000) ~ "190k-200k",
      three_end %in% c(200001:210000) ~ "200k-210k",
      three_end %in% c(210001:220000) ~ "210k-220k",
      three_end %in% c(220001:230000) ~ "220k-230k",
      three_end %in% c(230001:240000) ~ "230k-240k",
      three_end %in% c(240001:250000) ~ "240k-250k",
      three_end %in% c(250001:260000) ~ "250k-260k",
      three_end %in% c(260001:270000) ~ "260k-270k",
      three_end %in% c(270001:280000) ~ "270k-280k",
      three_end %in% c(280001:290000) ~ "280k-290k",
      three_end %in% c(290001:300000) ~ "290k-300k",
      three_end > 300000 ~ ">300k"
    )) %>% 
    select(geneID, Biotype, ExonNumber,DistGroup)

#count groupings for 3 end of exon
number_of_intron_annotations_3end <-  
DistGroup_3end %>% 
  group_by(DistGroup) %>% 
  summarise(annotation_count=n())

#make groupings for 5 end of exon
DistGroup_5end<-
  annoBed %>% 
  separate(geneID, into = c("geneID", "Biotype", "ExonNumber", "TotalNumberOfExons", "ExonSize", "ExonicDisanceToTSS", "GenomicDistToTSS","ExonStature", "GeneStructure"), sep=":::") %>%
  mutate(five_end = as.numeric(GenomicDistToTSS)) %>% 
    mutate(DistGroup = case_when(
      five_end < 10000 ~ "<10k",
      five_end %in% c(10000:20000) ~ "10k-20k",
      five_end %in% c(20001:30000) ~ "20k-30k",
      five_end %in% c(30001:40000) ~ "30k-40k",
      five_end %in% c(40001:50000) ~ "40k-50k",
      five_end %in% c(50001:60000) ~ "50k-60k",
      five_end %in% c(60001:70000) ~ "60k-70k",
      five_end %in% c(70001:80000) ~ "70k-80k",
      five_end %in% c(80001:90000) ~ "80k-90k",
      five_end %in% c(90001:100000) ~ "90k-100k",
      five_end %in% c(100001:110000) ~ "100k-110k",
      five_end %in% c(110001:120000) ~ "110k-120k",
      five_end %in% c(120001:130000) ~ "120k-130k",
      five_end %in% c(130001:140000) ~ "130k-140k",
      five_end %in% c(140001:150000) ~ "140k-150k",
      five_end %in% c(150001:160000) ~ "150k-160k",
      five_end %in% c(160001:170000) ~ "160k-170k",
      five_end %in% c(170001:180000) ~ "170k-180k",
      five_end %in% c(180001:190000) ~ "180k-190k",
      five_end %in% c(190001:200000) ~ "190k-200k",
      five_end %in% c(200001:210000) ~ "200k-210k",
      five_end %in% c(210001:220000) ~ "210k-220k",
      five_end %in% c(220001:230000) ~ "220k-230k",
      five_end %in% c(230001:240000) ~ "230k-240k",
      five_end %in% c(240001:250000) ~ "240k-250k",
      five_end %in% c(250001:260000) ~ "250k-260k",
      five_end %in% c(260001:270000) ~ "260k-270k",
      five_end %in% c(270001:280000) ~ "270k-280k",
      five_end %in% c(280001:290000) ~ "280k-290k",
      five_end %in% c(290001:300000) ~ "290k-300k",
      five_end > 300000 ~ ">300k"
    )) %>% 
    select(geneID, Biotype, ExonNumber,DistGroup)

#count groupings for 3 end of exon
number_of_intron_annotations_5end <-  
DistGroup_5end %>% 
  group_by(DistGroup) %>% 
  summarise(annotation_count=n())

#load 5end count file -----------------
print(paste0("loading 5end data",IN_FILE_PATH_5_END))

int_exon_junction_coverage_five_end<-read.table(IN_FILE_PATH_5_END, sep = "\t", header = F)
colnames(int_exon_junction_coverage_five_end)<-c("Sample","chr", "start", "end", "geneID", "DistToLandmark", "strand", "count")

print("wrangling coverage of five end")
int_exon_junction_coverage_five_end<-
  int_exon_junction_coverage_five_end %>% 
  wrangle_bed_counts_for_5end() %>%
  mutate(Position = "five_end")

#load 3end count file -----------------
print(paste0("loading 5end data",IN_FILE_PATH_3_END))

int_exon_junction_coverage_three_end<-read.table(IN_FILE_PATH_3_END, sep = "\t", header = F)
colnames(int_exon_junction_coverage_three_end)<-c("Sample","chr", "start", "end", "geneID", "DistToLandmark", "strand", "count")

print("wrangling coverage of three end")
int_exon_junction_coverage_three_end <-
  int_exon_junction_coverage_three_end %>% 
  wrangle_bed_counts_for_3end() %>%
  mutate(Position = "three_end")

#rbind the files together and write to output filename

OUTFILE<-rbind(int_exon_junction_coverage_three_end,int_exon_junction_coverage_five_end)

write.table(OUTFILE, OUTFILE_PATH, quote = F, sep = "\t", row.names = F)