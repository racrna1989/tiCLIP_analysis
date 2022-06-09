#!/usr/bin/env Rscript

suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
args = commandArgs(trailingOnly=TRUE)

ANNOTATION_BED_FILEPATH<-args[1]
EXPRESSION_VECTOR_FILEPATH<-args[2]
IN_FILE_PATH_3_END<-args[3]
IN_FILE_PATH_5_END<-args[4]
OUTFILE_PATH<-args[5]



#Load expression vector ---------------------------------------------------------------------------------------

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


#load annobed --------------------------------------------------------------------------------------------------------

annoBed<-read.table(ANNOTATION_BED_FILEPATH, sep = "\t", header = F) %>% 
setNames(c("chr", "start", "end", "geneID", "score", "strand")) %>% 
    separate(geneID, into = c("geneID", "Biotype", "ExonNumber", "TotalNumberOfExons", "ExonSize", "ExonicDistance", "ExonDistFromTSS","ExonStature", "GeneStructure"), sep=":::") %>%
  mutate_at(vars(ExonDistFromTSS,ExonicDistance,ExonSize,TotalNumberOfExons,ExonNumber), .funs = as.numeric)



# identify genes for analysis --------------------------------------------------------------------------------------
top_genes<-
  annoBed %>% 
  left_join(expression_vector) %>% 
  filter(TotalNumberOfExons == 1 & !grepl("rRNA|TR_C_gene|IG_C_pseudogene|miRNA|misc_RNA", Biotype) & as.numeric(ExonSize) >49 ) %>% 
  select(geneID, ctrl_RNAseq_expr) %>% 
  unique() %>%
  mutate(ctrl_RNAseq_expr = as.numeric(ctrl_RNAseq_expr)) %>%
  arrange(desc(ctrl_RNAseq_expr)) %>% 
  filter(ctrl_RNAseq_expr > 1 & geneID != "LINC00324")

#calculate number of annotations (must have same filter as wrangle bed counts) -----------------------------------------------------------

number_of_intron_annotations<-
  annoBed %>% 
  filter(geneID %in% top_genes$geneID) %>%
  group_by(GeneStructure) %>% 
  summarise(intron_count =n())
  
#gene count ---------------------------------------------------------------
GENECOUNT<-
   annoBed %>% 
  filter(geneID %in% top_genes$geneID) %>%
  select(geneID) %>% 
  unique() %>%
  summarise(geneCount =n())

#function to process count files -----------------------------------------------------------------------------

wrangle_bed_counts <- function(dataframe) {
  dataframe %>%
  setNames(c("Sample","chr", "start", "end", "geneID", "DistToLandmark", "strand", "count")) %>%
    separate(geneID, into = c("geneID", "Biotype", "ExonNumber", "TotalNumberOfExons", "ExonSize", "ExonicDistance", "ExonDistFromTSS","ExonStature", "GeneStructure"), sep=":::") %>% 
    mutate_at(vars(ExonDistFromTSS,ExonicDistance,ExonSize,TotalNumberOfExons,ExonNumber), .funs = as.numeric) %>% 
  	filter(!grepl("CBP20_3", Sample)) %>% 
  	
  	left_join(expression_vector) %>% 
    #this replaces NAs introduced by no value present in expression_vector, and replaces them with min value in expression_vector
  	mutate_at(vars(ctrl_RNAseq_expr), ~replace(., is.na(.), min(expression_vector$ctrl_RNAseq_expr))) %>% 
    mutate(norm_count = count/ctrl_RNAseq_expr) %>% 
    
     
    filter(geneID %in% top_genes$geneID) %>%
  group_by(Sample,GeneStructure, DistToLandmark) %>% 
  summarise(sum_RNAseq_norm_count_norm_annotation_number = sum(norm_count)) %>%
  left_join(number_of_intron_annotations) %>% 
  mutate(sum_RNAseq_norm_count_norm_annotation_number = sum_RNAseq_norm_count_norm_annotation_number/intron_count) %>%
  separate(Sample, c("Protein","Rep","Timepoint","readType"), sep = "_") %>% 
  filter(Timepoint != "negative") %>%
  mutate(Timepoint_f = case_when(
    Timepoint == "PBSDRB" ~ "t00", 
    TRUE ~ Timepoint )) %>%
  mutate(region = factor(region, levels = c("5end","3end")),
         Timepoint = factor(Timepoint, levels = c("PBSDRB", "t00", "t05", "t10", "t15", "t20", "t40", "t60", "DMSO")),
         Timepoint_f = factor(Timepoint_f, levels = c("t00", "t05", "t10", "t15", "t20", "t40", "t60", "DMSO")), 
         gene_n = GENECOUNT$geneCount ) 
}


#load 5end count file -----------------

int_exon_junction_coverage_five_end<-read.table(IN_FILE_PATH_5_END, sep = "\t", header = F)

int_exon_junction_coverage_five_end<-
  int_exon_junction_coverage_five_end %>% 
  wrangle_bed_counts() %>%
  mutate(Position = "5end")

#load 3end count file -----------------

int_exon_junction_coverage_three_end<-read.table(IN_FILE_PATH_3_END, sep = "\t", header = F)

int_exon_junction_coverage_three_end <-
  int_exon_junction_coverage_three_end %>% 
  wrangle_bed_counts() %>%
  mutate(Position = "3end")

#rbind the files together and write to output filename

OUTFILE<-rbind(int_exon_junction_coverage_three_end,int_exon_junction_coverage_five_end)

write.table(OUTFILE, OUTFILE_PATH, quote = F, sep = "\t", row.names = F)





