#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)


inFilePath<-args[1]
outFilePath<-args[2]


#load packages ------------
library(dplyr)
library(tidyr)
`%!in%` = Negate(`%in%`)

#functions ----------

wrangle_data <- function(all10kbbinrRNA) {
  all10kbbinrRNA %>%
    separate(Sample, into=c("Protein","Rep","Timepoint")) %>%
    filter(!(Protein == "CBP20" & Rep == "3")) %>%
    separate(ID, into = c("GeneID","bioType","TotalExonNumber", "size"), sep=":::") %>%
    mutate(geneType = case_when(
      TotalExonNumber > 1 ~ "multiExonic",
      TRUE ~ "monoExonic"
    )) %>%
    select(-c(chr, start,end, TotalExonNumber, size, strand)) %>%
    unite(Sample, c("Protein","Timepoint","Rep")) %>%
    unite(ID, c("GeneID","bioType","binNumber"), sep="|") %>%
    mutate(log2count = log2(count+1)) %>%
    select(Sample, ID, geneType,log2count) %>%
    separate(ID, into=c("ID", "biotype", "binNumber"), sep="\\|")
}


#preform data wrangling --------------------------------

all10kbbinrRNA<-read.csv(inFilePath, header= F, sep ="\t")

colnames(all10kbbinrRNA)<-c("Sample","chr","start","end","ID","binNumber","strand","count")

df_log2<-wrangle_data(all10kbbinrRNA)

write.table(df_log2, outFilePath, quote = F, sep = "\t", row.names = F)

