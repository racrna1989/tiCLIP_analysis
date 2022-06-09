#!/usr/bin/env Rscript

#get current wd
WD<-getwd()

#Make functions -----------------------
#wrangle data 
readAndWrangleDataFrame <- function(FILEPATH){
  print("reading in data")
  DF<-read.table(FILEPATH)
  colnames(DF)<-c("sample","chr","start","end","geneName","binNumber","strand","count")
  print(paste("read in ",FILEPATH))
  
  DF_1<-DF %>%
    select(-chr,-start,-end,-strand) %>%
    separate(sample, into=c("Protein","Rep","Timepoint")) %>%
    separate(geneName, into=c("geneName", "bioType", "exonID","genomicSize"), sep="\\:::") %>% 
    mutate(genomicSize = as.numeric(genomicSize)) %>%
    mutate(sizeRange = case_when( 
      genomicSize < 10000 ~ "<10k",
      genomicSize %in% c(10000:20000) ~ "10k-20k",
      genomicSize %in% c(20001:30000) ~ "20k-30k",
      genomicSize %in% c(30001:40000) ~ "30k-40k",
      genomicSize %in% c(40001:50000) ~ "40k-50k",
      genomicSize %in% c(50001:60000) ~ "50k-60k",
      genomicSize %in% c(60001:70000) ~ "60k-70k",
      genomicSize %in% c(70001:80000) ~ "70k-80k",
      genomicSize %in% c(80001:90000) ~ "80k-90k",
      genomicSize %in% c(90001:100000) ~ "90k-100k",
      genomicSize %in% c(100001:110000) ~ "100k-110k",
      genomicSize %in% c(110001:120000) ~ "110k-120k",
      genomicSize %in% c(120001:130000) ~ "120k-130k",
      genomicSize %in% c(130001:140000) ~ "130k-140k",
      genomicSize %in% c(140001:150000) ~ "140k-150k",
      genomicSize %in% c(150001:160000) ~ "150k-160k",
      genomicSize %in% c(160001:170000) ~ "160k-170k",
      genomicSize %in% c(170001:180000) ~ "170k-180k",
      genomicSize %in% c(180001:190000) ~ "180k-190k",
      genomicSize %in% c(190001:200000) ~ "190k-200k",
      genomicSize > 200000 ~ ">200k"
    ))
  
  return(DF_1)
  
}

#calculate coverage using 60kb to 300kb genes, bins 1-200, selecting genes that have over 1 count across all samples as judged by intronic reads. It also removes bins 1:3, and the last 3 bins. 
# mean trim 0.0025

calculateCoverage <- function(DF) {
  DF_1<- DF  %>%
    #remove CBP20-3 from dataset #only look at genes over 60k & only keep bins between 1 to 200 & add
    filter(!(Protein == "CBP20" & Rep == "3")) %>%
    #genomicSize %in% c(60000:300000) & 
    #binNumber %in% c(1:200) &
    #geneName %in% geneListOverZeroCounts$geneName) %>%
    #add total bins to table 
    left_join(mRNATotalBins) %>%
    #select the relevant columns (only one read type here) 
    select(Protein, Rep, Timepoint, geneName,TotalBins, binNumber,sizeRange, count) %>%
    #filter(binNumber %ni% c(1:3) &
    #         binNumber != TotalBins & 
    #         binNumber != TotalBins -1 &
    #         binNumber != TotalBins -2 )%>%
    #group together samples, sizeRange and genomic contexts along with binNumber # tested and mean 0.0025 is best
    group_by(Protein, Rep, Timepoint,sizeRange) %>%
    mutate(n=n_distinct(geneName)) %>%
    ungroup() %>%
    group_by(Protein, Rep, Timepoint,binNumber, sizeRange,n) %>%
    summarise(binSum=sum(count),
              mean  =mean(count),
              mean_trim0.0025 = mean(count, trim = 0.0025)) %>%
    ungroup() %>%
    #factorise and add levels
    mutate( binNumber = as.numeric(as.character(binNumber)),
            Timepoint = factor(Timepoint, levels=c("negative", "PBSDRB", "t00","t05", "t10", "t15", "t20", "t40", "t60", "DMSO")))
  
  return(DF_1)
}



#load libraries-------------------------------
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))

#define arguments
args = commandArgs(trailingOnly=TRUE)

DATA_FRAME_FILEPATH<-args[1]
OUTFILE_PATH<-args[2]

#Process data --------------------------

#STEP_1 load and wrangle the data
#read and wrangle data from filepath input
wrangled<-readAndWrangleDataFrame(DATA_FRAME_FILEPATH)

#STEP_2 calculate total bins for each gene
#calculate total bins for each gene
mRNATotalBins<-wrangled %>%
  select(geneName,binNumber, count) %>%
  group_by(geneName) %>%
  summarise(TotalBins = max(as.numeric(binNumber)))

#STEP_3 make a list of genes with over zero counts. this uses the total of all datasets. 
#make list of genes with over zero counts

geneListOverZeroCounts<-wrangled %>%
  select(Protein, Timepoint, geneName, count) %>%
  group_by(geneName) %>%
  summarise(totalCounts = sum(as.numeric(count))) %>% 
  filter(totalCounts > 0) %>% 
  select(geneName) %>% 
  unique()

#STEP_4 
#calculate coverage
wrangled_calculated<-calculateCoverage(wrangled)

write.table(wrangled_calculated,OUTFILE_PATH, quote = F, sep = "\t", row.names = F)

#-------------
