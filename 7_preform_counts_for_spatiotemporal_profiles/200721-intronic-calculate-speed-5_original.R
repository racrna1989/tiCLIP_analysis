#!/usr/bin/env Rscript


suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(ggplot2))
suppressMessages(library(scales))
suppressMessages(library(ggpubr))   


`%ni%` = Negate(`%in%`)

dir.create("figs")
dir.create"tables")
list.files()


##-------create functions -------

#wrangle data 
readAndWrangleDataFrame <- function(FILEPATH){
  
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
  DF_1 <- DF  %>%
    #remove CBP20-3 from dataset #only look at genes over 60k & only keep bins between 1 to 200 & add
    filter(!(Protein == "CBP20" & Rep == "3") &
             genomicSize %in% c(60000:300000) & 
             binNumber %in% c(1:200) &
             geneName %in% geneListOverZeroCounts$geneName) %>%
    #add total bins to table 
    left_join(mRNATotalBins) %>%
    #select the relevant columns (only one read type here) 
    select(Protein, Rep, Timepoint, geneName,TotalBins, binNumber, count) %>%
    filter(binNumber %ni% c(1:3) &
             binNumber != TotalBins & 
             binNumber != TotalBins -1 &
             binNumber != TotalBins -2 )%>%
    #group together samples, sizeRange and genomic contexts along with binNumber # tested and mean 0.0025 is best
    group_by(Protein, Rep, Timepoint,binNumber) %>%
    summarise(binSum=sum(count),
              mean  =mean(count),
              mean_trim0.0025 = mean(count, trim = 0.0025)) %>%
    ungroup() %>%
    #factorise and add levels
    mutate( binNumber = as.numeric(as.character(binNumber)),
            Timepoint = factor(Timepoint, levels=c("negative", "PBSDRB", "t00","t05", "t10", "t15", "t20", "t40", "t60", "DMSO")))
  
  return(DF_1)
}

# get the value for a given bin, this will be used to normalise the signal. This bin number should be number 4. As judged by testing of normalisation 
getBinMean <- function(DF, BINNUMBER) {
  DF <- DF %>% 
    filter(binNumber == BINNUMBER) %>% 
    select(Protein, Rep, Timepoint, mean_trim0.0025) %>%
    rename("meanAtBin" = mean_trim0.0025)
}

# make graph out put of sample. using 15 degrees of freedom polynormalisation. The normalisation it preformed in this function
makeGraph <- function(DF) {
  
  fig<-DF %>% 
    #filter(Timepoint %ni% c("negative","PBSDRB","t00","t05","DMSO")) %>%
    left_join(BinMean) %>%
    mutate(normToBin= mean_trim0.0025/meanAtBin) %>%
    ggplot() +
    geom_line(position="identity", stat="identity", alpha = 0.3, size=0.5, aes(x=binNumber, y=normToBin, linetype=Rep)) +
    geom_smooth(method="lm", aes(x=binNumber, y=normToBin, col=Protein, group=Protein), formula=y~poly(x, 15), alpha = 0.75, size=0.6) +
    #geom_line(position="identity", stat="identity", alpha = 0.2, size=0.2, aes(x=binNumber, y=normToMax, col=Timepoint)) +
    facet_grid(Timepoint~Protein, scales = "free_y" ) + 
    labs(col = "polynomal df=15") + 
    theme_bw()
  
  print(fig)
}

# make graph out put of sample. using 15 degrees of freedom polynormalisation. The normalisation it preformed in this function.

get_spatioTemporal_graph<- function(YAXISCOLUMN) {
  ggplot() +
    geom_line(position="identity", stat="identity", alpha = 0.3, size=0.5, aes_string(x="binNumber", y=YAXISCOLUMN, linetype="Rep")) +
    geom_smooth(method="lm", aes_string(x="binNumber", y=YAXISCOLUMN, col="Protein", group="Protein"), formula=y~poly(x, 15), alpha = 0.75, size=0.6) +
    #geom_line(position="identity", stat="identity", alpha = 0.2, size=0.2, aes(x=binNumber, y=normToMax, col=Timepoint)) +
    facet_grid(Timepoint~Protein) + 
    labs(col = "polynomal df=15") + 
    coord_cartesian(ylim=c(0,2)) +
    theme_bw()
}

#---process intronic data -----


#STEP_1 load and wrangle the data
#read and wrangle data from filepath input
mRNA_1kb_intronic_counts_wrangled<-readAndWrangleDataFrame("all.hg38_HeLa_Soren.1kbins.200409.sense.counts")

#STEP_2 calculate total bins for each gene
#calculate total bins for each gene
mRNATotalBins<-mRNA_1kb_intronic_counts_wrangled %>%
  select(geneName,binNumber, count) %>%
  group_by(geneName) %>%
  summarise(TotalBins = max(as.numeric(binNumber)))

#STEP_3 make a list of genes with over zero counts. this uses the total of all datasets. 
#make list of genes with over zero counts
geneListOverZeroCounts<-mRNA_1kb_intronic_counts_wrangled %>%
  select(Protein, Timepoint, geneName, count) %>%
  group_by(geneName) %>%
  summarise(totalCounts = sum(as.numeric(count))) %>% 
  filter(totalCounts > 0) %>% 
  select(geneName) %>% 
  unique()

#STEP_4 
#calculate coverage
mRNA_1kb_intronic_counts_wrangled_calculated<-calculateCoverage(mRNA_1kb_intronic_counts_wrangled)

#optional step 
#write down intermediate file
write.table(mRNA_1kb_intronic_counts_wrangled_calculated, "all.hg38_HeLa_Soren.1kbins.200409_calculatedCoverage.tab", quote = F)

#STEP5 calculate the bin mean at a given bin
BinMean<-getBinMean(mRNA_1kb_intronic_counts_wrangled_calculated, 4)

#STEP6 calculate the bin mean at a given bin
makeGraph(mRNA_1kb_intronic_counts_wrangled_calculated) + 
  labs(subtitle = "spatiotemporal RNA binding maps - intronic reads")

#------------------
#next steps are to preform the processing steps for the spatiotemporal maps, but to count the number of reads used for each map, and compare it with the number of reads excluded from the analysis (these are reads that over exons, or reads that fall within the first and last 3 bins)

#STEP1 COUNTS
#replicate the processing for above but for COUNTS

intronic_count_intermediate <-mRNA_1kb_intronic_counts_wrangled %>%
  #remove CBP20-3 from dataset #only look at genes over 60k & only keep bins between 1 to 200 & add
  filter(!(Protein == "CBP20" & Rep == "3") &
           genomicSize %in% c(60000:300000) & 
           binNumber %in% c(1:200) &
           geneName %in% geneListOverZeroCounts$geneName & 
           count > 0 ) %>% 
  left_join(mRNATotalBins)

#optional STEP for counts 

#label the ends that are removed, then sum the counts associated with the ends and the intermediat potion used 
write.table(intronic_count_intermediate, "intronic_count_intermediate.tab", quote = F, sep = "\t", row.names = F)


#STEP2 COUNTS
#segment the count file into internal counts and counts on the first and last 3 bins 

intronic_count<-intronic_count_intermediate %>%
  mutate(binPosition = case_when(
    binNumber %in% c(1:3) ~ "intronic-ends",
    binNumber < (TotalBins - 2) ~ "intronic-internal",
    binNumber >= (TotalBins - 2) ~ "intronic-ends"
  )) %>% 
  group_by(Protein, Rep, Timepoint, binPosition) %>% 
  summarise(count = sum(count))

#STEP3 COUNTS
#print the output in a graph that shows the distribution of intronic end reads over internal reads. 
intronic_count %>%  
  ggplot() + 
  geom_bar(aes(x=Timepoint, y=count, fill = binPosition), position = "fill", stat = "identity") + 
  facet_grid(Protein ~ Rep) + 
  labs(subtitle = "intronic reads summed (rRNAScaled)")


############EXONIC COUNTS ################
# make a vector of all the file paths 

# STEP 1 for exonic counts is to read in each proteins exonic counts into a loop that loads, and wrangles the file using the same readAndWrangeDataFrame function
#This loop also calculates max bins for each gene. NOTE no function for calculating gene list. THe gene list is used from intronic reads. 

files<-list.files("../200713-200810-exonic-mRNA-2-1kb", include.dirs = TRUE, full.names = TRUE)
exonic_list<-c()

for(FILEPATH in files){
  
  mRNA_1kb_exonic_counts_wrangled<-readAndWrangleDataFrame(FILEPATH)
  
  #calculate total bins for each gene
  mRNATotalBins<-mRNA_1kb_exonic_counts_wrangled %>%
    select(geneName,binNumber, count) %>%
    group_by(geneName) %>%
    summarise(TotalBins = max(as.numeric(binNumber)))
  
  exonic_list[[FILEPATH]]<-mRNA_1kb_exonic_counts_wrangled
  
}  

mRNA_1kb_exonic_counts_wrangled<-rbindlist(exonic_list)


#STEP 2 replicate the processing for above COUNTS

exonic_count_intermediate <- mRNA_1kb_exonic_counts_wrangled %>%
  #remove CBP20-3 from dataset #only look at genes over 60k & only keep bins between 1 to 200 & add
  filter(!(Protein == "CBP20" & Rep == "3") &
           genomicSize %in% c(60000:300000) & 
           binNumber %in% c(1:200) &
           geneName %in% geneListOverZeroCounts$geneName & 
           count > 0 ) %>% 
  left_join(mRNATotalBins)

#Optional step to save file 
write.table(exonic_count_intermediate, "exonic_count_intermediate.tab", quote = F, sep = "\t", row.names = F)

#STEP 3 label 
#label the ends that are removed, then sum the counts associated with the ends and the intermediat potion used 

exonic_count<-exonic_count_intermediate %>%
  mutate(binPosition = "exonic") %>%
  group_by(Protein, Rep, Timepoint, binPosition) %>% 
  summarise(count = sum(count))

exonic_count %>%
  ggplot() + 
  geom_bar(aes(x=Timepoint, y=count, fill = binPosition), position = "fill", stat = "identity") + 
  facet_grid(Protein ~ Rep) + 
  labs(subtitle = "exonic reads summed (rRNAScaled)")


#STEP 4 combine the intronic and exonic counts and output graphs displaying raw numbers (rRNA adjusted), the stacked bar chart, percentage of reads, and tables. 
pdf("spatioTemporalCoverage_counts_induv.pdf")
rbind(exonic_count, intronic_count) %>% 
  mutate(Timepoint = factor(Timepoint, levels=c("negative", "PBSDRB", "t00","t05", "t10", "t15", "t20", "t40", "t60", "DMSO"))) %>%
  ggplot() + 
  geom_bar(aes(x=Timepoint, y=count/1000000, fill = binPosition), position = "stack", stat = "identity") + 
  facet_grid(Protein ~ Rep, scales = "free_y") + 
  labs(fill = "read position",
       subtitle = "rRNA adjusted counts of reads covering genes used for spatiotemporal analysis") +
  ylab("counts (10^6)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

rbind(exonic_count, intronic_count) %>% 
  mutate(Timepoint = factor(Timepoint, levels=c("negative", "PBSDRB", "t00","t05", "t10", "t15", "t20", "t40", "t60", "DMSO"))) %>%
  ggplot() + 
  geom_bar(aes(x=Timepoint, y=count/1000000, fill = binPosition), position = "dodge", stat = "identity") + 
  facet_grid(Protein ~ Rep, scales = "free_y") + 
  labs(fill = "read position",
       subtitle = "rRNA adjusted counts of reads covering genes used for spatiotemporal analysis") + 
  ylab("counts (10^6)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

rbind(exonic_count, intronic_count) %>% 
  mutate(Timepoint = factor(Timepoint, levels=c("negative", "PBSDRB", "t00","t05", "t10", "t15", "t20", "t40", "t60", "DMSO"))) %>%
  ggplot() + 
  geom_bar(aes(x=Timepoint, y=count, fill = binPosition), position = "fill", stat = "identity") + 
  facet_grid(Protein ~ Rep, scales = "free_y") + 
  labs(fill = "read position",
       subtitle = "rRNA adjusted counts of reads covering genes used for spatiotemporal analysis") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off() 

#-----------------------------------------------------------
pdf("spatioTemporalCoverage_counts_all.pdf")
rbind(exonic_count, intronic_count) %>% 
  mutate(Timepoint = factor(Timepoint, levels=c("negative", "PBSDRB", "t00","t05", "t10", "t15", "t20", "t40", "t60", "DMSO"))) %>%
  ggplot() + 
  geom_bar(aes(x=Timepoint, y=count/1000000, fill = binPosition), position = "stack", stat = "summary") + 
  facet_grid(Protein ~ ., scales = "free_y") + 
  labs(fill = "read position",
       subtitle = "rRNA adjusted counts of reads covering genes used for spatiotemporal analysis") +
  ylab("counts (10^6)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

rbind(exonic_count, intronic_count) %>% 
  mutate(Timepoint = factor(Timepoint, levels=c("negative", "PBSDRB", "t00","t05", "t10", "t15", "t20", "t40", "t60", "DMSO"))) %>%
  ggplot() + 
  geom_bar(aes(x=Timepoint, y=count/1000000, fill = binPosition), position = "dodge", stat = "summary") + 
  facet_grid(Protein ~ ., scales = "free_y") + 
  labs(fill = "read position",
       subtitle = "rRNA adjusted counts of reads covering genes used for spatiotemporal analysis") + 
  ylab("counts (10^6)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

rbind(exonic_count, intronic_count) %>% 
  mutate(Timepoint = factor(Timepoint, levels=c("negative", "PBSDRB", "t00","t05", "t10", "t15", "t20", "t40", "t60", "DMSO"))) %>%
  ggplot() + 
  geom_bar(aes(x=Timepoint, y=count, fill = binPosition), position = "fill", stat = "summary") + 
  facet_grid(Protein ~ ., scales = "free_y") + 
  labs(fill = "read position",
       subtitle = "rRNA adjusted counts of reads covering genes used for spatiotemporal analysis") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

dev.off() 


## wrangle data to output tables 
#
#dataForTable<-rbind(exonic_count, intronic_count) %>% 
#  group_by(Protein, Rep, Timepoint) %>%
#  mutate(totalReads = sum(count)) %>% ungroup() %>%
#  unite("Protein-Rep", c(Protein, Rep), sep ="-") %>%
#  mutate(count = (count / totalReads)) %>% 
#  spread(binPosition, count)
#
#
#sample_list<-unique(dataForTable$`Protein-Rep`)
#
#for(SAMPLE in sample_list){
#  tab<-dataForTable %>% 
#    mutate(Timepoint = factor(Timepoint, levels=c("negative", "PBSDRB", "t00","t05", "t10", "t15", "t20", "t40", "t60", "DMSO"))) %>%
#    arrange(Timepoint) %>%
#    filter(`Protein-Rep` == SAMPLE) %>%
#    select(-`Protein-Rep`) %>%
#    gt() %>%
#    fmt_percent(columns = vars("exonic",'intronic-ends', 'intronic-internal'), decimal = 0) %>%
#    fmt_number(columns = vars('totalReads'), decimal = 0) %>%
#    tab_header(
#      title = paste0(SAMPLE)
#    )
#  gtsave(tab, paste0("tables/",SAMPLE,"-reads-used-for-spatiotemporal-coverage.png"))
#  print(tab)
#}
#




