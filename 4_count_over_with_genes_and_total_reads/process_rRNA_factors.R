library(dplyr)
library(tidyr)

  args <- commandArgs(trailingOnly = TRUE)
  inCountFile = args[1]
  outFile     = args[2]
    
  wrangled_table <- read.table(inCountFile) %>% 
    setNames(c("Sample", "chr","start","stop","name","dot", "strand", "count")) %>%
    select(-chr,-start,-stop,-dot,-strand) %>% 
    separate(Sample, into=c("Protein", "Rep", "Timepoint")) %>% 
    filter(grepl("rRNA:::", name)) %>%
    separate(name, into=c("geneID","biotype", "geneName"), sep ="\\:::")%>% 
    filter(Timepoint != "val")
  
  rRNA.list<-c("ENSG00000199839",
               "ENSG00000202264",
               "ENSG00000199523",
               "ENSG00000199480",
               "ENSG00000199415",
               "ENSG00000201059",
               "ENSG00000278189",
               "ENSG00000200558",
               "ENSG00000201321",
               "ENSG00000199994",
               "ENSG00000200408",
               "ENSG00000272435",
               "ENSG00000274917",
               "ENSG00000272351",
               "ENSG00000201185",
               "ENSG00000210082",
               "ENSG00000211459",
               "ENSG00000275215",
               "ENSG00000275757",
               "ENSG00000276700")
  
 rRNACounts <- wrangled_table %>%
  filter(geneID %in% rRNA.list) %>%
  group_by(Protein, Rep, Timepoint) %>% 
  summarise(rRNA.counts = sum(count)) %>% 
  arrange(desc(rRNA.counts)) %>% 
  mutate(rRNAFactor = 30000/rRNA.counts) %>% 
  unite("Sample", c(Protein:Timepoint), sep = "_") %>% 
  select(Sample, rRNAFactor)


write.table(rRNACounts, outFile, quote = F, col.names = F, row.names =F)
