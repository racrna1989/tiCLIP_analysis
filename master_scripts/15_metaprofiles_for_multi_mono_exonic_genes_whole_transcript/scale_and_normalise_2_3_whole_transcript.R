ScaledNormalise_2_3 <-function(DF, MINSIZE=1, MAXSIZE=10000000, annotation_filepath, rRNA_factor_filepath, total_counts_filepath, EXPRESSION_VECTOR_FILEPATH=NULL){
  
  print("starting analysis")

  #1_load annotation file and make metadata file called exons
  annoDF<-read.csv(annotation_filepath, sep = "\t", header =F) %>%
    setNames(c("chr","start","stop","ID","score","strand"))
    
    
    print("loaded annotation file")
  
  exons <-annoDF %>%
    separate(ID, c("geneName","biotype","totalExons","geneSize"), sep =":::") %>%
    select(geneName, biotype, totalExons,geneSize) %>% 
    mutate(exonDesc = ifelse(totalExons > 1, "multiExonic", "monoExonic" ),
           codingPotential = case_when(
             grepl("protein_coding|histone_coding", biotype) ~ "protein_coding",
             !grepl("protein_coding|histone_coding", biotype) & geneSize > 200 ~ "non_coding>200nt",
             !grepl("protein_coding|histone_coding", biotype) & geneSize <= 200 ~ "non_coding<200nt"
           ))
  
  
  #2_wrangle dataframe and make relative position to TSS, relative to the size of the RNA transcript
  DF_rel_pos<-DF %>%
    mutate(rel.pos   = as.numeric(rel.pos),
           geneSize = as.numeric(geneSize)) %>%
    filter(geneSize %in% c(MINSIZE:MAXSIZE)) %>%
    mutate(rel.pos.2 = rel.pos/geneSize)
  
  print("corrected relative position")
  
  #3_assigning bin numbers to the rel.pos.2. Here the min position (0) is bumped up to be included in our binning regime.
  
  ALL <-DF_rel_pos %>% mutate(rel.pos.2 = ifelse(rel.pos.2 == 0, rel.pos.2 + 0.0000001, rel.pos.2))
  max(ALL$rel.pos.2)
  min(ALL$rel.pos.2)
  ALL$new_distBin <- cut(ALL$rel.pos.2, seq(0, 1, by=0.01), labels=seq("1", "100", by = 1 ) ) 
  
  print("made bins")
  
  #4_ this step normalises the bin count to the size of that bin
  combined<-ALL %>%
    separate(Sample, into=c("Protein", "Rep", "Timepoint"), sep="_") %>%
    group_by(geneName, biotype, Protein, Timepoint, Rep, geneSize, new_distBin) %>%
    summarise(tally = n()) %>%
    mutate(tally_n = tally/(geneSize/100)) #step to normalise count to bin size. 
  

  
  print("load library sizes and normalisation factors")
  
  #5_total library sizes, and normalisation factors loaded
  RNAFactor_filePath=rRNA_factor_filepath
  totalCounts_filePath=total_counts_filepath
  
  totalLibrarySizes<-read.table(RNAFactor_filePath) %>% 
    setNames(c("Sample","rRNAFactor")) %>%
    separate(Sample, into =c("Protein","Rep","Timepoint")) %>%
    mutate(Rep = as.character(Rep)) %>%
    
    left_join(read.csv(totalCounts_filePath, sep = "\t", header = F) %>%
                setNames(c("Sample", "totalCounts")) %>%
                separate(Sample, into=c("Protein", "Rep","Timepoint")) %>%
                mutate(RPMFactor = totalCounts/1000000) %>% 
                select(-totalCounts)
    )
  
  print("calculate RPM and rRNA normalisation")
  
  #6_calculate both the RPM normalised and rRNA normalised 
  #this step taken from norm profile 3
  combined_RPM<-
    combined %>%    
    left_join(totalLibrarySizes) %>%
    ungroup() %>%
    mutate(tally_n.RPM = tally_n * RPMFactor) %>%
    select(-RPMFactor, -rRNAFactor) %>%
    group_by(geneName, biotype, Protein, Timepoint, Rep, ) %>%
    mutate(totals_unNorm    = sum(tally), 
           pct              = tally_n.RPM/sum(tally_n.RPM)*100) %>%
    ungroup() %>%
    select(Protein, Timepoint, Rep, geneName, biotype, new_distBin, totals_unNorm, pct) %>%
    spread(new_distBin, pct, fill = 0) %>%
    gather("distBin", "value", c('1':'100') )%>%
    mutate(distBin = factor(distBin, levels = c(1:100)),
           NORM = "RPM")
  
  combined_rRNA<-
    combined %>%    
    left_join(totalLibrarySizes) %>%
    mutate(tally_n.rRNA = tally_n * rRNAFactor) %>%
    select(-RPMFactor) %>%
    group_by(geneName, biotype, Protein, Timepoint, Rep,  geneSize) %>%
    mutate(totals_unNorm    = sum(tally), 
           pct              = tally_n.rRNA/sum(tally_n.rRNA)*100) %>%
    ungroup() %>%
    select(Protein, Timepoint, Rep, geneName, biotype, new_distBin, totals_unNorm, pct) %>%
    spread(new_distBin, pct, fill = 0) %>%
    gather("distBin", "value", c('1':'100') )%>%
    mutate(distBin = factor(distBin, levels = c(1:100)),
           NORM = "rRNA")
  print(names(combined_rRNA))
  print(names(combined_RPM))
  print("rbind both normalisation and calculate sum, mean coverage and number of genes used")
  #7_rbind both data frame and return
  average_coverage_normalised_for_transcript_size<-rbind(combined_RPM,combined_rRNA) %>%
    left_join(exons) %>%
    filter(!(Protein == "CBP20" & Rep == "3")) %>%
    filter(totals_unNorm > 5) %>%
    ungroup() %>%
    group_by(Protein, Timepoint, Rep, NORM, exonDesc, distBin) %>%
    summarise(count.sum = sum(value),
              count.mean = mean(value),
              count.n = n())%>% 
    ungroup() 
  
  return(average_coverage_normalised_for_transcript_size)
  
}
