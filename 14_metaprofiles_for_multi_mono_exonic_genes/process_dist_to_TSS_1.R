process_distToTSS_1 <- function(data_frame_filepath, annotation_filepath){
  
  #1_add colnames to vectors
  col.names<-c("Sample","chrAnno","startAnno", "stopAnno", "IDAnno", "scoreAnno", "strandAnno","chrRead","startRead","stopRead","IDRead","scoreRead","strandRead")
  col.numeric<- c("chrAnno","startAnno","stopAnno","totalExons","exonSize","cumSumExons","distToTSS","startRead","stopRead")  
  
  
  #2_load counts file using dfFilePath input. 
  df<-read.csv(data_frame_filepath, sep = "\t", header =F)
  
  print("wrangling data") 
  
  #3_name columns, separate relevant feilds, and make defined columns numeric
  df1<-df %>%
    setNames(col.names) %>%
    select(-chrRead, -IDRead, -scoreRead, -scoreAnno) %>%
    separate(IDAnno, c("geneName","biotype","exonID", "totalExons", "exonSize","cumSumExons", "distToTSS", "exonDesc","geneDesc"), sep =":::") %>%
    select(-exonDesc) 
  
  df1[col.numeric] <- sapply(df1[col.numeric],as.numeric)
  
  #4_calculate distance from xlink to TSS using the cumulative sum of upstream exons. 
  print("processing rel distance from TSS")
  df2<-df1 %>% 
    mutate(rel.pos = case_when(
      strandAnno == "+" ~ ( (startRead-startAnno) + (cumSumExons-exonSize) ), 
      strandAnno == "-"  ~ ( (stopAnno-stopRead) + (cumSumExons-exonSize) )
    ))%>%
    select(Sample, geneName, biotype, geneDesc, rel.pos)
  
  
  
  #5_add genomic and mature RNA sizes to protein coding genes. 
  print("loading annotation file")
  print("making table with total sizes of RNAs")
  
  annoDF<-read.csv(annotation_filepath, sep = "\t", header =F)
  
  totalSizes<-annoDF %>%
    setNames(c("chr","start","stop","ID","score","strand")) %>%
    separate(ID, c("geneName","biotype","exonID", "totalExons", "exonSize","cumSumExons", "distToTSS", "exonDesc","geneDesc"), sep =":::") %>%
    filter(exonID == totalExons) %>%
    mutate(matureRNA = as.numeric(cumSumExons),
           geneSize = as.numeric(ifelse(totalExons > 1, distToTSS, exonSize))) %>%
    select(geneName, biotype, matureRNA, geneSize)
  
  df3<-df2 %>%
    left_join(totalSizes)
  
  return(df3)
  
}
