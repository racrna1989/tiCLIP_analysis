process_distToTSS_1 <- function(data_frame_filepath, annotation_filepath=NULL){
  
  #1_add colnames to vectors
  col.names<-c("Sample","chrAnno","startAnno", "stopAnno", "IDAnno", "scoreAnno", "strandAnno","chrRead","startRead","stopRead","IDRead","scoreRead","strandRead")
  col.numeric<- c("chrAnno","startAnno","stopAnno","totalExons","geneSize","startRead","stopRead")
  
  
  #2_load counts file using dfFilePath input. 
  df<-read.csv(data_frame_filepath, sep = "\t", header =F)
  
  print("wrangling data") 
  
  #3_name columns, separate relevant feilds, and make defined columns numeric
  df1<-df %>%
    setNames(col.names) %>%
    select(-chrRead, -IDRead, -scoreRead, -scoreAnno) %>%
    separate(IDAnno, c("geneName","biotype","totalExons","geneSize"), sep =":::")
  
  df1[col.numeric] <- sapply(df1[col.numeric],as.numeric)
  
  #4_calculate distance from xlink to TSS using the cumulative sum of upstream exons. 
  
  print("processing rel distance from TSS")
  
  df2<-df1 %>%
    mutate(rel.pos = case_when(
      strandAnno == "+" ~ (startRead-startAnno),
      strandAnno == "-" ~ (stopAnno-stopRead)
    ))%>%
    select(Sample, geneName, biotype, geneSize, totalExons, rel.pos)
    
    print("done")

  df3<-df2
  
  return(df3)
  
}
