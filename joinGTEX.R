library(readr)
library(dplyr)

getSubjID<-function(sampIds){
  res<-c()
  for(s in sampIds){
    #print(s)
    res<-c(res,paste(unlist(strsplit(s,"-"))[1],unlist(strsplit(s,"-"))[2],sep="-"))
  }
  return(res)
}

#join GTEX tables

GTEXSA <- read_delim("D:/MOGdata/mog_testdata/cancer/raw/combined/tcgaMetadataALL/gtex/GTEx_v7_Annotations_SampleAttributesDS.txt", 
                                                     "\t", escape_double = FALSE, trim_ws = TRUE)
GTEXSP <-read_delim("D:/MOGdata/mog_testdata/cancer/raw/combined/tcgaMetadataALL/gtex/GTEx_v7_Annotations_SubjectPhenotypesDS.txt", 
           "\t", escape_double = FALSE, trim_ws = TRUE)

#mutate GTEXSA to add SUBJID
x<-sapply(GTEXSA[,1],FUN = getSubjID)
colnames(x)<-"SUBJID"
GTEXSA<-cbind(GTEXSA[,],x)
t<-GTEXSA[,c("SUBJID","SAMPID")]

#join
joined<-left_join(GTEXSP,GTEXSA)

write_tsv(joined,"GTEX_merged.txt")


#add column analyte type R or D
library(readxl)
gtexMerged <- read_excel("gtexMerged.xlsx")

gtexMerged<-gtexMerged%>%mutate(portions.analytes.analyte_type=ifelse(grepl("RNA",SMNABTCHT),"RNA",ifelse(grepl("DNA",SMNABTCHT),"DNA","NA") ))
