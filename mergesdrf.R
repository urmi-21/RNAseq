library(readr)
library(magrittr)
library(dplyr)
library(purrr)

makeasChar<-function(listOfDF){
  t<-length(listOfDF)
  
  for(i in 1:t){
    listOfDF[[i]]<- listOfDF[[i]] %>% mutate_all(as.character)
  }
  return(listOfDF)
}


#megre sdrf files downloaded from arrayexpress
df_sdrf<-list.files(full.names = TRUE,path = "sdrfDir/", pattern = "*.txt") %>% lapply(read_tsv) %>% makeasChar %>% bind_rows()

#rename col
colnames(df_sdrf)[which(colnames(df_sdrf)=="Array Data File")] <- "Array_Data_File"
df_sdrf[is.na(df_sdrf)]<-"NA"

#check all rows have Array Data File
sum(is.na(df_sdrf$Array_Data_File))
df_sdrf$`Source Name`[which(df_sdrf$Array_Data_File=="NA")]

#read cel files in the data
sample_to_accession_human <- read_delim("~/Downloads/Immuno-Navigator_datasets/human/sample_to_accession_human.txt", 
                                        "\t", escape_double = FALSE, col_names = FALSE, 
                                        trim_ws = TRUE)
#check if all cel files are in sdrf data
not_in_sdrf<- sample_to_accession_human %>% filter(!(X1 %in% df_sdrf$Array_Data_File))
#result is zero

#reduce sdrf to the data
df_sdrf_inData<-df_sdrf %>% filter(Array_Data_File %in% sample_to_accession_human$X1)

#which cel name are duplcated
df_sdrf_inData$Array_Data_File[duplicated(df_sdrf_inData$Array_Data_File)]
#E-GEOD-13501.sdrf.txt seems to have multiple rows for same cel files


#consolidate rows based on cel files
df_sdrf_inData_noRep <- aggregate( .~ Array_Data_File, df_sdrf_inData, function(x) paste0(unique(x),collapse = ";"))


#extract useful info
keyWords<-c("organism","cell","disease","treatment","arrayexpress","Array_Data_File")

colsToKeep<-grepl(paste(keyWords,collapse = "|"),colnames(df_sdrf_inData_noRep),ignore.case = T)


df_sdrf_final<-df_sdrf_inData_noRep[,colsToKeep]

#order by col name
df_sdrf_final<-df_sdrf_final[,colnames(df_sdrf_final)[order(colnames(df_sdrf_final))]]

#merge column information
df_sdrf_final[df_sdrf_final=="NA"]<-""

#init new df
test<-data.frame()
test<-as.data.frame(df_sdrf_final[,"Characteristics [cell group]"])


for(keyw in keyWords){
print(keyw)
thisColName<-paste(keyw,"Info",sep = "_")
thisCols<-colnames(df_sdrf_final)[grepl(keyw,colnames(df_sdrf_final),ignore.case = T)]
tempdf<-as.data.frame(df_sdrf_final[ , thisCols ])
test[,thisColName] <- apply( tempdf , 1 , function(row) paste(row[nzchar(row)], collapse = ";") )

}

test<-test[,2:ncol(test)]

#join accession
colnames(sample_to_accession_human)<-c("Array_Data_File","Accession")
sample_to_accession_human<-sample_to_accession_human[,c("Accession","Array_Data_File")]
colnames(test)[which(colnames(test)=="Array_Data_File_Info")]<-"Array_Data_File"
test<-full_join(sample_to_accession_human,test)
#open batch information
sample_to_batches_human <- read_delim("~/Downloads/Immuno-Navigator_datasets/human/sample_to_batches_human.txt", 
                                      "\t", escape_double = FALSE, col_names = FALSE, 
                                      trim_ws = TRUE)
colnames(sample_to_batches_human)<-c("Array_Data_File","Batch")
test<-full_join(test,sample_to_batches_human)
#open celltype info
sample_to_cell_types_human <- read_delim("~/Downloads/Immuno-Navigator_datasets/human/sample_to_cell_types_human.txt", 
                                         "\t", escape_double = FALSE, col_names = FALSE, 
                                         trim_ws = TRUE)
colnames(sample_to_cell_types_human)<-c("Array_Data_File","CellType")
test<-full_join(test,sample_to_cell_types_human)

write_tsv(test,"combined_srdf_Summary.tsv")














