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
