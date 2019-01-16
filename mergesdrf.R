library(readr)
library(magrittr)
library(dplyr)
library(purrr)
library(data.table)

#megre sdrf files downloaded from arrayexpress
df_sdrf<-list.files(full.names = TRUE,path = "sdrfDir/", pattern = "*.txt") %>% lapply(read_tsv) %>% makeasChar %>% bind_rows()


makeasChar<-function(listOfDF){
  t<-length(listOfDF)
  
  for(i in 1:t){
    listOfDF[[i]]<- listOfDF[[i]] %>% mutate_all(as.character)
  }
  return(listOfDF)
}

