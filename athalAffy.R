library("ArrayExpress", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.5")
library(affy)
library(dplyr)
library(gcrma)
library(arrayQualityMetrics)

#query
ATsets = queryAE(species = "Arabidopsis+thaliana")
#keep ones with raw data
ATsetsfiltered<-ATsets %>% filter(tolower(Species)=="arabidopsis thaliana") %>% filter(tolower(Raw)=="yes")


smallATset<-head(ATsetsfiltered)

#download all .cel files
failed<-c()
for(id in smallATset$ID){
  print(id)
  atSamp = getAE(id, type = "raw")
  thisRawList<-atSamp$rawFiles
  #if there are no .cel files
  celInd<-grep("\\.CEL$",thisRawList)
  if(length(celInd)<1){
    print("Err")
    failed<-c(failed,id)
  }else{
    #read data
    affy.data = ReadAffy(filenames = thisRawList)
    #apply gcrma
    eset <- gcrma(affy.data)
    #get metadata
    thisPD<-pData(eset)
    #get expression
    thisExp<-exprs(eset)
  }
}

files = list.files(pattern = "\\.CEL$", full.names = TRUE)
affy.data = ReadAffy(filenames = files)

#apply gcrma
eset <- gcrma(affy.data)
