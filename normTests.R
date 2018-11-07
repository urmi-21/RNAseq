library(tximport)
library(readr)
library(edgeR)
library(tidyverse)
library(magrittr)
setwd("C:/Users/mrbai/Desktop/normrandomTest")
memory.limit(size=56000)

#real un normalized tpms an plot boxplots
changeNames<-function (thisList){
  thissuff<-paste(c("_", fnames$fname[ind]), collapse = "")
  #ind=ind+1
  assign("ind",ind+1,envir = globalenv())
  #print(ind)
  #print(thissuff)
  names(thisList)<-c("Name", paste(c("Length", thissuff), collapse = ""), paste(c("EffectiveLength", thissuff), collapse = ""),paste(c("TPM", thissuff), collapse = ""),paste(c("NumReads", thissuff), collapse = ""))
  if(ind>2){
    thisList[,3:5]
  }else{
    thisList
  }
  
}
ind<<-1
fnames<-as.data.frame(list.files(full.names = TRUE,pattern = "SRR*"))
names(fnames)<-c("fname")
fnames<-fnames %>% mutate(fname=str_replace_all(fname, "[^[:alnum:]]", ""))

df <- list.files(full.names = TRUE,pattern = "SRR*") %>% lapply(read_tsv) %>% lapply(changeNames) %>% reduce(cbind)

tpmCols<- names(df)[grepl("Name|TPM",names(df))]
df_TPM<-df[,tpmCols]

#plot histograms for tpms





read_tsv("SRR073758")%>%set_colnames(c("cyl", "disp_mean", "hp_mean"))

which(is.na(df))


gene <- read.csv("transcript_names.txt")
#Define the gene name if your transcript name is different from gene name
tx2gene <- data.frame(GENEID=gene,TXNAME=gene)
#import run ID (SRR)
run <- read.csv("runs.txt",header=F)
IDLIST <- as.character(run$V1)
file <- paste(IDLIST,"",sep="")
names(file) <- run$V1
txi.salmon <- tximport(file, type = "salmon", txOut = TRUE, tx2gene = tx2gene, importer = function(x) read_tsv(x,col_types = list(col_character(), col_double(), col_double(), col_double(), col_double()  )))
