library(tximport)
library(readr)
library(edgeR)
library(tidyverse)
library(magrittr)
library(ggplot2)
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
#ggplot(stack(df_TPM[1:194]), aes(x = ind, y = values)) +  geom_boxplot()






# function to produce summary statistics (mean and +/- sd), as required for ggplot2
data_summary <- function(x) {
  mu <- mean(x)
  med<-median(x)
  #print(paste("mean is",log(mu)))
  return(c(y=log(mu),ymax=log(med),ymin=log(med)))
}

getlogmax<- function(x){
  return(c(y=log(max(x))))
}

getlogmin<- function(x){
  return(c(y=log(min(x))))
}

newdf<-df_TPM[,2:3]+1

ggplot(data=stack(newdf), aes(x=ind, y=(values))) + geom_crossbar(stat="summary", fun.y=data_summary, fun.ymax=getlogmax, fun.ymin=getlogmin) 

ggplot(data=stack(newdf), aes(x=ind, y=log(values))) +   geom_violin() 


plotlist = list()
k=1
for(i in seq(2, dim(df_TPM)[2], by = 200)){
  df_s<-df_TPM[,c(i:min(i+200,dim(df_TPM)[2]))]+1
   #p<-ggplot(stack(df_s), aes(x = ind, y = log(values))) +  geom_boxplot(aes(middle = mean(values)),coef = 6) +theme(axis.text.x = element_text(angle = 90, hjust = 1))
   p<-ggplot(data=stack(df_s), aes(x=ind, y=(values))) + geom_crossbar(stat="summary", fun.y=data_summary, fun.ymax=getlogmax, fun.ymin=getlogmin) +theme(axis.text.x = element_text(angle = 90, hjust = 1))
   plotlist[[k]]=p
  k=k+1
}

pdf("RawTPMPlots.pdf")
for (i in 1:length(plotlist)){
  write(i, stderr())
  print(plotlist[[i]])
}
dev.off()

#sort by depth and plot
bases <- read_delim("bases.txt", "\t", escape_double = FALSE,trim_ws = TRUE)
bases<-bases %>% mutate(srr=paste("TPM_",run_accession,sep = ""))
bases <- bases[bases$srr %in% names(df_TPM),]
bases <- bases[order(bases$bases),]

df_TPM_sorted<-df_TPM[,bases$srr]

plotlist = list()
k=1
for(i in seq(2, dim(df_TPM_sorted)[2], by = 50)){
  df_s<-df_TPM_sorted[,c(i:min(i+50,dim(df_TPM_sorted)[2]))]+1 #add 1 to avoid inf
  p<-ggplot(data=stack(df_s), aes(x=ind, y=(values))) + geom_crossbar(stat="summary", fun.y=data_summary, fun.ymax=getlogmax, fun.ymin=getlogmin) +theme(axis.text.x = element_text(angle = 90, hjust = 1))
  plotlist[[k]]=p
  k=k+1
}

pdf("RawTPMPlotsSortedBase_50.pdf")
for (i in 1:length(plotlist)){
  write(i, stderr())
  print(plotlist[[i]])
}
dev.off()

###separate tpms of Ens genes and others

enstNames<-df_TPM$Name[(grepl("ENST",df_TPM$Name))]

df_TPM_enst<-df_TPM[df_TPM$Name %in% enstNames,]

df_TPM_nonenst<-df_TPM[!(df_TPM$Name %in% enstNames),]
#sort by bases
df_TPM_enst<-df_TPM_enst[,c('Name',bases$srr)]
df_TPM_nonenst<-df_TPM_nonenst[,c('Name',bases$srr)]

plotlist = list()
k=1
for(i in seq(2, dim(df_TPM_enst)[2], by = 200)){
  df_s<-df_TPM_enst[,c(i:min(i+200,dim(df_TPM_enst)[2]))]+1
  p<-ggplot(data=stack(df_s), aes(x=ind, y=(values))) + geom_crossbar(stat="summary", fun.y=data_summary, fun.ymax=getlogmax, fun.ymin=getlogmin) +theme(axis.text.x = element_text(angle = 90, hjust = 1))
  plotlist[[k]]=p
  k=k+1
}
pdf("RawTPMPlots_Enst.pdf")
for (i in 1:length(plotlist)){
  write(i, stderr())
  print(plotlist[[i]])
}
dev.off()
#plot non enst
plotlist = list()
k=1
for(i in seq(2, dim(df_TPM_nonenst)[2], by = 200)){
  df_s<-df_TPM_nonenst[,c(i:min(i+200,dim(df_TPM_nonenst)[2]))]+1
  p<-ggplot(data=stack(df_s), aes(x=ind, y=(values))) + geom_crossbar(stat="summary", fun.y=data_summary, fun.ymax=getlogmax, fun.ymin=getlogmin) +theme(axis.text.x = element_text(angle = 90, hjust = 1))
  plotlist[[k]]=p
  k=k+1
}
pdf("RawTPMPlots_NonEnst.pdf")
for (i in 1:length(plotlist)){
  write(i, stderr())
  print(plotlist[[i]])
}
dev.off()



############normalize###################
gene <- read.csv("transcript_names.txt")
#Define the gene name if your transcript name is different from gene name
tx2gene <- data.frame(GENEID=gene,TXNAME=gene)
#import run ID (SRR)
run <- read.csv("runs.txt",header=F)
IDLIST <- as.character(run$V1)
file <- paste(IDLIST,"",sep="")
names(file) <- run$V1
txi.salmon <- tximport(file, type = "salmon", txOut = TRUE, tx2gene = tx2gene, importer = function(x) read_tsv(x,col_types = list(col_character(), col_double(), col_double(), col_double(), col_double()  )))
