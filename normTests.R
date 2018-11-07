library(tximport)
library(readr)
library(edgeR)
library(tidyverse)
library(magrittr)
library(ggplot2)

library('seluth')
library('wasabi')
library("NOISeq", lib.loc="~/R/win-library/3.5")

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
  if(is.na(mu)){
    mu<-1
  }
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

ggplot(data=stack(newdf), aes(x=ind, y=(values))) + geom_crossbar(stat="summary", fun.y=data_summary, fun.ymax=getlogmax, fun.ymin=getlogmin,inherit.aes = TRUE)+ geom_hline(aes(yintercept = log(mean(values))),size=1)  

#toplot meanlines
meanlines<-log(as.data.frame((apply((newdf),MARGIN=2,FUN=mean))))
names(meanlines)<-c("value")
meanlines<-meanlines%>%mutate(num=row_number(value))

ggplot(data=stack(newdf), aes(x=ind, y=(values))) + geom_crossbar(stat="summary", fun.y=data_summary, fun.ymax=getlogmax, fun.ymin=getlogmin)+geom_segment(data=meanlines,aes(x=num-0.45,xend=num+0.45,y=value,yend=value),inherit.aes=FALSE,color="Red",size=1.5)

ggplot(data=stack(newdf), aes(x=ind, y=log(values))) +   geom_violin() 

##################Plot to File#################################

plotlist = list()
k=1
for(i in seq(2, dim(df_TPM)[2], by = 200)){
  df_s<-df_TPM[,c(i:min(i+200,dim(df_TPM)[2]))]+1
  #toplot meanlines
  meanlines<-log((as.data.frame((apply((df_s),MARGIN=2,FUN=mean)))))
  names(meanlines)<-c("value")
  meanlines[is.na(meanlines)] <- 0
  meanlines<-meanlines%>%mutate(num=row_number(value))
  meanlines$num<-as.numeric(rownames(meanlines))
  p<-ggplot(data=stack(df_s), aes(x=ind, y=(values))) + geom_crossbar(stat="summary", fun.y=data_summary, fun.ymax=getlogmax, fun.ymin=getlogmin)+geom_segment(data=meanlines,aes(x=num-0.25,xend=num+0.25,y=value,yend=value),inherit.aes=FALSE,color="Red",size=1.5) +theme(axis.text.x = element_text(angle = 90, hjust = 1))
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


df_TPM_sorted<-df_TPM[,c('Name',bases$srr)]

plotlist = list()
k=1
for(i in seq(2, dim(df_TPM_sorted)[2], by = 200)){
  df_s<-df_TPM_sorted[,c(i:min(i+200,dim(df_TPM_sorted)[2]))]+1 #add 1 to avoid inf
  meanlines<-log((as.data.frame((apply((df_s),MARGIN=2,FUN=mean)))))
  names(meanlines)<-c("value")
  meanlines[is.na(meanlines)] <- 0
  meanlines<-meanlines%>%mutate(num=row_number(value))
  meanlines$num<-as.numeric(rownames(meanlines))
  p<-ggplot(data=stack(df_s), aes(x=ind, y=(values))) + geom_crossbar(stat="summary", fun.y=data_summary, fun.ymax=getlogmax, fun.ymin=getlogmin)+geom_segment(data=meanlines,aes(x=num-0.25,xend=num+0.25,y=value,yend=value),inherit.aes=FALSE,color="Red",size=1.5) +theme(axis.text.x = element_text(angle = 90, hjust = 1))
  plotlist[[k]]=p
  k=k+1
}

pdf("RawTPMPSortedBase.pdf")
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
for(i in seq(2, dim(df_TPM_enst)[2], by = 50)){
  df_s<-df_TPM_enst[,c(i:min(i+50,dim(df_TPM_enst)[2]))]+1
  meanlines<-log((as.data.frame((apply((df_s),MARGIN=2,FUN=mean)))))
  names(meanlines)<-c("value")
  meanlines[is.na(meanlines)] <- 0
  meanlines<-meanlines%>%mutate(num=row_number(value))
  meanlines$num<-as.numeric(rownames(meanlines))
  p<-ggplot(data=stack(df_s), aes(x=ind, y=(values))) + geom_crossbar(stat="summary", fun.y=data_summary, fun.ymax=getlogmax, fun.ymin=getlogmin)+geom_segment(data=meanlines,aes(x=num-0.25,xend=num+0.25,y=value,yend=value),inherit.aes=FALSE,color="Red",size=1.5) +theme(axis.text.x = element_text(angle = 90, hjust = 1))
  plotlist[[k]]=p
  k=k+1
}
pdf("RawTPMPlots_Enst_50.pdf")
for (i in 1:length(plotlist)){
  write(i, stderr())
  print(plotlist[[i]])
}
dev.off()

#plot non enst
plotlist = list()
k=1
for(i in seq(2, dim(df_TPM_nonenst)[2], by = 50)){
  df_s<-df_TPM_nonenst[,c(i:min(i+50,dim(df_TPM_nonenst)[2]))]+1
  meanlines<-log((as.data.frame((apply((df_s),MARGIN=2,FUN=mean)))))
  names(meanlines)<-c("value")
  meanlines[is.na(meanlines)] <- 0
  meanlines<-meanlines%>%mutate(num=row_number(value))
  meanlines$num<-as.numeric(rownames(meanlines))
  p<-ggplot(data=stack(df_s), aes(x=ind, y=(values))) + geom_crossbar(stat="summary", fun.y=data_summary, fun.ymax=getlogmax, fun.ymin=getlogmin)+geom_segment(data=meanlines,aes(x=num-0.25,xend=num+0.25,y=value,yend=value),inherit.aes=FALSE,color="Red",size=1.5) +theme(axis.text.x = element_text(angle = 90, hjust = 1))
  plotlist[[k]]=p
  k=k+1
}
pdf("RawTPMPlots_NonEnst_50.pdf")
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
#remove with na values
run<- as.data.frame(run$V1[ !(run$V1 %in% c("SRR363862","SRR363862"))])
names(run)<-c("V1")
IDLIST <- as.character(run$V1)
file <- paste(IDLIST,"",sep="")
names(file) <- run$V1
txi.salmon <- tximport(file, type = "salmon", txOut = TRUE, tx2gene = tx2gene, importer = function(x) read_tsv(x,col_types = list(col_character(), col_double(), col_double(), col_double(), col_double()  )))

##using egder
cts <- txi.salmon$counts
normMat <- txi.salmon$length
normMat <- normMat/exp(rowMeans(log(normMat)))

o <- log(calcNormFactors(cts/normMat)) + log(colSums(cts/normMat))
y <- DGEList(cts)
y <- scaleOffset(y, t(t(log(normMat)) + o))
#store srr by lib size
libsizesOrderd<-sm[order(sm$lib.size),]
# filtering
keep <- filterByExpr(y)
y <- y[keep, ]
# y is now ready for estimate dispersion functions see edgeR User's Guide
y$samples

yNF <- calcNormFactors(y)
head(yNF$samples)
head(y$samples)
yNFwzp <- calcNormFactors(y,method="TMMwzp")
head(yNFwzp$samples)
yNFuq <- calcNormFactors(y,method="upperquartile")
head(yNFuq$samples)

##example of norffactors
sampexp <- matrix( rpois(10000, lambda=5), nrow=100 )
snf<-calcNormFactors(sampexp)
sampexp_norm<-sampexp %*% diag(snf)
boxplot(sampexp)
boxplot(sampexp_norm)





#using noiseq
head(cts)
myfactors<-run
nsData<-readData(data = cts,factors = myfactors)
#### ENST00000378936 is duplicated##########
length(rownames(cts))
length(unique(rownames(cts)))
#sum the repeated transcript
t(sapply(by(cts,rownames(cts),colSums),identity))
















