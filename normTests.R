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


########################Function defns############################################################

plotAM<-function(twoCtsdf){
 
  names(twoCtsdf)<-c("c1","c2")
  N1<-sum(twoCtsdf[,1])
  N2<-sum(twoCtsdf[,2])
  twoCtsdf$rn<-rownames(ctsdf)
  
  
  
  z1=which(twoCtsdf$c1==0)
  z2=which(twoCtsdf$c2==0)
  z1<-union(z1,z2)
  twoCtsdf<-twoCtsdf[-z1,]
  
  twoCtsdf<-twoCtsdf %>% mutate(M = log2((c1/N1)/(c2/N2)))
  twoCtsdf<-twoCtsdf %>% mutate(A = 0.5*log2((c1/N1)*(c2/N2)))
  
  twoCtsdf<-twoCtsdf %>% mutate(source = ifelse(grepl("ENST", rn),"ENST","NONENST"))
  twoCtsdf<-twoCtsdf %>% group_by(rn) %>% mutate(id= unlist(strsplit(rn, "[.]"))[1])
  twoCtsdf<-twoCtsdf %>% mutate(HS = ifelse(id %in% humanHS$`Transcript stable ID`,"HS", ifelse(source == "ENST" ,"NONHS", "NONENST")))
  
  ggplot(data=twoCtsdf%>%filter(HS=="NONHS"|HS=="NONENST"),aes(x=A,y=M,color=HS))+geom_point(alpha = 0.2)
}

plotnSave<-function(expdf,start,x,fname){
  #expdf df object
  #start start from this column, start=1 if there are no non data columns
  #x num of plots on one page
  #fname filename to save 
  plotlist = list()
  k=1
  for(i in seq(start, dim(expdf)[2], by = x)){
   
    colrange<-c(i:min(i+x,dim(expdf)[2]))
   # colrange
    df_s<-expdf[,c(colrange)]+1
    meanlines<-log((as.data.frame((apply((df_s),MARGIN=2,FUN=mean)))))
    names(meanlines)<-c("value")
    meanlines[is.na(meanlines)] <- 0
    meanlines<-meanlines%>%mutate(num=row_number(value))
    meanlines$num<-as.numeric(rownames(meanlines))
    #3rd qt
    q3Lines<-log(as.data.frame((apply((df_s),MARGIN=2,FUN=qt3))))
    names(q3Lines)<-c("value")
    q3Lines[is.na(q3Lines)] <- 0
    q3Lines<-q3Lines%>%mutate(num=row_number(value))
    q3Lines$num<-as.numeric(rownames(q3Lines))
    
    p<-ggplot(data=stack(df_s), aes(x=ind, y=(values))) + geom_crossbar(stat="summary", fun.y=data_summary, fun.ymax=getlogmax, fun.ymin=getlogmin)+geom_segment(data=meanlines,aes(x=num-0.45,xend=num+0.45,y=value,yend=value),inherit.aes=FALSE,color="Red",size=1.5)+geom_segment(data=q3Lines,aes(x=num-0.45,xend=num+0.45,y=value,yend=value),inherit.aes=FALSE,color="Green",size=1.5)+scale_y_continuous(breaks=seq(0,20,1))+theme(axis.text.x = element_text(angle = 90, hjust = 1))
    
    plotlist[[k]]=p
    k=k+1
  }
  
  pdf(fname)
  for (i in 1:length(plotlist)){
    write(i, stderr())
    print(plotlist[[i]])
  }
  dev.off()
}

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
qt3<-function(x){
  return(quantile(x,na.rm=T)[4])
}

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

########################End function defn############################################################

##################################Read raw data files##############################
ind<<-1
fnames<-as.data.frame(list.files(full.names = TRUE,pattern = "SRR*"))
names(fnames)<-c("fname")
fnames<-fnames %>% mutate(fname=str_replace_all(fname, "[^[:alnum:]]", ""))

df <- list.files(full.names = TRUE,pattern = "SRR*") %>% lapply(read_tsv) %>% lapply(changeNames) %>% reduce(cbind)

tpmCols<- names(df)[grepl("Name|TPM",names(df))]
df_TPM<-df[,tpmCols]


##################Plot test#############################################################
newdf<-df_TPM[,2:4]+1

ggplot(data=stack(newdf), aes(x=ind, y=(values))) + geom_crossbar(stat="summary", fun.y=data_summary, fun.ymax=getlogmax, fun.ymin=getlogmin,inherit.aes = TRUE)+ geom_hline(aes(yintercept = log(mean(values))),size=1)  

#toplot meanlines
meanlines<-log((as.data.frame((apply((df_s),MARGIN=2,FUN=mean)))))
names(meanlines)<-c("value")
meanlines[is.na(meanlines)] <- 0
meanlines<-meanlines%>%mutate(num=row_number(value))
meanlines$num<-as.numeric(rownames(meanlines))
#3rd qt
q3Lines<-log(as.data.frame((apply((newdf),MARGIN=2,FUN=qt3))))
names(q3Lines)<-c("value")
q3Lines[is.na(q3Lines)] <- 0
q3Lines<-q3Lines%>%mutate(num=row_number(value))
q3Lines$num<-as.numeric(rownames(q3Lines))
ggplot(data=stack(newdf), aes(x=ind, y=(values))) + geom_crossbar(stat="summary", fun.y=data_summary, fun.ymax=getlogmax, fun.ymin=getlogmin)+geom_segment(data=meanlines,aes(x=num-0.45,xend=num+0.45,y=value,yend=value),inherit.aes=FALSE,color="Red",size=1.5)+geom_segment(data=q3Lines,aes(x=num-0.45,xend=num+0.45,y=value,yend=value),inherit.aes=FALSE,color="Green",size=1.5)+scale_y_continuous(breaks=seq(0,20,1))
######################################################################################

##################Plot to File#################################

#sort by depth and plot,read bases info
bases <- read_delim("bases.txt", "\t", escape_double = FALSE,trim_ws = TRUE)
bases<-bases %>% mutate(srr=paste("TPM_",run_accession,sep = ""))
bases <- bases[bases$srr %in% names(df_TPM),]
#bases <- bases[bases$srr %in% names(ctsdf),]
bases <- bases[order(bases$bases),]

df_TPM<-df_TPM[,c('Name',bases$srr)]
#plot
plotnSave(df_TPM_sorted,2,200,"RawTPMPSortedBase.pdf")
plotnSave(df_TPM_sorted,2,50,"RawTPMPSortedBase_50.pdf")
###separate tpms of Ens genes and others
enstNames<-df_TPM$Name[(grepl("ENST",df_TPM$Name))]
df_TPM_enst<-df_TPM[df_TPM$Name %in% enstNames,]
df_TPM_nonenst<-df_TPM[!(df_TPM$Name %in% enstNames),]
#sort by bases
df_TPM_enst<-df_TPM_enst[,c('Name',bases$srr)]
df_TPM_nonenst<-df_TPM_nonenst[,c('Name',bases$srr)]
plotnSave(df_TPM_enst,2,200,"RawTPMPlots_Enst.pdf")
plotnSave(df_TPM_enst,2,50,"RawTPMPlots_Enst_50.pdf")
#plot non enst
plotnSave(df_TPM_nonenst,2,200,"RawTPMPlots_NonEnst.pdf")
plotnSave(df_TPM_nonenst,2,50,"RawTPMPlots_NonEnst_50.pdf")


#################################################################normalize##################################################################
gene <- read.csv("transcript_names.txt")
#Define the gene name if your transcript name is different from gene name
tx2gene <- data.frame(GENEID=gene,TXNAME=gene)
#import run ID (SRR)
run <- read.csv("runs.txt",header=F)
#remove runs with na values
run<- as.data.frame(run$V1[ !(run$V1 %in% c("SRR363862","SRR2723895"))])
names(run)<-c("V1")
IDLIST <- as.character(run$V1)
file <- paste(IDLIST,"",sep="")
names(file) <- run$V1
txi.salmon <- tximport(file, type = "salmon", txOut = TRUE, tx2gene = tx2gene, importer = function(x) read_tsv(x,col_types = list(col_character(), col_double(), col_double(), col_double(), col_double()  )))

##using egder
cts <- txi.salmon$counts
lengths<-txi.salmon$length
length(rownames(cts))
length(unique(rownames(cts)))
#from cts aggregate repeated value of ENST00000378936
ctsNew<- as.matrix(t(colSums(cts[rownames(cts)%in% c("ENST00000378936"),])))
rownames(ctsNew)<-c("ENST00000378936")
cts_removed<-cts[!(rownames(cts)%in% c("ENST00000378936")),]
cts<-rbind(cts_removed,ctsNew)
cts_removed<-NULL
ctsNew<-NULL
#remove extra len frm lengths
lensNew<- as.matrix(t(colSums(lengths[rownames(lengths)%in% c("ENST00000378936"),])))
rownames(lensNew)<-c("ENST00000378936")
lens_removed<-lengths[!(rownames(lengths)%in% c("ENST00000378936")),]
lengths<-rbind(lens_removed,lensNew)
lens_removed<-NULL
lensNew<-NULL
dim(lengths)

#plot dist of raw counts by bases
ctsdf<-as.data.frame(cts)

#remove excluded cols from bases
#if bases not in memory
bases <- read_delim("bases.txt", "\t", escape_double = FALSE,trim_ws = TRUE)
bases<-bases %>% mutate(srr=paste("TPM_",run_accession,sep = ""))
bases <- bases[bases$run_accession %in% names(ctsdf),]
bases <- bases[order(bases$bases),]
basesRemoved<-bases[which(!(bases$run_accession %in% c("SRR363862","SRR2723895"))),]
#sort by bases
ctsdf<-ctsdf[,c(basesRemoved$run_accession)]

#plot estimated counts
plotnSave(ctsdf,1,200,"Counts.pdf")
plotnSave(ctsdf,1,50,"Counts_50.pdf")

enstNames<-rownames(ctsdf)[(grepl("ENST",rownames(ctsdf)))]
nonenstNames<-rownames(ctsdf)[(grepl("ENST",rownames(ctsdf)))]
###################################################################################################

###########Start normalization#####################################################################
normMat <- lengths
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
snf<-calcNormFactors(sampexp,refColumn=1,logratioTrim=0.5,sumTrim=0.5,doWeighting = TRUE, Acutoff = -1e+10)
sampexp_norm<-sampexp %*% diag(sqrt(snf))
par(mfrow=c(1,2))
boxplot(sampexp)
boxplot(sampexp_norm)
par(mfrow=c(1,1))




#using noiseq
head(cts)
myfactors<-run
#read gene lengths
mylength <- read_delim("SRR073758", "\t", escape_double = FALSE, col_types = cols(EffectiveLength = col_skip(), NumReads = col_skip(), TPM = col_skip()), 
                        trim_ws = TRUE)
ml<-as.data.frame(t(mylength[,-1]))
colnames(ml)<-mylength$Name
head(ml[,1:4])


nsData<-readData(data = cts,factors = myfactors)

myTMM <- tmm(assayData(nsData)$exprs, long = ml, lc = 0)
myTMM_default <- tmm(assayData(nsData)$exprs, long = 1000, lc = 0, k = 0, refColumn = 1, logratioTrim = 0.3, sumTrim = 0.05, doWeighting = TRUE, Acutoff = -1e+10)
myUQUA = uqua(assayData(nsData)$exprs)
myRPKM = rpkm(assayData(nsData)$exprs, long = ml, k = 0, lc = 1)

dim(myTMM)
dim(myTMM_default)
dim(myUQUA)
dim(myRPKM)
head(myTMM,2)
head(myUQUA,2)
head(myRPKM,2)



#plot graphs
#plot non enst
thisDF<-as.data.frame(myTMM_default)
#sort by bases
basesRemoved<-bases$run_accession[which(!(bases$run_accession %in% c("SRR363862","SRR2723895")))]
thisDF<-thisDF[,c(basesRemoved)]

plotnSave(thisDF,1,200,"TMMdefaulf.pdf")
plotnSave(thisDF,1,50,"TMMdefaulf_50.pdf")

#plot enst and other
thisDF_enst<-thisDF[rownames(thisDF) %in% enstNames,]
thisDF_nonenst<-thisDF[!(rownames(thisDF) %in% enstNames),]
plotnSave(thisDF_enst,1,200,"TMMdefaulf_ENST.pdf")
plotnSave(thisDF_enst,1,50,"TMMdefaulf_ENST_50.pdf")
plotnSave(thisDF_nonenst,1,200,"TMMdefaulf_nonENST.pdf")
plotnSave(thisDF_nonenst,1,50,"TMMdefaulf_nonENST_50.pdf")

log(mean(dthisDF_nonenst$SRR2735916+1))
log(mean(dthisDF_enst$SRR2121909))
log(mean(df_s$SRR2121909))



#plot A vs M values
plotAM(ctsdf[,190:191])
plotAM(ctsdf[,c(19,71)])






#read human HS genes/tx names
humanHS <- read_delim("humanHS.txt", "\t",escape_double = FALSE, trim_ws = TRUE)

