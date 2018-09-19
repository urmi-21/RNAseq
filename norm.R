library(tximport)
library(readr)
library(edgeR)
setwd("~/Desktop/mog_testdata/humandata/sample_named_files/quant_files_named")
#import the transcripts list
#gene <- read.csv("gene.csv")$gene
gene <- read.csv("transcript_names.txt")
#Define the gene name if your transcript name is different from gene name
tx2gene <- data.frame(GENEID=gene,TXNAME=gene)
#import run ID (SRR)
run <- read.csv("runs.txt",header=F)
IDLIST <- as.character(run$V1)
file <- paste(IDLIST,"",sep="")
names(file) <- run$V1
#import salmon file
#gives error now fixed with importer function
txi.salmon <- tximport(file, type = "salmon", txOut = TRUE, tx2gene = tx2gene, importer = function(x) read_tsv(x,col_types = list(col_character(), col_double(), col_double(), col_double(), col_double()  )))
txi.salmon_scaledTPM <- tximport(file, countsFromAbundance="scaledTPM", type = "salmon",txOut = TRUE, tx2gene = tx2gene, importer = function(x) read_tsv(x,col_types = list(col_character(), col_double(), col_double(), col_double(), col_double()  )))
txi.salmon_lenscaledTPM <- tximport(file, countsFromAbundance="lengthScaledTPM", type = "salmon",txOut = TRUE, tx2gene = tx2gene, importer = function(x) read_tsv(x,col_types = list(col_character(), col_double(), col_double(), col_double(), col_double()  )))
df <- txi.salmon$counts
df_scaledTPM <- txi.salmon$counts
df_lenscaledTPM <- txi.salmon$counts

#make boxplots
 
par(mfrow=c(2,1))
boxplot(df,outline=F,xaxt="n",range=0.5)
boxplot(df_scaledTPM,outline=F,xaxt="n",range=0.5)
boxplot(df_lenscaledTPM,outline=F,xaxt="n",range=0.5)
#boxplot(df)
#boxplot(df_scaledTPM)
#boxplot(df_lenscaledTPM)
#import raw counts into edgeR
library(edgeR)
cts <- txi.salmon$counts
normMat <- txi.salmon$length
normMat <- normMat/exp(rowMeans(log(normMat)))
o <- log(calcNormFactors(cts/normMat)) + log(colSums(cts/normMat))
y <- DGEList(cts)
y_norm<-calcNormFactors(y)
y$offset <- t(t(log(normMat)) + o)
boxplot(y$counts,outline=F,xaxt="n",range=0.5)
#after genes are read as DGEList object we should filter very low expressed genes Ref:https://www.bioconductor.org/packages/devel/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf
keep <- rowSums(cpm(y_norm)>1) >= 2
y2 <- y[keep, , keep.lib.sizes=FALSE]
boxplot(y2$counts,outline=F,xaxt="n",range=0.5)

boxplot(y_norm$counts,outline=F,xaxt="n",range=0.5)
##how to make groups of SRRs in a study, e.g. make by studyid
groups<-c("a","B","b")
groups<-c()
cols<-colnames(cts)
for (c in cols){
  #print (c)
  temp <- study_runs_mapping[which(study_runs_mapping$run_accession == c),]$study_accession
  groups<-c(groups,temp)
}
y_grp <- DGEList(cts,group=groups)
boxplot(y_grp$counts,outline=F,xaxt="n",range=0.5)
y_grp_norm<-calcNormFactors(y_grp)
boxplot(y_grp_norm$counts,outline=F,xaxt="n",range=0.5)
cpm <- cpm(y_grp_norm,normalized.lib.sizes=T)
logcpm <- log(cpm+1)

no <-  rowSums(cpm(y)>1) < 2
nocols <- y[no, keep.lib.sizes=FALSE]


cpm <- cpm(y_norm,normalized.lib.sizes=T)
logcpm <- log(cpm+1)
boxplot(logcpm,outline=F,xaxt="n",range=0.5)
boxplot(cpm,outline=F,xaxt="n",range=0.5)

##final edgeR aftr tximport
library(edgeR)
cts <- txi.salmon$counts
normMat <- txi.salmon$length
normMat <- normMat/exp(rowMeans(log(normMat)))
o <- log(calcNormFactors(cts/normMat)) + log(colSums(cts/normMat))
##how to make groups of SRRs in a study, e.g. make by studyid
groups<-c("a","B","b")
groups<-c()
cols<-colnames(cts)
for (c in cols){
  #print (c)
  temp <- study_runs_mapping[which(study_runs_mapping$run_accession == c),]$study_accession
  groups<-c(groups,temp)
}
y_grp <- DGEList(cts,group=groups)
y_norm<-calcNormFactors(y)
y_norm$offset <- t(t(log(normMat)) + o)
boxplot(y_grp_norm$counts,outline=F,xaxt="n",range=0.5)
cpm <- cpm(y_norm,normalized.lib.sizes=T)
logcpm <- log(cpm+1)
boxplot(logcpm,outline=F,xaxt="n",range=0.5)
                                    
 #urmi 14 sept 2018
#Ref: https://bioconductor.org/packages/devel/bioc/vignettes/tximport/inst/doc/tximport.html
# We can avoid gene-level summarization by setting txOut=TRUE, giving the original transcript level estimates as a list of matrices.
txi.salmon <- tximport(file, type = "salmon",countsFromAbundance="scaledTPM", txOut = TRUE)                         
#use importer function for parsing failiures
txi.salmon_scaledTPM <- tximport(file, countsFromAbundance="scaledTPM", type = "salmon",txOut = TRUE, importer = function(x) read_tsv(x,col_types = list(col_character(), col_double(), col_double(), col_double(), col_double() )))
           
