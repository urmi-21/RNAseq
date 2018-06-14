library(tximport)
library(readr)
library(edgeR)
#setwd("~/Desktop/mog_testdata/humandata/sample_named_files/quant_files_named")
setwd("C:/Users/mrbai/Desktop/sample_named_files/quant_files_named")
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
#import study runs mapping
study_runs_mapping <- read_delim("study_runs.mapping.txt", 
                                 "\t", escape_double = FALSE, trim_ws = TRUE)

#import salmon file
#gives error now fixed with importer function
txi.salmon <- tximport(file, type = "salmon", txOut = TRUE, tx2gene = tx2gene, importer = function(x) read_tsv(x,col_types = list(col_character(), col_double(), col_double(), col_double(), col_double()  )))
cts <- txi.salmon$counts
#get tpms
tpms<-txi.salmon$abundance
par(mfrow=c(2,1))
boxplot(cts,outline=F,xaxt="n",range=0.5)
boxplot(tpms,outline=F,xaxt="n",range=0.5)
#txi.salmon_scaledTPM <- tximport(file, countsFromAbundance="scaledTPM", type = "salmon",txOut = TRUE, tx2gene = tx2gene, importer = function(x) read_tsv(x,col_types = list(col_character(), col_double(), col_double(), col_double(), col_double()  )))
#txi.salmon_lenscaledTPM <- tximport(file, countsFromAbundance="lengthScaledTPM", type = "salmon",txOut = TRUE, tx2gene = tx2gene, importer = function(x) read_tsv(x,col_types = list(col_character(), col_double(), col_double(), col_double(), col_double()  )))
#df <- txi.salmon$counts
#df_scaledTPM <- txi.salmon$counts
#df_lenscaledTPM <- txi.salmon$counts
#make boxplots
#par(mfrow=c(2,1))
#boxplot(cts,outline=F,xaxt="n",range=0.5)
#boxplot(df_scaledTPM,outline=F,xaxt="n",range=0.5)
#boxplot(df_lenscaledTPM,outline=F,xaxt="n",range=0.5)
##final edgeR aftr tximport
normMat <- txi.salmon$length
normMat <- normMat/exp(rowMeans(log(normMat)))
o <- log(calcNormFactors(cts/normMat)) + log(colSums(cts/normMat))
##how to make groups of SRRs in a study, e.g. make by studyid
groups<-c()
cols<-colnames(cts)
for (c in cols){
  #print (c)
  temp <- study_runs_mapping[which(study_runs_mapping$run_accession == c),]$study_accession
  groups<-c(groups,temp)
}

#do using tpms
y_tpms <- DGEList(tpms)
y_tpms_norm <- calcNormFactors(y_tpms)
#multiply counts with norfactor
counts<-y_tpms_norm$counts
norms<-y_tpms_norm$samples
cols<-colnames(counts)
counts_df <- as.data.frame(counts)
for(col in cols){
  #print(col)
  #nm1 <- as.symbol(col)
  #counts_df[,(col) := eval(nm1)*100]
  #print(counts_df[,(col)])
  counts_df[,(col)]<-counts_df[,(col)]*(norms[col,]$norm.factors)
}


boxplot(y_tpms_norm$counts,outline=F,xaxt="n",range=0.5)
boxplot(counts_df,outline=F,xaxt="n",range=0.5)

y_grp <- DGEList(cts,group=groups)
y_norm<-calcNormFactors(y_grp)
y_norm_uq<-calcNormFactors(y_grp,method="upperquartile")
y_norm_rle<-calcNormFactors(y_grp,method="RLE")
y_norm$offset <- t(t(log(normMat)) + o)
#after genes are read as DGEList object we should filter very low expressed genes Ref:https://www.bioconductor.org/packages/devel/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf
cpm_cutoff <- 1
keep <- rowSums(cpm(y_norm)>cpm_cutoff) >= 2
removed<- rowSums(cpm(y_norm)>cpm_cutoff) < 2
y_norm_keep <- y_norm[keep, , keep.lib.sizes=FALSE]
y_norm_rem <- y_norm[removed, , keep.lib.sizes=FALSE]
print(paste("removed:",(dim(y_norm_rem)[1]),"transcripts for having low exp"))
boxplot(y_norm$counts,outline=F,xaxt="n",range=0.5)
boxplot(y_norm_keep$counts,outline=F,xaxt="n",range=0.5)
boxplot(y_norm_rem$counts,outline=F,xaxt="n",range=0.5)

par(mfrow=c(2,1))
cpm <- cpm(y_grp,normalized.lib.sizes=T)
logcpm <- log(cpm+1)
cpm_norm <- cpm(y_norm,normalized.lib.sizes=T)
logcpm_norm <- log(cpm_norm+1)
cpm_uq <- cpm(y_norm_uq,normalized.lib.sizes=T)
logcpm_uq <- log(cpm_uq+1)
#these plots are different
boxplot(logcpm,outline=F,xaxt="n",range=0.5)
boxplot(logcpm_norm,outline=F,xaxt="n",range=0.5)

#multiply counts with norfactor
counts<-y_norm$counts
norms<-y_norm$samples
cols<-colnames(counts)
counts_df <- as.data.frame(counts)
for(col in cols){
  #print(col)
  #nm1 <- as.symbol(col)
  #counts_df[,(col) := eval(nm1)*100]
  #print(counts_df[,(col)])
  counts_df[,(col)]<-counts_df[,(col)]*(norms[col,]$norm.factors)
}

boxplot(counts,outline=F,xaxt="n",range=0.5)
boxplot(counts_df,outline=F,xaxt="n",range=0.5)

