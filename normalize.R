library(tximport)
library(readr)
library(edgeR)
#import the transcripts list
#gene <- read.csv("gene.csv")$gene
gene <- read.csv("transcript_names.txt")
#Define the gene name if your transcript name is different from gene name
tx2gene <- data.frame(GENEID=gene,TXNAME=gene)
#import run ID (SRR)


run <- read.csv("runs.txt",header=F)
#library(readr)
#runids <- read_delim("~/Desktop/mog_testdata/humandata/runids.txt", 
#                    "\t", escape_double = FALSE, col_names = FALSE, 
 #                   trim_ws = TRUE)
#run <- t(runids[-2,])
#run <- as.data.frame(run)
IDLIST <- as.character(run$V1)
file <- paste(IDLIST,"",sep="")
names(file) <- run$V1
#import the output from kallisto. It include raw counts, efficient length, and tpm 
#txi.kallisto.tsv <- tximport(file, type = "kallisto", tx2gene = tx2gene)
#import salmon file
txi.salmon <- tximport(file, type = "salmon", tx2gene = tx2gene)
txi.salmon <- tximport("quant.sf", type = "salmon", tx2gene = tx2gene)
#df <- txi.kallisto.tsv$counts
df <- txi.salmon$counts

#import raw counts into edgeR
y <- DGEList(df)
#may have runs without any expression, just discard them
keep <- y$samples$lib.size != 0
y <- y[ ,keep]
no <-  y$samples$lib.size == 0
colnames(y[,no])
#calculate TMM scaling factor, after that, the matrix is still raw count, it just add a column of scaling factor in the y$sample
y <- calcNormFactors(y)
#calculate cpm, the matrix is after TMM and cpm
cpm <- cpm(y,normalized.lib.sizes=T)
logcpm <- log(cpm+1)
#boxplot
png(filename="anno_del_normalizaed.png",width=3000,height=1000,type="cairo")
boxplot(logcpm,outline=F,xaxt="n",range=0.5)
boxplot(cpm,outline=F,xaxt="n",range=0.5)
par(mfrow=c(4,1))
boxplot(logcpm[1:6692,1:864],outline=F,xaxt="n",range=0.5)
boxplot(logcpm[1:6692,865:1728],outline=F,xaxt="n",range=0.5)
boxplot(logcpm[1:6692,1729:2592],outline=F,xaxt="n",range=0.5)
boxplot(logcpm[1:6692,2593:3457],outline=F,xaxt="n",range=0.5)
dev.off()
#expression plot
cpm <- as.matrix(cpm)
quantile(cpm,probs = seq(0,1,0.1))
png(filename="heatmap.png",width=3000,height=1500,type="cairo")
image(t(cpm2[36046:1,]),col = rev(heat.colors(4)), breaks=c(0.000000001,0.309,2.712,18.030,1.13e7),axes=F, main="Expression level for all CDF (cpm)", cex.main=2)
dev.off()
