library("ArrayExpress", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.5")
library(affy)
library(dplyr)

#query
ATsets = queryAE(species = "Arabidopsis+thaliana")
#keep ones with raw data
ATsetsfiltered<-ATsets %>% filter(tolower(Raw)=="yes")

#download
atSamp = getAE("E-MTAB-4130", type = "raw")

#load downloaded files
atrawset= ae2bioc(mageFiles = atSamp)
atrawset_expr<-exprs(atrawset)
colnames(atrawset_expr)
rownames(atrawset_expr)
atpheno = pData(atrawset)

#read the downloaded file using affy
affy.data = ReadAffy(filenames = atSamp$rawFiles)
#converts an instance of AffyBatch into an instance of ExpressionSet; normalize data
eset.mas5 = mas5(affy.data,normalize = T)
#get expression in matrix
exprSet.nologs = exprs(eset.mas5)
# List the column (chip) names
colnames(exprSet.nologs)
rownames(exprSet.nologs)
dim(exprSet.nologs)
atpheno = pData(eset.mas5)


####using affy http://rstudio-pubs-static.s3.amazonaws.com/14438_729bb8430be242fba106c8ae3968458f.html

files = list.files("celfiles/", full.names = TRUE)
affy.data = ReadAffy(filenames = files)

#converts an instance of AffyBatch into an instance of ExpressionSet; normalize data
eset.mas5 = mas5(affy.data,normalize = F)
#get expression in matrix
exprSet.nologs = exprs(eset.mas5)
# List the column (chip) names
colnames(exprSet.nologs)
rownames(exprSet.nologs)
dim(exprSet.nologs)
heatmap(exprSet.nologs, main = "Normalized ME matrix for brain, liver, N=8")
