library(readr)
library(magrittr)
library(dplyr)
library(purrr)
#join all data files
setwd("~/Downloads/RNAseqDB/RNAseqDB/data/normalized/data_files")
#full join to keep all the genes
df <- list.files(full.names = TRUE,pattern = "*.txt") %>% lapply(read_tsv) %>% reduce(full_join)
which(is.na(df))


(as.character(df[which(df$Hugo_Symbol=="GJA10"),]))

df[which(df$Hugo_Symbol=="SQLE"),"GTEX-WI4N-2426-SM-4OOSC"]

#replace NA with 0
df[is.na(df)] <- 0
which(is.na(df))

#read ensembl info and add to df
ens_mart_genes <- read_delim("~/Downloads/RNAseqDB/RNAseqDB/data/normalized/data_files/ens_mart_genes.txt", 
                             "\t", escape_double = FALSE, trim_ws = TRUE)

ens_mart_genes_clean<-aggregate(.~`Gene stable ID`,ens_mart_genes, function(x) toString(unique(x)))
