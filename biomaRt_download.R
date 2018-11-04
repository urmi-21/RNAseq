
#if (!requireNamespace("BiocManager", quietly = TRUE)){
#  install.packages("BiocManager")
#}
#  
#BiocManager::install("biomaRt", version = "3.8")
library(readr)
library(biomaRt)
allNormalizedTCGA <- read_csv("allNormalizedTCGA.csv")


geneSymbols<-allNormalizedTCGA$Hugo_Symbol
write.csv(geneSymbols,"genesSymbols.csv",row.names=F)

mart = useMart("ensembl", dataset="hsapiens_gene_ensembl")
#see attribute and filter list
attlist<-listAttributes(mart = mart)
filterlist<-listFilters(mart)

go=c("GO:0051330","GO:0000080","GO:0000114","GO:0000082")
chrom=c(17,20,"Y")
#list of attributes to get back
myAttlist<-c("ensembl_gene_id","description","chromosome_name","start_position","end_position","strand","transcript_count","percentage_gene_gc_content","gene_biotype","source","go_id","name_1006","definition_1006","go_linkage_type","namespace_1003","goslim_goa_accession","goslim_goa_description","arrayexpress","hgnc_symbol","kegg_enzyme")
geneMetadata<-getBM(attributes= myAttlist,filters=c("hgnc_symbol"),values=list(geneSymbols), mart=mart)

#aggregate repeated cols
geneMetadata_agg<-aggregate(geneMetadata[,2:20], list(geneMetadata[,1]), function(x) paste0(unique(x)))
