
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
geneMetadata_agg<-aggregate(geneMetadata[,2:20], list(geneMetadata[,1]), function(x) paste0(unique(x),collapse = ";"))

                            
#read hgnc data
hgnc_complete_set <- read_delim("hgnc_complete_set.txt", 
                                "\t", escape_double = FALSE, col_types = cols(alias_symbol = col_skip(),  rna_central_ids = col_skip()), trim_ws = TRUE)
#hgnc columns to keep
toKeep<-c("hgnc_id","symbol","name","locus_group",	"locus_type"	,"gene_family_id","ensembl_gene_id",	"refseq_accession",	"ccds_id",	"uniprot_ids",	"pubmed_id",	"omim_id")
hgnc_complete_set<-hgnc_complete_set[,toKeep]

names(hgnc_complete_set)
names(hgnc_complete_set)[2]<-"hgnc_symbol"
names(allNormalizedTCGA)[1]<-"hgnc_symbol"
joinedDF<-plyr::join(allNormalizedTCGA,hgnc_complete_set,type="left")
#rearrange names
newOrder<-names(joinedDF)[1:2]
newOrder<-c(head(names(joinedDF),2),tail(names(joinedDF),46),names(joinedDF)[3:(length(names(joinedDF))-46)])
joinedDF<-joinedDF[,newOrder]

fwrite(joinedDF,file ="allNormalizedTCGA_hgncdata.csv", row.names = F)
