library(readr)
library(magrittr)
library(dplyr)
library(purrr)
setwd("~/Downloads/RNAseqDB/RNAseqDB/data/expected_count")

brca_rsem_fpkm_tcga_t <- read_delim("~/Downloads/RNAseqDB/RNAseqDB/data/normalized/brca-rsem-fpkm-tcga-t.txt","\t", escape_double = FALSE, trim_ws = TRUE)
brca_rsem_fpkm_tcga<- read_delim("~/Downloads/RNAseqDB/RNAseqDB/data/normalized/brca-rsem-fpkm-tcga.txt","\t", escape_double = FALSE, trim_ws = TRUE)
brca_gtex <- read_delim("~/Downloads/RNAseqDB/RNAseqDB/data/normalized/breast-rsem-fpkm-gtex.txt","\t", escape_double = FALSE, trim_ws = TRUE)
joined<- plyr::join(brca_rsem_fpkm_tcga_t,brca_rsem_fpkm_tcga)
joined2<-plyr::join(joined,brca_gtex)

write.csv(joined2,"brca_combined.csv",row.names = F)

cerv<-read_delim("~/Downloads/RNAseqDB/RNAseqDB/data/normalized/cervix-rsem-fpkm-gtex.txt","\t", escape_double = FALSE, trim_ws = TRUE)

joined3<-plyr::join(joined2,cerv)

write.csv(joined3,"brca_cerv_combined.csv",row.names = F)

setwd("~/Downloads/RNAseqDB/RNAseqDB/data/normalized/data_files")

#df <- list.files(full.names = TRUE) %>% lapply(read_tsv)  %>% reduce(full_join)



#plot hist of all the samples
df <- list.files(full.names = TRUE,pattern = "*.txt") %>% lapply(read_tsv) %>% reduce(inner_join)
which(is.na(df))

##remove duplicated columns. these cols are ambiguous. first two matches are false
toremove <- names(df)[grepl("_",names(df))]
toremove<-toremove[3:length(toremove)]
tr<-c()
for(s in toremove){
  tr<-c(tr,unlist(strsplit(s,"_"))[1])
}
tr<-unique(tr)
df2<-df[,!(names(df) %in% toremove)]
df2<-df2[,!(names(df2) %in% tr)]

##remove GTEX columns
toremove <- names(df2)[grepl("GTEX",names(df2))]
df2<-df[,!(names(df) %in% toremove)]

##only GTEX data

dfgtex<-df[,(names(df) %in% c(names(df[1:2]),toremove))]

#write.csv(df,"allCombined.csv",row.names = F)
#faster
library("data.table")
fwrite(dfgtex,file ="allGTEXExpectedCts.csv", row.names = F)


#plot hist of all the samples
df <- list.files(full.names = TRUE,pattern = "*.txt") %>% lapply(read_tsv) %>% reduce(inner_join)
which(is.na(df))

#remove info cols
df2<-df[,c(3:dim(df)[2])]
#order names
df2<-df2[,order(names(df2))]
#ggplot(stack(df_s), aes(x = ind, y = log(values))) +  geom_boxplot()+theme(axis.text.x = element_text(angle = 90, hjust = 1))


plotlist = list()
k=1
for(i in seq(3, dim(df2)[2], by = 50)){
  df_s<-df2[,c(i:min(i+50,dim(df2)[2]))]
  #print(colnames(df[,c(i:min(i+50,dim(df)[2]))]))
  p<-ggplot(stack(df_s), aes(x = ind, y = log(values))) +  geom_boxplot()+theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  plotlist[[k]]=p
  k=k+1
}

pdf("countsPlots.pdf")
for (i in 1:length(plotlist)){
  write(i, stderr())
  print(plotlist[[i]])
}
dev.off()

