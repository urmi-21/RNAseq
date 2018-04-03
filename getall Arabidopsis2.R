#IMPORTANT: Before running make sure that none of the objects exist in current work space.
rm(list=ls())
library(tidyr) #for splitting gds table
library(plyr) #for join
setwd("~/work/urmi/geoscripts")
library(SRAdb)
#update sradb file
#getSRAdbFile()
sra_con <- dbConnect(SQLite(), 'SRAmetadb.sqlite')
dbListTables(sra_con)
dbGetQuery(sra_con, paste("select * from metaInfo", sep = " "))
library(GEOmetadb)
#getSQLiteFile()#update db
congeo <- dbConnect(SQLite(), 'GEOmetadb.sqlite')
dbListTables(congeo)
dbGetQuery(congeo, paste("select * from metaInfo", sep = ""))

#select all data from SRA with Arabidopsis
#taxon_id for Athaliana: 3702
#maize: 4577
#other filters taxon_id in (3702); platform='ILLUMINA';library_strategy='RNA-SEQ';library_layout LIKE 'PAIRED%';library_source='TRANSCRIPTOMIC'
#for human #taxon ID 9606
all_sra_data <-
  dbGetQuery(
    sra_con,
    paste(
      "select * from sra WHERE (taxon_id IN ('3702') AND UPPER(platform) = 'ILLUMINA' AND UPPER(library_strategy) = 'RNA-SEQ' AND UPPER(library_layout) LIKE 'PAIRED%' AND UPPER(library_source) = 'TRANSCRIPTOMIC')",
      sep = " "
    )
  )

#all_sra_data <-
#  dbGetQuery(
#    sra_con,
#    paste(
#      "select * from sra WHERE (taxon_id IN ('4577','381124','112001','381124','334825','4579','76912','1980702','368615') AND UPPER(platform) = 'ILLUMINA' AND UPPER(library_strategy) = 'RNA-SEQ' AND UPPER(library_layout) LIKE 'PAIRED%' AND UPPER(library_source) = 'TRANSCRIPTOMIC')",
#      sep = " "
#    )
#  )

#retain only GSE ids
gse_in_sra <- c()
gsm_in_sra <- c()
#add new column to results with clean GSM ids this will be used to join data later
#Gsm can be in run_alias or sample_alias
for (i in 1:length(all_sra_data$run_alias)) {
  
  if (grepl("GSE", all_sra_data$study_alias[i])) {
    
    if (grepl("GSM", all_sra_data$run_alias[i])) {
      s <- all_sra_data$run_alias[i]
      all_sra_data$gsm[i] <-
        unlist(strsplit(s, split = '_', fixed = TRUE))[1]
      gsm_in_sra<-c(gsm_in_sra,all_sra_data$gsm[i])
      all_sra_data$gse[i]<-all_sra_data$study_alias[i]
      gse_in_sra<-c(gse_in_sra,all_sra_data$study_alias[i])    
      
    }
    else if (grepl("GSM", all_sra_data$sample_alias[i])) {
      s <- all_sra_data$sample_alias[i]
      all_sra_data$gsm[i] <-
        unlist(strsplit(s, split = '_', fixed = TRUE))[1]
      gsm_in_sra<-c(gsm_in_sra,all_sra_data$gsm[i])
      all_sra_data$gse[i]<-all_sra_data$study_alias[i]
      gse_in_sra<-c(gse_in_sra,all_sra_data$study_alias[i])
    }
  }
  else{
    #print(all_sra_data$run_alias[i])
    all_sra_data$gsm[i] <- NA
    all_sra_data$gse[i]<-NA
  }
  
  
}
#keep only uniqe values
gsm_in_sra <-unique(gsm_in_sra)
gse_in_sra <- unique(gse_in_sra)
my_gse_gsm <- as.data.frame(all_sra_data[,c("gse","gsm")])
#remove duplicate
my_gse_gsm <- unique(my_gse_gsm)
#remove NAs
my_gse_gsm<-my_gse_gsm[rowSums(is.na(my_gse_gsm)) == 0,]
#Many entries in SRA have GSE id in study but no coressponding sample ids
#i will ignore such rows

#Get all data from geo where gse is in sra results
all_gse <-
  dbGetQuery(congeo,
             paste(
               "select * from gse WHERE gse IN(",
               paste(shQuote(gse_in_sra), collapse = ", "),
               ")",
               sep = " "
             ))

#change col names in gse
for (i in 1:length(names(all_gse))) {
  #don't change is col name is gse
  if (colnames(all_gse)[i] == "gse") {
    colnames(all_gse)[i] <- paste("gse")
  }
  else{
    colnames(all_gse)[i] <- paste("gse.", colnames(all_gse)[i],sep="")
  }
}

#get all GSM which has gsm ids found in sra data
all_gsm <-
  dbGetQuery(congeo,
            paste(
               "select * from gsm WHERE gsm IN(",
               paste(shQuote(gsm_in_sra), collapse = ", "),
               ") AND UPPER(molecule_ch1) LIKE ('%POLYA%')",
               sep = " "
             ))

all_gsm <-
  dbGetQuery(congeo,
             paste(
               "select * from gsm WHERE gsm IN(",
               paste(shQuote(gsm_in_sra), collapse = ", "),
               ")",
               sep = " "
             ))

#change col names in gse
for (i in 1:length(names(all_gsm))) {
  #don't change is col name is gsm
  if (colnames(all_gsm)[i] == "gsm") {
    colnames(all_gsm)[i] <- paste("gsm")
  }
  else{
    colnames(all_gsm)[i] <- paste("gsm.", colnames(all_gsm)[i],sep="")
  }
}

#step1 join gsm data with all_sra_data
#we want inner join keep only those rows which are in both sra and geo and keep the filters
#join_sra_gsm <- plyr::join(all_sra_data, all_gsm,by="gsm",type ="inner")
#join type default (left) keep all sra data even with no match in geo
join_sra_gsm <- plyr::join(all_sra_data, all_gsm,by="gsm")

#step2 join gse data with all_sra_data
#join_sra_gsm_gse <- plyr::join(join_sra_gsm, all_gse,by="gse",type ="inner")
join_sra_gsm_gse <- plyr::join(join_sra_gsm, all_gse,by="gse")

##replace all \t in each column of join_sra_gsm_gse
for (col in colnames(join_sra_gsm_gse)){
  print(col)
  for (j in 1:length(join_sra_gsm_gse[[col]])){
    #print (join_sra_gsm_gse[j,col])
    join_sra_gsm_gse[j,col]=gsub("\t",";",join_sra_gsm_gse[j,col])
  }
    
}

#filter out unwanted columns
keep <-
  c(
    "study_accession",
    "study_title",
    "study_abstract",
    
    
    "sample_accession",		
    "description",
    "sample_attribute",
    
    "experiment_accession",
    "experiment_title",
    "design_description",
    "library_name",
    "library_selection",
    "library_layout",
    "library_construction_protocol",
    "instrument_model",
    "experiment_attribute",
    
    
    "run_accession",
    "bases",
    
    
    "gse",
    "gse.title",
    "gse.summary",
    "gse.type",
    "gse.overall_design",
    
    
    "gsm",
    "gsm.source_name_ch1",
    "gsm.characteristics_ch1",
    "gsm.molecule_ch1",
    "gsm.treatment_protocol_ch1",
    "gsm.extract_protocol_ch1",
    "gsm.description"
    
  )


filterfinal<-join_sra_gsm_gse[,keep]
#remove . from col names otherwise MOG wont read because of No2 library
this_cols<-colnames(filterfinal)
for (j in 1:length(this_cols)){
  #print (join_sra_gsm_gse[j,col])
  this_cols[j]=gsub('[.]','_',this_cols[j],fixed = F)
}
colnames(filterfinal)<-this_cols
#write to csv file
write.csv(filterfinal, "US_AT_removedtabs.csv", row.names=FALSE)
write_tsv(filterfinal, "US_AT_removedtabs.txt")
#write.table(filterfinal, "final.csv", row.names=FALSE)

#debug
rw<-dbGetQuery(sra_con, paste("select * from sra WHERE run_accession in ('ERR1407268')", sep=" "))
##below is a problem
#gse_gsm table doesn't map all gse_gsm belo gsm isn't in gse_gsm fin a solution to this.

rs<-dbGetQuery(congeo, paste("select * from gsm WHERE gsm in ('GSM1868811')", sep=" "))
rs2<-dbGetQuery(congeo, paste("select * from gse_gsm WHERE gsm in ('GSM1868811')", sep=" "))

