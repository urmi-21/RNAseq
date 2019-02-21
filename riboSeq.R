#IMPORTANT: Before running make sure that none of the objects exist in current work space.
rm(list=ls())
library(dplyr) #for splitting gds table
library(plyr) #for join
library(readr)
setwd("~/work/urmi/geoscripts")
library(SRAdb)
library(GEOmetadb)
#update sradb file
#getSRAdbFile()
sra_con <- dbConnect(SQLite(), 'SRAmetadb.sqlite')
dbListTables(sra_con)
dbGetQuery(sra_con, paste("select * from metaInfo", sep = " "))
#update geodb
#getSQLiteFile()
congeo <- dbConnect(SQLite(), 'GEOmetadb.sqlite')
dbListTables(congeo)
dbGetQuery(congeo, paste("select * from metaInfo", sep = ""))

#select all data from SRA with Arabidopsis
#taxon_id for Athaliana: 3702
#maize: 4577
#other filters taxon_id in (3702); platform='ILLUMINA';library_strategy='RNA-SEQ';library_layout LIKE 'PAIRED%';library_source='TRANSCRIPTOMIC'
#for human #taxon ID 9606

q1 <-
  dbGetQuery(
    sra_con,
    paste(
      "select * from sra WHERE ( taxon_id IN ('9606') AND  UPPER(platform) = 'ILLUMINA' AND 
      ( UPPER(description) LIKE '%RIBOSOME%' OR UPPER(study_description) LIKE '%RIBOSOME%' OR 
      UPPER(study_abstract) LIKE '%RIBOSOME%' OR UPPER(study_title) LIKE '%RIBOSOME%') )",
      sep = " "
    )
  )

write_tsv(q1,"US_HS_Riboseq.tsv")

q1_filter<- q1 %>% filter(library_source == "OTHER")



q2 <-
  dbGetQuery(
    sra_con,
    paste(
      "select * from sra WHERE ( taxon_id IN ('9606') )",
      sep = " "
    )
  )


qSRA061778<-dbGetQuery(
  sra_con,
  paste(
    "select * from sra WHERE ( experiment_accession = 'SRX740748' )",
    sep = " "
  )
)
