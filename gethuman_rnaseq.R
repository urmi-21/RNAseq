setwd("~/work/urmi/geoscripts")
library(SRAdb)
sra_con <- dbConnect(SQLite(),'SRAmetadb.sqlite')
dbListTables(sra_con)
library(GEOmetadb)
congeo <- dbConnect(SQLite(),'GEOmetadb.sqlite')
rw <- dbGetQuery(congeo, paste("select * from gse LIMIT 5",sep=""))


# #get human RNA-seq from SRA
# sra_data_human_rnaseq <-dbGetQuery(sra_con, paste("select s.*, e.*, m.* from study s INNER JOIN experiment e ON s.study_accession = e.study_accession INNER JOIN sample m ON e.sample_accession = m.sample_accession  where (e.platform = 'ILLUMINA' AND e.library_strategy = 'RNA-Seq' AND e.library_layout LIKE 'PAIRED%' AND e.library_selection in('PolyA', 'cDNA') AND e.library_source = 'TRANSCRIPTOMIC' AND  (UPPER(m.scientific_name) = 'HOMO SAPIENS'))", sep=" "))
# 
# GSE_inhuman_sra_poly_cdna <- c()
# for (i in gse_from_sra_cdna){
#   #print(noquote(i))
#   k=noquote(i)
#   #print (k)
#   l=noquote(paste("'",k,"' ",sep=""))
#   print (l)
#   GSE_inhuman_sra_poly_cdna <- c(paste(GSE_inhuman_sra_poly_cdna),l)
# }
# gl <- noquote(GSE_inhuman_sra_poly_cdna)
# GE_list_final=toString(gl)
# 
# hrs <-dbGetQuery(sra_con, paste("select s.*, e.*, m.* from study s INNER JOIN experiment e ON s.study_accession = e.study_accession INNER JOIN sample m ON e.sample_accession = m.sample_accession  where (e.platform = 'ILLUMINA' AND e.library_strategy = 'RNA-Seq' AND e.library_layout LIKE 'PAIRED%' AND e.library_selection in('PolyA', 'cDNA') AND s.study_alias in (",GE_list_final,") AND e.library_source = 'TRANSCRIPTOMIC' AND  (UPPER(m.scientific_name) = 'HOMO SAPIENS'))", sep=""))
# 
# sra_data_GSE_list = unique(sra_data_human_rnaseq$study_alias)
# GSE_inhuman_sra_poly_cdna <- c()
# for (i in sra_data_GSE_list){
#   if (!is.na(i) ){
#     #print (i)
#     v=startsWith(i,"GSE")
#     #print (v)
#     if (v==TRUE){
#       #print(noquote(i))
#       k=noquote(i)
#       #print (k)
#       l=noquote(paste("'",k,"' ",sep=""))
#       print (l)
#       GSE_inhuman_sra_poly_cdna <- c(paste(GSE_inhuman_sra_poly_cdna),l)
#     }
#   }
# }
# gl <- noquote(GSE_inhuman_sra_poly_cdna)
# GE_list_final=toString(gl)
# 
# geo_data_human_rnaseq<-dbGetQuery(congeo, paste("select a.ID, a.gse, a.title, c.gsm,  c.series_id, c.gpl, a.status, a.submission_date, a.last_update_date, a.type, c.source_name_ch1, c.organism_ch1, c.characteristics_ch1, c.molecule_ch1, c.treatment_protocol_ch1, c.extract_protocol_ch1, c.description, c.data_processing, a.contact, a.supplementary_file, a.pubmed_id, a.summary, a.contributor, a.web_link, a.overall_design from gse a INNER JOIN gse_gsm gg ON a.gse = gg.gse INNER JOIN gsm c ON gg.gsm = c.gsm WHERE a.gse IN (",GE_list_final,")", sep=""))
# 
# 
# #New SRA query with specific fields
# sra_data_human_rnaseq_onlyGSE <-dbGetQuery(sra_con, paste("select DISTINCT s.submission_accession,s.sradb_updated,s.study_accession,s.study_alias,s.study_title,s.study_type,s.center_project_name,s.study_description,s.study_url_link, s.study_entrez_link, s. study_attribute, s.related_studies,s.primary_study,e.experiment_alias, e.experiment_accession,e.title,e.study_name,e.sample_name,e.design_description,e.library_name,e.library_source,e.library_selection,e.library_layout,e.library_construction_protocol,e.read_spec,e.platform,e.instrument_model,e.platform_parameters,e.sequence_space,e.base_caller,e.quality_scorer,e.number_of_levels,e.multiplier,e.qtype,e.experiment_url_link,e.experiment_entrez_link,e.experiment_attribute,m.sample_alias, m.sample_accession,m.taxon_id,m.scientific_name,m.common_name,m.anonymized_name,m.individual_name,m.description,m.sample_url_link,m.sample_entrez_link,m.sample_attribute from study s INNER JOIN experiment e ON s.study_accession = e.study_accession INNER JOIN sample m ON e.sample_accession = m.sample_accession  where (e.platform = 'ILLUMINA' AND e.library_strategy = 'RNA-Seq' AND e.library_layout LIKE 'PAIRED%' AND e.library_selection in('PolyA', 'cDNA') AND e.library_source = 'TRANSCRIPTOMIC' AND  (UPPER(m.scientific_name) = 'HOMO SAPIENS') AND s.study_alias LIKE '%GSE%')", sep=" "))
# #New SRA query with specific fields only experiment details
# sra_data_human_rnaseq_expfields_onlyGSE <-dbGetQuery(sra_con, paste("select DISTINCT s.study_accession,s.study_alias,s.study_title,e.experiment_alias, e.experiment_accession,e.title,e.study_name,e.sample_name,e.design_description,e.library_name,e.library_source,e.library_selection,e.library_layout,e.library_construction_protocol,e.read_spec,e.platform,e.instrument_model,e.platform_parameters,e.sequence_space,e.base_caller,e.quality_scorer,e.number_of_levels,e.multiplier, m.sample_accession, m.taxon_id,m.scientific_name,m.description from study s INNER JOIN experiment e ON s.study_accession = e.study_accession INNER JOIN sample m ON e.sample_accession = m.sample_accession  where (e.platform = 'ILLUMINA' AND e.library_strategy = 'RNA-Seq' AND e.library_layout LIKE 'PAIRED%' AND e.library_source = 'TRANSCRIPTOMIC' AND  (UPPER(m.scientific_name) = 'HOMO SAPIENS') AND s.study_alias LIKE '%GSE%')", sep=" "))
# 
# 
# #test queries for polyA
# test_onlyGSE <-dbGetQuery(sra_con, paste("select DISTINCT s.study_accession,s.study_alias,s.study_title,e.experiment_alias, e.experiment_accession,e.title,e.study_name,e.sample_name,e.design_description,e.library_name,e.library_source,e.library_selection,e.library_layout,e.library_construction_protocol,e.read_spec,e.platform,e.instrument_model,e.platform_parameters,e.sequence_space,e.base_caller,e.quality_scorer,e.number_of_levels,e.multiplier, m.sample_accession, m.taxon_id,m.scientific_name,m.description from study s INNER JOIN experiment e ON s.study_accession = e.study_accession INNER JOIN sample m ON e.sample_accession = m.sample_accession  where (e.library_selection in('PolyA'))", sep=" "))
# test_onlyGSE <-dbGetQuery(sra_con, paste("select DISTINCT * from experiment e  where (e.library_selection in('PolyA'))", sep=" "))
# test_onlyGSE <-dbGetQuery(sra_con, paste("select DISTINCT s.submission_accession,s.sradb_updated,s.study_accession,s.study_alias,s.study_title,s.study_type,s.center_project_name,s.study_description,s.study_url_link, s.study_entrez_link, s. study_attribute, s.related_studies,s.primary_study,e.experiment_alias, e.experiment_accession,e.title,e.study_name,e.sample_name,e.design_description,e.library_name,e.library_source,e.library_selection,e.library_layout,e.library_construction_protocol,e.read_spec,e.platform,e.instrument_model,e.platform_parameters,e.sequence_space,e.base_caller,e.quality_scorer,e.number_of_levels,e.multiplier,e.qtype,e.experiment_url_link,e.experiment_entrez_link,e.experiment_attribute,m.sample_alias, m.sample_accession,m.taxon_id,m.scientific_name,m.common_name,m.anonymized_name,m.individual_name,m.description,m.sample_url_link,m.sample_entrez_link,m.sample_attribute from study s INNER JOIN experiment e ON s.study_accession = e.study_accession INNER JOIN sample m ON e.sample_accession = m.sample_accession  where (s.study_alias = 'GSE53031')", sep=" "))
# test <-dbGetQuery(sra_con, paste("select DISTINCT s.study_accession,s.study_alias,e.experiment_alias, e.experiment_accession, m.sample_accession, m.sample_alias, r.run_accession, r.run_alias from study s INNER JOIN experiment e ON s.study_accession = e.study_accession INNER JOIN sample m ON e.sample_accession = m.sample_accession INNER JOIN run r ON e.experiment_accession=r.experiment_accession  where (e.library_selection in('cDNA','PolyA') AND e.library_source = 'TRANSCRIPTOMIC' AND  (UPPER(m.scientific_name) = 'HOMO SAPIENS') AND s.study_alias LIKE '%GSE%')", sep=" "))


#get rows SRAID and GSEID from SRA for all human DATA; no crieteria for library selection right now (PolyA or CDNA)
sra_data_human_rnaseq_1 <-dbGetQuery(sra_con, paste("select DISTINCT s.study_accession, s.study_alias from study s INNER JOIN experiment e ON s.study_accession = e.study_accession INNER JOIN sample m ON e.sample_accession = m.sample_accession  where (e.platform = 'ILLUMINA' AND e.library_strategy = 'RNA-Seq' AND e.library_layout LIKE 'PAIRED%' AND e.library_source = 'TRANSCRIPTOMIC' AND  (UPPER(m.scientific_name) = 'HOMO SAPIENS') AND s.study_alias LIKE '%GSE%')", sep=" "))
#step 2 get all gse, gsm data from geo with GSE ids from above query
gseids <- unique(sra_data_human_rnaseq_1$study_alias)
temp <- c()
for (i in gseids){
  #print(noquote(i))
  k=noquote(i)
  #print (k)
  l=noquote(paste("'",k,"' ",sep=""))
  print (l)
  temp <- c(paste(temp),l)
}
gl <- noquote(temp)
GE_list_final=toString(gl)

#select information from GEO where GSE is in SRAdb from above query
geo_human_rnaseq_2<-dbGetQuery(congeo, paste("select a.*,c.* from gse a INNER JOIN gse_gsm gg ON a.gse = gg.gse INNER JOIN gsm c ON gg.gsm = c.gsm WHERE (NOT( (UPPER(c.molecule_ch1) NOT like '%POLYA%') AND (UPPER(c.description) NOT like '%MRNA%')  AND (UPPER(c.description) NOT like '%POLYA%') )) AND  UPPER(c.organism_ch1) LIKE 'HOMO SAPIENS' AND c.molecule_ch1 LIKE '%RNA%' AND a.gse IN (",GE_list_final,")", sep=""))

#some gse ids are not in gse_gsm table; search directly in gsm table series id
geo_human_rnaseq_2_2<-dbGetQuery(congeo, paste("select a.*, c.* from gse a INNER JOIN gsm c ON a.gse = c.series_id WHERE (NOT( (UPPER(c.molecule_ch1) NOT like '%POLYA%') AND (UPPER(c.description) NOT like '%MRNA%')  AND (UPPER(c.description) NOT like '%POLYA%') )) AND UPPER(c.organism_ch1) LIKE 'HOMO SAPIENS' AND c.molecule_ch1 LIKE '%RNA%' AND a.gse IN (",GE_list_final,")", sep=""))

#comibne above two res into one removing redundant rows
geo_human_rnaseq <- unique(rbind(geo_human_rnaseq_2,geo_human_rnaseq_2_2))

#store all the gsms
gsm_list<-unique(geo_human_rnaseq$gsm)

#given list of GSM make query from SRA
#query looks like select * from run where run_alias  LIKE 'GSM618465%' OR run_alias  LIKE 'GSM618466%'"
#sample_gsm_list <- c('GSM618465','GSM618466','GSM618479','GSM618489')
#gsm_list<-read.table('GSM_list_fromGEO')

# temp_gsm <- c()
# for (i in gsm_list){
#   k=noquote(i)
#   #print (k)
#   l=noquote(paste("run_alias LIKE '",k,"%'",sep=""))
#   #print (l)
#   temp_gsm <- c(paste(temp_gsm),l)
# }
# gsm <- noquote(temp_gsm)
# gsm=toString(gsm)
# gsm<-noquote(gsub(","," OR",gsm))
# #finally quey
# temp <-dbGetQuery(sra_con, paste("select * from run where ",gsm, sep=""))

#above method wont work too many conditions method #2
temp_gsm <- c()
for (i in gsm_list){
  k=noquote(i)
  #print (k)
  l=noquote(paste("'",k,"'",sep=""))
  #print (l)
  temp_gsm <- c(paste(temp_gsm),l)
}
gsm <- noquote(temp_gsm)
gsm=toString(gsm)
#finally query
#temp <-dbGetQuery(sra_con, paste("select * from run where substr(run_alias,1,9) IN (",gsm,")", sep=""))
#slower method
#temp3 <-dbGetQuery(sra_con, paste("select * from run where substr(run_alias,1,9) IN (",gsm,") OR substr(run_alias,1,10) IN (",gsm,") OR substr(run_alias,1,11) IN (",gsm,") OR substr(run_alias,1,12) IN (",gsm,") OR substr(run_alias,1,13) IN (",gsm,") OR substr(run_alias,1,14) IN (",gsm,") OR substr(run_alias,1,8) IN (",gsm,")", sep=""))

#get final SRA DATA; 
#sra_data_human_rnaseq_3 <-dbGetQuery(sra_con, paste("select DISTINCT s.study_accession, s.study_alias from study s INNER JOIN experiment e ON s.study_accession = e.study_accession INNER JOIN sample m ON e.sample_accession = m.sample_accession  INNER JOIN run r ON e.experiment_accession = r.experiment_accession where (e.platform = 'ILLUMINA' AND e.library_strategy = 'RNA-Seq' AND e.library_layout LIKE 'PAIRED%' AND e.library_source = 'TRANSCRIPTOMIC' AND  (UPPER(m.scientific_name) = 'HOMO SAPIENS') AND s.study_alias LIKE '%GSE%' AND r.run_alias IN (select run_alias from run where substr(run_alias,1,9) IN (",gsm,") OR substr(run_alias,1,10) IN (",gsm,") OR substr(run_alias,1,11) IN (",gsm,") OR substr(run_alias,1,12) IN (",gsm,") OR substr(run_alias,1,13) IN (",gsm,") OR substr(run_alias,1,14) IN (",gsm,") OR substr(run_alias,1,8) IN (",gsm,")))", sep=" "))

##finally extract all SRA data with run alias in the gsm list
sra_data_human_rnaseq_3 <-dbGetQuery(sra_con, paste("select DISTINCT s.*, e.*,  m.*, r.* from study s INNER JOIN experiment e ON s.study_accession = e.study_accession INNER JOIN sample m ON e.sample_accession = m.sample_accession  INNER JOIN run r ON e.experiment_accession = r.experiment_accession where (e.platform = 'ILLUMINA' AND e.library_strategy = 'RNA-Seq' AND e.library_layout LIKE 'PAIRED%' AND e.library_source = 'TRANSCRIPTOMIC' AND  (UPPER(m.scientific_name) = 'HOMO SAPIENS') AND s.study_alias LIKE '%GSE%' AND r.run_alias IN (select run_alias from run where substr(run_alias,1,9) IN (",gsm,") OR substr(run_alias,1,10) IN (",gsm,") OR substr(run_alias,1,11) IN (",gsm,") OR substr(run_alias,1,12) IN (",gsm,") OR substr(run_alias,1,13) IN (",gsm,") OR substr(run_alias,1,14) IN (",gsm,") OR substr(run_alias,1,8) IN (",gsm,")))", sep=" "))

##join SRA and GEO data by GSM ids
##first split run_alias in SRA by '_' and add GSM part as a column to SRA_data_human
all_run_alias<-sra_data_human_rnaseq_3$run_alias
library(stringr)
run_alias_2<-c()
for (i in all_run_alias){
  
  run_alias_2<-c(paste(run_alias_2),str_split_fixed(i,"_",2)[1])
  
}
SRA_data_human$run_alias_2 <- run_alias_2

SRA_GEO_merged_data_human <- merge(SRA_data_human,geo_human_rnaseq,by.x = c("run_alias_2"),by.y = c("gsm"))

##write required cols to file
write.csv(SRA_GEO_merged_data_human[,c("submission_accession","sradb_updated","study_accession","study_type","study_type","center_project_name","study_description","study_url_link","study_entrez_link","study_attribute","related_studies","primary_study","experiment_alias","experiment_accession","title.x","study_name","sample_name","design_description","library_name","library_source","library_selection","library_layout","library_construction_protocol","read_spec","platform","instrument_model","instrument_name","platform_parameters","sequence_space","base_caller","quality_scorer","multiplier","qtype","experiment_url_link","experiment_entrez_link","experiment_attribute","sample_alias","sample_accession","taxon_id","common_name","anonymized_name","individual_name","description.x","sample_url_link","sample_entrez_link","sample_attribute","run_alias","run_accession","run_date","sradb_updated.3","spot_length","base_caller","run_center","experiment_name","run_url_link","run_entrez_link","run_attribute","study_abstract","gse","title.y","run_alias_2","series_id","gpl","status","submission_date","last_update_date","type","source_name_ch1","organism_ch1","molecule_ch1","label_ch1","treatment_protocol_ch1","extract_protocol_ch1","description.y","data_processing","contact","supplementary_file")],"s3.csv")
