library(TCGAbiolinks)
library(dplyr)
library(DT)
library(data.table)
library(plyr)
packageVersion("TCGAbiolinks")


query <- GDCquery(project = "TCGA-CHOL",  data.category = "Clinical", file.type = "xml")
GDCdownload(query)
clinical <- GDCprepare_clinic(query, clinical.info = "patient")

queryB <- GDCquery(project = "TCGA-CHOL",  data.category = "Biospecimen", file.type = "xml")
GDCdownload(queryB)

aliquot <- GDCprepare_clinic(queryB, clinical.info = c("aliquot"))
sample <- GDCprepare_clinic(queryB, clinical.info = c("sample"))
bio_patient <- GDCprepare_clinic(queryB, clinical.info = c("bio_patient"))
analyte <- GDCprepare_clinic(queryB, clinical.info = c("analyte"))
portion <- GDCprepare_clinic(queryB, clinical.info = c("portion"))
protocol <- GDCprepare_clinic(queryB, clinical.info = c("protocol"))
slide <- GDCprepare_clinic(queryB, clinical.info = c("slide"))

hbp<-as.data.frame(names(bio_patient))
hanalyte<-as.data.frame(names(analyte))
hportion<-as.data.frame(names(portion))
hprot<-as.data.frame(names(protocol))
hslid<-as.data.frame(names(slide))


join(aliquot,bio_patient)

#remove duplicated rows from all tables
