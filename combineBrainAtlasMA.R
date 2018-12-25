library(dplyr)
library(readr)





MicroarrayExpression <- read_csv("~/Downloads/humanBrainAtlas/microarray/normalized_microarray_donor9861/MicroarrayExpression.csv", col_names = FALSE)

Probes <- read_csv("~/Downloads/humanBrainAtlas/microarray/normalized_microarray_donor9861/Probes.csv")

SampleAnnot <- read_csv("~/Downloads/humanBrainAtlas/microarray/normalized_microarray_donor9861/SampleAnnot.csv")
#well id is uniq
colnames(MicroarrayExpression)[1]<-"probe_id"
colnames(MicroarrayExpression)[2:ncol(MicroarrayExpression)]<-SampleAnnot$well_id

MicroarrayExpression_joined<-inner_join(Probes,MicroarrayExpression)

sum(duplicated(colnames(MicroarrayExpression_joined)))
names(MicroarrayExpression_joined)<-colnames(MicroarrayExpression_joined)
write_tsv(MicroarrayExpression_joined,"JoinedData.tsv",col_names = TRUE)


