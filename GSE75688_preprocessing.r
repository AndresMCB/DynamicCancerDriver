#################################################################
# Script to pre-process single cell data GSE75688
# Download final_sample_information and raw_TPM_matrix from:
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE75688
#

if(!require(tidyverse))
  install.packages("tidyverse")

library(tidyverse)
library(AMCBGeneUtils)

#---- Load TPM data (it is assumed file is in the work directory) ----

GSE75688_TPM <-
  read_delim("./GSE75688_GEO_processed_Breast_Cancer_raw_TPM_matrix.txt"
             , "\t"
             , escape_double = FALSE
             , trim_ws = TRUE)
GSE75688_TPM <- as.data.frame(GSE75688_TPM)
gene_id <- GSE75688_TPM$gene_id

# remove "gene_id" "gene_name"  and  "gene_type"
GSE75688_TPM <- data.matrix(GSE75688_TPM[,-c(1:3), drop = F])
row.names(GSE75688_TPM) <- gene_id
GSE75688_TPM <- t(GSE75688_TPM)

#---- Load sample information (it is assumed file is in the work directory) ----
GSE75688_sample_information <-
  read_delim("./GSE75688_final_sample_information.txt"
             , "\t"
             , escape_double = FALSE
             , trim_ws = TRUE)

#---- remove genes not expressed in at least 20% of the samples ----
GSE75688_TPM_tumor <-
  GSE75688_TPM[,colSums(GSE75688_TPM>0)>(0.2*nrow(GSE75688_TPM))
                     , drop=F]

#---- Keep only samples that are single cell from tumors ----
TumorSamples <- GSE75688_sample_information%>%
  filter(type =="SC",index=="Tumor")
index <- row.names(GSE75688_TPM)%in%TumorSamples$sample
GSE75688_TPM_tumor <- GSE75688_TPM_tumor[index
                                     ,,drop=FALSE]

#---- remove version from Ensembl.ID ----
temp <- colnames(GSE75688_TPM_tumor)
temp <-  gsub("\\..*","",temp)
colnames(GSE75688_TPM_tumor)<- temp

#---- remove any column with no equivalence un HGNC.symbol ----
IdSource <- GeneIdSource(colnames(GSE75688_TPM_tumor))
genes <- changeGeneId(colnames(GSE75688_TPM_tumor),from = IdSource)
GSE75688_TPM_tumor <- GSE75688_TPM_tumor[,!is.na(genes$HGNC.symbol),drop=F]

save(GSE75688_TPM_tumor, file = "GSE75688_TPM_tumor.rda")



