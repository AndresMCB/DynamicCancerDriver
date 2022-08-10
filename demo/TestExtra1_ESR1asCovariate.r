#####---Additional experiment 1: ESR1 (a BRCA Driver) as path covariate ---#####

library(tidyverse)
library(DynamicCancerDriver)

#####---------  Loading Datasets (TCGA-BRCA andGSE75688_TPM_tumor) ---------#####
# Dataset downloaded and normalised using TCGAbiolinks (March 2021).
# Only samples from primary tumor were downloaded
# Included in DynamicCancerDriver package as TCGA_BRCA_TP_NormCounts.rda
# ----- Single Cell Data ------
# pre-processed Single Cell data, GSE75688
# Genes not expressed in a least 20% of the dataset were removed.
# afterwards, only samples from tumor were kept

data("GSE75688_TPM_tumor", package = "DynamicCancerDriver")
data("TCGA_BRCA_TP_NormCounts", package = "DynamicCancerDriver")
# remove genes with no Hugo Symbol
genes <- AMCBGeneUtils::changeGeneId(colnames(TCGA_BRCA_TP_NormCounts)
                                     , from = "Ensembl.ID")
TCGA_BRCA <- TCGA_BRCA_TP_NormCounts[,!is.na(genes$HGNC.symbol),drop=F]

# remove genes not expressed in at least 20% of the samples
TCGA_BRCA <-TCGA_BRCA[,colSums(TCGA_BRCA>0)>(0.2*nrow(TCGA_BRCA))
                      , drop=F]

# For comparison purposes, we keep only samples from the same patients
# analysed in CBNA paper

data("CBNApaper.patients", package = "DynamicCancerDriver")
patients <- str_sub(row.names(TCGA_BRCA), start = 1, end = 12)
index <- patients%in%CBNApaper.patients
TCGA_BRCA <- TCGA_BRCA[index,, drop=FALSE]


#  Find Dynamic Cancer Drivers (Single cell), PPI top 40%
DCD.ESR1time_SC <- findDCD(GeneExpression = GSE75688_TPM_tumor
                           , pathCovariate = "ESR1"
                           , PPItop = 0.4
                           , findEvent = TRUE
                           , alpha = 0.05)

DCD.ESR1time_SC$res$summary
data(CGC.driverNames,package = "DynamicCancerDriver")
intersect(DCD.ESR1time_SC$res$CDinfer$Ensembl.ID
          ,CGC.driverNames$Ensembl.ID)
length(intersect(DCD.ESR1time_SC$res$CDinfer$Ensembl.ID
                 ,CGC.driverNames$Ensembl.ID))/
  length(DCD.ESR1time_SC$res$CDinfer$Ensembl.ID)
write.csv(DCD.ESR1time_SC$res$CDinfer
          ,file = "Supplementary table 12 - dynamic cancer drivers ESR1time(SC).csv")

#  Find Dynamic Cancer Drivers (BULK data), PPI top 40%
DCD.ESR1time_Bulk <- findDCD(GeneExpression = TCGA_BRCA
                             , pathCovariate = "ESR1"
                             , PPItop = 0.4
                             , findEvent = TRUE
                             , project = "BRCA")

top <- c("50", "100","150", "200")
ESR1BulkPerfomance <- matrix(nrow =1, ncol = 4)
colnames(ESR1BulkPerfomance) <- top
aux <- DCD.ESR1time_Bulk$res$CDinfer
for (i in top) {
  index <- as.numeric(1:i)
  ESR1BulkPerfomance[1,i] <-
    length(intersect(CGC.driverNames$Ensembl.ID
                     ,aux$Ensembl.ID[index]))
}
ESR1BulkPerfomance
write.csv(DCD.ESR1time_Bulk$res$CDinfer
          ,file = "supplementary table 13 - dynamic cancer drivers ESR1time(Bulk).csv")
