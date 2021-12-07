#### ----- Script to test DynamicCancerDriver package ------ ####
#
#  This script follows the procedure described in the
#  Bioinformatic-Oxford Paper.
#
#

rm(list = ls())
if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")

if (!requireNamespace("tidyverse", quietly = TRUE))
  install.packages("tidyverse")

if (!requireNamespace("AMCBGeneUtils", quietly = TRUE))
  devtools::install_github(repo = "AndresMCB/AMCBGeneUtils")

if (!requireNamespace("DynamicCancerDriver", quietly = TRUE))
  devtools::install_github(repo = "AndresMCB/DynamicCancerDriver")

library(DynamicCancerDriver)
library(tidyverse)

####---- Load TCGA-BRCA Data -----####
# Dataset downloaded and normalised using TCGAbiolinks (March 2021).
# Only samples from primary tumor were downloaded
# Included in DynamicCancerDriver package as TCGA_BRCA_TP_NormCounts.rda

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

#----- Find Dynamic Cancer Drivers (BULK data), PPI top 40% -----


DCD.HER2time_Bulk <- findDCD(GeneExpression = TCGA_BRCA
                           , pathCovariate = "HER2"
                           , PPItop = 0.4
                           , findEvent = TRUE
                           , project = "BRCA")
write.csv(DCD.VIMtime_Bulk$res$CDinfer
          , file =  "supplementary table 8 - dynamic cancer drivers HER2time(Bulk).csv")


DCD.VIMtime_Bulk <- findDCD(GeneExpression = TCGA_BRCA
                          , pathCovariate = "VIM"
                          , PPItop = 0.4
                          , findEvent = TRUE
                          , project = "BRCA")
write.csv(DCD.VIMtime_Bulk$res$CDinfer
          , file =  "supplementary table 9 - dynamic cancer drivers VIMtime(Bulk).csv")






