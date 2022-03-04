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


#### ----- Load Single Cell Data ------ ####
# pre-processed Single Cell data, GSE75688
# Genes not expressed in a least 20% of the dataset were removed.
# afterwards, only samples from tumor were kept

data("GSE75688_TPM_tumor", package = "DynamicCancerDriver")

#----- Find Dynamic Cancer Drivers, PPI top 40% -----
DCD.HER2time_SC <- findDCD(GeneExpression = GSE75688_TPM_tumor
                           , pathCovariate = "HER2"
                           , PPItop = 0.4
                           , findEvent = TRUE
                           , alpha = 0.05)


write.csv(DCD.HER2time_SC$res$CDinfer
         , file =  "supplementary table 1 - dynamic cancer drivers HER2time(SC).csv")

DCD.VIMtime_SC <- findDCD(GeneExpression = GSE75688_TPM_tumor
                           , pathCovariate = "VIM"
                           , PPItop = 0.4
                           , findEvent = TRUE)

write.csv(DCD.VIMtime_SC$res$CDinfer
          , file =  "supplementary table 2 - dynamic cancer drivers VIMtime(SC).csv")

# "Combined" as the union of inferred DCD
# keeping the highest relative Causal Impact
unionSet <- union(DCD.HER2time_SC$res$CDinfer$Ensembl.ID
                          ,DCD.VIMtime_SC$res$CDinfer$Ensembl.ID)
Combined <-
  rbind(DCD.HER2time_SC$res$CDinfer,DCD.VIMtime_SC$res$CDinfer)%>%
  filter(Ensembl.ID %in% unionSet) %>%
  arrange(desc(abs(RelEffect)))%>%
  distinct(Ensembl.ID, .keep_all = T)
