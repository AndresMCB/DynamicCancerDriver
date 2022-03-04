#### ----- Script to test DynamicCancerDriver package ------ ####
#
#  This script follows the procedure described in the
#  Bioinformatic-Oxford Paper.
#
#

#rm(list = ls())
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
HER2_SC_4_02 <- findDCD(GeneExpression = GSE75688_TPM_tumor
                           , pathCovariate = "HER2"
                           , PPItop = 0.4
                           , findEvent = F
                           , alpha = 0.02)

#--- trying to reduce false positive rate

putativeDCD <- HER2_SC_4_02$res$CDinfer$Ensembl.ID
GeneExpression = as.tibble(GSE75688_TPM_tumor)
temp <- vector(length = length(putativeDCD))
names(temp) <- putativeDCD
for (i in putativeDCD[4]) {
  aux <-DCD.findEvent(z = HER2_SC_4_02$z
                       , GeneExpression = GeneExpression
                       , pathCovariate.name = i
                         , Step=10)
  temp[i] <- aux$eventAt
}





