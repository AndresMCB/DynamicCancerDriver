#### ----- Script with additional experiments suggested by reviewers ------ ####
#
#  This script follows the procedure described in the
#  Briefings in Functional Genomics - Oxford Paper.
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

if (!requireNamespace("SCORPIUS", quietly = TRUE))
  install.packages("SCORPIUS")

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

# #----- Find Dynamic Cancer Drivers (BULK data), PPI top 40% -----
#
#
# DCD.ESR1time_Bulk <- findDCD(GeneExpression = TCGA_BRCA
#                              , pathCovariate = "ESR1"
#                              , PPItop = 0.4
#                              , findEvent = TRUE
#                              , project = "BRCA")
#
# DCD.AKT1time_Bulk <- findDCD(GeneExpression = TCGA_BRCA
#                                , pathCovariate = "AKT1"
#                                , PPItop = 0.4
#                                , findEvent = TRUE
#                                , project = "BRCA")
#
#
# write.csv(DCD.ESR1time_Bulk$res$CDinfer
#           , file =  "supplementary table 10 - dynamic cancer drivers ESR1time(Bulk).csv")
#

#### ----- Load Single Cell Data ------ ####
# pre-processed Single Cell data, GSE75688
# Genes not expressed in a least 20% of the dataset were removed.
# afterwards, only samples from tumor were kept

data("GSE75688_TPM_tumor", package = "DynamicCancerDriver")

#----- Find Dynamic Cancer Drivers, PPI top 40% -----
DCD.ESR1time_SC <- findDCD(GeneExpression = GSE75688_TPM_tumor
                           , pathCovariate = "ESR1"
                           , PPItop = 0.4
                           , findEvent = TRUE
                           , alpha = 0.05)

DCD.AKT1time_SC <- findDCD(GeneExpression = GSE75688_TPM_tumor
                           , pathCovariate = "AKT1"
                           , PPItop = 0.4
                           , findEvent = TRUE
                           , alpha = 0.05)


#------ dynamic driver genes that are not differentially
#------ expressed (cancer vs normal samples)
# Note: the SC data set does not have normal samples



#----- Test stability to a different pseudotime ----
library(SCORPIUS)
library(monocle3)

# Single cell
space <- reduce_dimensionality(GSE75688_TPM_tumor, "spearman", ndim = 317)
set.seed(1)
traj <- SCORPIUS::infer_trajectory(space)

 draw_trajectory_plot(
   space,
   #progression_group = group_name,
   path = traj$path,
   contour = TRUE
 )

DCD.HER2scorpius_SC <- findDCD(GeneExpression = GSE75688_TPM_tumor
                             , pathCovariate = "HER2"
                             , z = scale(traj$time)
                             , PPItop = 0.4
                             , findEvent = T
                             , project = "BRCA")

DCD.VIMscorpius_SC <- findDCD(GeneExpression = GSE75688_TPM_tumor
                               , pathCovariate = "VIM"
                               , z = scale(traj$time)
                               , PPItop = 0.4
                               , findEvent = T
                               , project = "BRCA")


# Monocle 3
cds <- new_cell_data_set(t(GSE75688_TPM_tumor))
set.seed(1)
## Step 1: Normalize and pre-process the data
cds <- preprocess_cds(cds, num_dim = 300)
## Step 3: Reduce the dimensions using UMAP
cds <- reduce_dimension(cds)
## Step 4: Cluster the cells
cds <- cluster_cells(cds)
## Step 5: Learn a graph
cds <- learn_graph(cds)
## Step 6: Order cells
cds <- order_cells(cds)

plot_cells(cds)
MonoclePseudotime <- cds@principal_graph_aux@listData[["UMAP"]][["pseudotime"]]


DCD.HER2monocle_SC <- findDCD(GeneExpression = GSE75688_TPM_tumor
                               , pathCovariate = "HER2"
                               , z = scale(MonoclePseudotime)
                               , PPItop = 0.4
                               , findEvent = T
                               , project = "BRCA")

DCD.VIMmonocle_SC <- findDCD(GeneExpression = GSE75688_TPM_tumor
                              , pathCovariate = "VIM"
                              , z = scale(MonoclePseudotime)
                              , PPItop = 0.4
                              , findEvent = T
                              , project = "BRCA")
