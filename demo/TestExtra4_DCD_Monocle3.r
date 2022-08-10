#### ----- Script with additional experiments suggested by reviewers ------ ####
#
#  This script follows the procedure described in the
#  Briefings in Functional Genomics - Oxford Paper.
#

#####---------  Loading required packages ---------#####


if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")

if (!requireNamespace("tidyverse", quietly = TRUE))
  install.packages("tidyverse")

if (!requireNamespace("AMCBGeneUtils", quietly = TRUE))
  devtools::install_github(repo = "AndresMCB/AMCBGeneUtils")

if (!requireNamespace("DynamicCancerDriver", quietly = TRUE))
  devtools::install_github(repo = "AndresMCB/DynamicCancerDriver")

if (!require("monocle3", quietly = TRUE))
  BiocManager::install("monocle3")


library(monocle3)
library(DynamicCancerDriver)
library(tidyverse)



# Function for computing Jaccard Similarity
jaccard_similarity <- function(A, B) {
  intersection = length(intersect(A, B))
  union = length(A) + length(B) - intersection
  return (intersection/union)
}


# ----- Single Cell Data ------
# pre-processed Single Cell data, GSE75688
# Genes not expressed in a least 20% of the dataset were removed.
# afterwards, only samples from tumor were kept

data("GSE75688_TPM_tumor", package = "DynamicCancerDriver")







#####---Additional experiment 4: Performance for a provided pseudotime ---#####
# Monocle 3

cds <- new_cell_data_set(t(GSE75688_TPM_tumor))
## Step 1: Normalize and pre-process the data
cds <- preprocess_cds(cds, num_dim = 100)
## Step 3: Reduce the dimensions using UMAP
cds <- reduce_dimension(cds, reduction_method = "UMAP")
## Step 4: Cluster the cells
# Note: for consistency in results, please choose
# the rightmost node as zero time
cds <- cluster_cells(cds)
## Step 5: Learn a graph
cds <- learn_graph(cds)
## Step 6: Order cells
cds <- order_cells(cds)

plot_cells(cds)
MonoclePseudotime <- cds@principal_graph_aux@listData[["UMAP"]][["pseudotime"]]
sum(duplicated(MonoclePseudotime))

DCD.HER2monocle_SC <- findDCD(GeneExpression = GSE75688_TPM_tumor
                              , pathCovariate = "HER2"
                              , z =MonoclePseudotime
                              , PPItop = 0.4
                              , findEvent = T)
#, project = "BRCA")
write.csv(DCD.HER2monocle_SC$res$CDinfer
          , file = "Supplementary table 17 - DCDs(HER2)monocle3_SC.csv")

DCD.VIMmonocle_SC <- findDCD(GeneExpression = GSE75688_TPM_tumor
                             , pathCovariate = "VIM"
                             , z = MonoclePseudotime
                             , PPItop = 0.4
                             , findEvent = T)
#, project = "BRCA")
write.csv(DCD.VIMmonocle_SC$res$CDinfer
          , file = "Supplementary table 18 - DCDs(VIM)monocle3_SC.csv")

jaccard_similarity(DCD.VIMmonocle_SC$res$CDinfer$Ensembl.ID
                   ,DCD.HER2monocle_SC$res$CDinfer$Ensembl.ID)

aux1 <- intersect(DCD.HER2monocle_SC$res$CDinfer$Ensembl.ID
                  , CGC.driverNames$Ensembl.ID)
aux2 <- intersect(DCD.VIMmonocle_SC$res$CDinfer$Ensembl.ID
                  , CGC.driverNames$Ensembl.ID)
jaccard_similarity(aux1,aux2)

# PhenoPath
DCD.HER2time_SC <- findDCD(GeneExpression = GSE75688_TPM_tumor
                           , pathCovariate = "HER2"
                           , PPItop = 0.4
                           , findEvent = TRUE
                           , alpha = 0.05)

DCD.VIMtime_SC <- findDCD(GeneExpression = GSE75688_TPM_tumor
                          , pathCovariate = "VIM"
                          , PPItop = 0.4
                          , findEvent = TRUE)


jaccard_similarity(DCD.HER2time_SC$res$CDinfer$Ensembl.ID
                   ,DCD.VIMtime_SC$res$CDinfer$Ensembl.ID)

aux1 <- intersect(DCD.HER2time_SC$res$CDinfer$Ensembl.ID
                  , CGC.driverNames$Ensembl.ID)
aux2 <- intersect(DCD.VIMtime_SC$res$CDinfer$Ensembl.ID
                  , CGC.driverNames$Ensembl.ID)
jaccard_similarity(aux1,aux2)
