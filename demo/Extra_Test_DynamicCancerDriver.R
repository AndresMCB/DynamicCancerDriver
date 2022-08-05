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

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!require("TCGAbiolinks", quietly = TRUE))
  BiocManager::install("TCGAbiolinks")

library(monocle3)
library(DynamicCancerDriver)
library(tidyverse)
library(TCGAbiolinks)


# Function for computing Jaccard Similarity
jaccard_similarity <- function(A, B) {
  intersection = length(intersect(A, B))
  union = length(A) + length(B) - intersection
  return (intersection/union)
}


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


#####---Additional experiment 1: BRCA Driver as path covariate ---#####

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
          ,file = "supplementary table 10 - dynamic cancer drivers ESR1time(SC).csv")

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
          ,file = "supplementary table 11 - dynamic cancer drivers ESR1time(Bulk).csv")

#####---Additional experiment 2: DEG analysis (normal-cancer) ---#####
# Please download the full TCGA-BRCA data from
# http://4llab.net/Bioinformatics/TCGA.BRCA.rda

# Load the dataset using file/open file...
# you can use the following code (the dataset needs to be in your wd
#           wdir <- getwd()
#           TCGA.BRCA <- load(paste0(wdir,"/TCGA.BRCA.rda"))

# Diff.expr.analysis (DEA)
#  NOTE: TCGABiolinks requires that each row represents a gene,
#  and each column represents a sample with Cond_type

GE <- TCGA.BRCA$HTSeq_Norm_Counts
# remove genes not expressed in at least 20% of the samples
GE <-GE[,colSums(GE>0)>(0.2*nrow(GE))
                      , drop=F]
index <- TCGA.BRCA$clinical%>%
  dplyr::filter(shortLetterCode %in% "NT")%>%
  dplyr::select(barcode)
normal <- GE[index$barcode,,drop = F]

index <- TCGA.BRCA$clinical%>%
  dplyr::filter(shortLetterCode %in% "TP")%>%
  dplyr::select(barcode)
tumour <- GE[index$barcode,,drop = F]

dataDEGs <- TCGAanalyze_DEA(mat1 = t(normal),
                            mat2 = t(tumour),
                            Cond1type = "Normal",
                            Cond2type = "Primary",
                            fdr.cut = 0.01 ,
                            logFC.cut = 1,
                            method = "glmLRT")

# loading the discovered DCD from our original experiments
 wdir <- getwd()
 DCD.HER2_bulk <-
   read.csv(paste0(wdir,"/supplementary table 8 - dynamic cancer drivers HER2time(Bulk).csv"))
 DCD.VIM_bulk <-
   read.csv(paste0(wdir,"/supplementary table 9 - dynamic cancer drivers VIMtime(Bulk).csv"))

# Analysing top 100 inferred DCD
DCD.HER2noDEG <- setdiff(DCD.HER2_bulk$HGNC.symbol[1:100]
                         ,dataDEGs$gene_name)
DCD.VIMnoDEG <- setdiff(DCD.VIM_bulk$HGNC.symbol[1:100]
                         ,dataDEGs$gene_name)

intersect(DCD.HER2noDEG, CGC.driverNames$HGNC.symbol)
intersect(DCD.VIMnoDEG, CGC.driverNames$HGNC.symbol)

data(CGC.Breast)

intersect(DCD.HER2noDEG, CGC.Breast$HGNC.symbol)
intersect(DCD.VIMnoDEG, CGC.Breast$HGNC.symbol)


#####---Additional experiment 3: Regulatory relationships in DCD top 100 ---#####
# loading highly confident grn from http://www.grndb.com/
 wdir <- getwd()
 BRCA_TCGA.regulons <- read.delim(paste0(wdir,"/BRCA_TCGA-regulons.txt"))

# keep only confidence == "High"
BRCA_TCGA.regulons <- BRCA_TCGA.regulons%>%
  dplyr::filter(Confidence == "High")%>%
  dplyr::select(TF, gene, NES, Confidence)

#bulk
aux <- intersect(DCD.HER2_bulk$HGNC.symbol, CGC.Breast$HGNC.symbol)
aux <- union (DCD.HER2_bulk$HGNC.symbol[1:200], aux)
regulons.HER2_bulk <- BRCA_TCGA.regulons%>%
  dplyr::filter(TF %in% aux, gene %in% aux)%>%
  mutate(TF.isCGC = TF %in% CGC.Breast$HGNC.symbol, .before = 2)%>%
  mutate(gene.isCGC = gene %in% CGC.Breast$HGNC.symbol, .before = 4)

# Single Cell
aux <- intersect (DCD.HER2time_SC$res$CDinfer$HGNC.symbol
              , CGC.Breast$HGNC.symbol)
aux <- union (DCD.HER2time_SC$res$CDinfer$HGNC.symbol[1:200]
                  , aux)

regulons.HER2_SC <- BRCA_TCGA.regulons%>%
  dplyr::filter(TF %in% aux, gene %in% aux)%>%
  mutate(TF.isCGC = TF %in% CGC.Breast$HGNC.symbol, .before = 2)%>%
  mutate(gene.isCGC = gene %in% CGC.Breast$HGNC.symbol, .before = 4)


write.csv(regulons.HER2_bulk, file = "regulons.HER2_bulk.csv")
write.csv(regulons.HER2_SC, file = "regulons.HER2_SC.csv")


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
                               , findEvent = T
                               , project = "BRCA")

DCD.VIMmonocle_SC <- findDCD(GeneExpression = GSE75688_TPM_tumor
                              , pathCovariate = "VIM"
                              , z = MonoclePseudotime
                              , PPItop = 0.4
                              , findEvent = T
                              , project = "BRCA")

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
