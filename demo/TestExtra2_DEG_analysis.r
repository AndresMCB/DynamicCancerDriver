
#####---Additional experiment 2: DEG analysis (normal-cancer) ---#####
#
# Please download the full TCGA-BRCA data from
# http://4llab.net/Bioinformatics/TCGA.BRCA.rda
# Load the dataset using file/open file...
# you can use the following code (the dataset needs to be in your wd
#           wdir <- getwd()
#           TCGA.BRCA <- load(paste0(wdir,"/TCGA.BRCA.rda"))

# Diff.expr.analysis (DEA)
#  NOTE: TCGABiolinks requires that each row represents a gene,
#  and each column represents a sample with Cond_type

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

library(DynamicCancerDriver)
library(tidyverse)
library(TCGAbiolinks)


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

write.csv(dataDEGs
          ,file = "supplementary table 14 - TCGA-BRCA DEGs.csv")
# loading the discovered DCD from our original experiments
aux <- system.file("Supplementary/"
                   , "supplementary table 8 - dynamic cancer drivers HER2time(Bulk).csv"
                   , package = "DynamicCancerDriver")
DCD.HER2_bulk <- read.csv(aux)

aux <- system.file("Supplementary/"
                   , "/supplementary table 9 - dynamic cancer drivers VIMtime(Bulk).csv"
                   , package = "DynamicCancerDriver")
DCD.VIM_bulk <- read.csv(aux)

# Analysing top 100 inferred DCD
data(CGC.Breast)
DCD.HER2noDEG <- setdiff(DCD.HER2_bulk$HGNC.symbol[1:100]
                         ,dataDEGs$gene_name)

DCD.HER2noDEG <- AMCBGeneUtils::changeGeneId(DCD.HER2noDEG)

DCD.HER2noDEG <- DCD.HER2noDEG%>%
  transmute(Ensembl.ID, HGNC.symbol
            , is.CGC = HGNC.symbol%in%CGC.driverNames$HGNC.symbol
            , is.BRCA_CGC =HGNC.symbol%in%CGC.Breast$HGNC.symbol)

DCD.VIMnoDEG <- setdiff(DCD.VIM_bulk$HGNC.symbol[1:100]
                        ,dataDEGs$gene_name)

DCD.VIMnoDEG <- AMCBGeneUtils::changeGeneId(DCD.VIMnoDEG)

DCD.VIMnoDEG <- DCD.VIMnoDEG%>%
  transmute(Ensembl.ID, HGNC.symbol
            , is.CGC = HGNC.symbol%in%CGC.driverNames$HGNC.symbol
            , is.BRCA_CGC =HGNC.symbol%in%CGC.Breast$HGNC.symbol)

write.csv(DCD.HER2noDEG
          ,file = "supplementary table 15 - DCD.HER2noDEG.csv")
write.csv(DCD.VIMnoDEG
          ,file = "supplementary table 16 - DCD.VIMnoDEG.csv")


sum(DCD.HER2noDEG$is.CGC)
sum(DCD.HER2noDEG$is.BRCA_CGC)
sum(DCD.VIMnoDEG$is.CGC)
sum(DCD.VIMnoDEG$is.BRCA_CGC)
