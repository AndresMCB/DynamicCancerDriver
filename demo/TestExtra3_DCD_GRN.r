#####---Additional experiment 3: Regulatory relationships in DCD top 100 ---#####
#* regulatory relationships between the top 200 of DCDs
#* and the set (not restricted to the top 200 DCDs) of
#* CGC Breast Cancer drivers detected by our method
#* (when using “HER2time(SC)” and “HER2time(Bulk)”
#* as the path covariates respectively).
#

library(tidyverse)

# loading highly confident grn  (BRCA_TCGA from http://www.grndb.com/)

aux <- system.file("extdata/"
                   , "BRCA_TCGA-regulons.txt"
                   , package = "DynamicCancerDriver")
BRCA_TCGA.regulons <- read.delim(aux)

# keep only confidence == "High"
BRCA_TCGA.regulons <- BRCA_TCGA.regulons%>%
  dplyr::filter(Confidence == "High")%>%
  dplyr::select(TF, gene, NES, Confidence)


# loading the discovered DCD from our original experiments
aux <- system.file("Supplementary/"
                           , "supplementary table 1 - dynamic cancer drivers HER2time(SC).csv"
                           , package = "DynamicCancerDriver")

DCD.HER2_SC <- read.csv(aux)


aux <- system.file("Supplementary/"
                   , "supplementary table 8 - dynamic cancer drivers HER2time(Bulk).csv"
                   , package = "DynamicCancerDriver")
DCD.HER2_bulk <- read.csv(aux)
rm(aux)
data(CGC.Breast)

#bulk
aux <- intersect(DCD.HER2_bulk$HGNC.symbol, CGC.Breast$HGNC.symbol)
aux <- union (DCD.HER2_bulk$HGNC.symbol[1:200], aux)
regulons.HER2_bulk <- BRCA_TCGA.regulons%>%
  dplyr::filter(TF %in% aux, gene %in% aux)%>%
  mutate(TF.isCGC = TF %in% CGC.Breast$HGNC.symbol, .before = 2)%>%
  mutate(gene.isCGC = gene %in% CGC.Breast$HGNC.symbol, .before = 4)

# Single Cell
aux <- intersect (DCD.HER2_SC$HGNC.symbol
                  , CGC.Breast$HGNC.symbol)
aux <- union (DCD.HER2_SC$HGNC.symbol[1:200]
              , aux)

regulons.HER2_SC <- BRCA_TCGA.regulons%>%
  dplyr::filter(TF %in% aux, gene %in% aux)%>%
  mutate(TF.isCGC = TF %in% CGC.Breast$HGNC.symbol, .before = 2)%>%
  mutate(gene.isCGC = gene %in% CGC.Breast$HGNC.symbol, .before = 4)

write.csv(regulons.HER2_SC
          , file = "Supplementary table 10 - DCD(HER2_SC) Reg. relationships.csv")
write.csv(regulons.HER2_bulk
          , file = "Supplementary table 11 - DCD(HER2_Bulk) Reg. relationships.csv")
