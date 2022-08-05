#####---Additional experiment: Regulatory relationships in DCD top 100 ---#####

library(tidyverse)

# loading highly confident grn (from http://www.grndb.com/)

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
aux <- intersect (DCD.HER2time_SC$res$CDinfer$HGNC.symbol
                  , CGC.Breast$HGNC.symbol)
aux <- union (DCD.HER2time_SC$res$CDinfer$HGNC.symbol[1:200]
              , aux)

regulons.HER2_SC <- BRCA_TCGA.regulons%>%
  dplyr::filter(TF %in% aux, gene %in% aux)%>%
  mutate(TF.isCGC = TF %in% CGC.Breast$HGNC.symbol, .before = 2)%>%
  mutate(gene.isCGC = gene %in% CGC.Breast$HGNC.symbol, .before = 4)


write.csv(regulons.HER2_bulk
          , file = "Supplementary table 10 –regulons.HER2_bulk.csv")
write.csv(regulons.HER2_SC
          , file = "Supplementary table 11regulons.HER2_SC.csv")
