#rm(list = ls())
library(R.utils)
library(CausalImpact)
library(phenopath)
library(parallel)
library(AMCBGeneUtils)
library(tidyverse)
library(readr)

#*************************************************************************
#*            Load Bulk Data

defaultDir <- getwd()
dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(dir)
sourceDirectory("./functions", modifiedOnly = F)
loadData()
setwd(defaultDir)

primaryTumor <- data.matrix(BRCA.TP.tibble[,5:ncol(BRCA.TP.tibble)])
row.names(primaryTumor) <-BRCA.TP.tibble$barcode
genes <- changeGeneId(colnames(primaryTumor), from = "Ensembl.ID")

# remove genes with no Hugo Symbol
primaryTumor <- primaryTumor[,!is.na(genes$HGNC.symbol),drop=F]
# remove genes not expressed in at least 20% of the samples
primaryTumor <- 
  primaryTumor[,colSums(primaryTumor>0)>(0.2*nrow(primaryTumor))
           , drop=F]
primaryTumor <- cbind(BRCA.TP.tibble[,1:4],primaryTumor)

TP_patients <- primaryTumor$patient

 All <- data.matrix(BRCA.tibble[,5:ncol(BRCA.tibble)])
 row.names(All) <-BRCA.tibble$barcode
 genes <- changeGeneId(colnames(All))

 # remove genes with no Hugo Symbol
 All <- All[,!is.na(genes$HGNC.symbol),drop=F]
 # remove genes not expressed in at least 20% of the samples
 All <- 
   All[,colSums(All>0)>(0.2*nrow(All))
                , drop=F]
 All <- cbind(BRCA.tibble[,1:4],All)
 

rm(BRCA.TP.tibble,BRCA.tibble, genes)
gc()

load("C:/Users/cifam001/OneDrive - University of South Australia/PhD UNISA/DATA/Hoang's dataset/BRCA_matchedData_full.RData")

HoangPatients <- rownames(BRCA_matchedData$mRNAs)
PseudoHoangData <- primaryTumor%>%
  filter(patient%in%HoangPatients)


#*************************************************************************
#*            Load Single Cell Data


GSE75688_GEO_processed_Breast_Cancer_raw_TPM_matrix <- read_delim("C:/Users/cifam001/OneDrive - University of South Australia/PhD UNISA/DATA/2017 - GSE75688 - Breast Cancer SC/GSE75688_GEO_processed_Breast_Cancer_raw_TPM_matrix.txt", 
                                                                  "\t", escape_double = FALSE, trim_ws = TRUE)
SingleCell <- GSE75688_GEO_processed_Breast_Cancer_raw_TPM_matrix[,-c(1:3)]

GSE75688_final_sample_information <- read_delim("C:/Users/cifam001/OneDrive - University of South Australia/PhD UNISA/DATA/2017 - GSE75688 - Breast Cancer SC/GSE75688_final_sample_information.txt", 
                                                "\t", escape_double = FALSE, trim_ws = TRUE)

TumorSamples <- GSE75688_final_sample_information%>%
#  filter(type =="SC")
 filter(type =="SC",index=="Tumor")

GSE75688_series_matrix <- read_delim("C:/Users/cifam001/OneDrive - University of South Australia/PhD UNISA/DATA/2017 - GSE75688 - Breast Cancer SC/GSE75688_series_matrix.txt", 
                                     "\t", escape_double = FALSE, trim_ws = TRUE, 
                                     skip = 30)

row.names(SingleCell) <- GSE75688_GEO_processed_Breast_Cancer_raw_TPM_matrix$gene_id
SingleCell <-t(SingleCell[,,drop=F])
# remove genes not expressed in at least 30% of the samples
SingleCell <- 
  SingleCell[,colSums(SingleCell>0)>(0.2*nrow(SingleCell))
             , drop=F]

temp <- colnames(SingleCell)
temp <-  gsub("\\..*","",temp)
colnames(SingleCell)<- temp

SingleCell <- SingleCell[TumorSamples$sample,,drop=F]

genes <- changeGeneId(colnames(SingleCell),from = "Ensembl.ID")
# remove genes with no Hugo Symbol
SingleCell <- SingleCell[,!is.na(genes$HGNC.symbol),drop=F]

########################################################################
#*************************************************************************
#*           FIND EVENT = TRUE

#system.time(HER2_95_TP <- findCDbyCI(CancerDS =  primaryTumor,ptimeCov = "HER2"
#                                     , PPIquantile = 0.95,returnModel = F, findEvent=T))

#system.time(Hoang_VIM_95_TP <- findCDbyCI(CancerDS =  PseudoHoangData,ptimeCov = "VIM"
#                                     , PPIquantile = 0.95, findEvent = T))

############ Primary Tumor#####
#----0.7
system.time(TP_HER2_07True <- findCDbyCI(CancerDS =  primaryTumor,ptimeCov = "HER2"
                             , PPIquantile = 0.7,findEvent = T))
save(TP_HER2_07True,file=paste0(dir,"/Experiments/TP_HER2_07True.rda"))

# dummyGeneData <- primaryTumor[,c(10, 1500)]
# dummyGeneData <- dummyGeneData%>%
#   mutate(z=TP_HER2_07True$z, .before = 1)%>%
#   arrange(z)
# 
# dummyCI <- CausalImpact(data = dummyGeneData
#                         , pre.period = c(1,TP_HER2_07True$eventAt-1)
#                         ,post.period = c(TP_HER2_07True$eventAt,1101)
#                         )

system.time(TP_FGF2_07True <- findCDbyCI(CancerDS =  primaryTumor,ptimeCov = "FGF2"
                                         , PPIquantile = 0.7,findEvent = T))
save(TP_FGF2_07True,file=paste0(dir,"/Experiments/TP_FGF2_07True.rda"))

system.time(TP_VIM_07True <- findCDbyCI(CancerDS =  primaryTumor,ptimeCov = "VIM"
                                         , PPIquantile = 0.7,findEvent = T))
save(TP_VIM_07True,file=paste0(dir,"/Experiments/TP_VIM_07True.rda"))

#----0.5
system.time(TP_HER2_05True <- findCDbyCI(CancerDS =  primaryTumor,ptimeCov = "HER2"
                                         , PPIquantile = 0.5,findEvent = T))
save(TP_HER2_05True,file=paste0(dir,"/Experiments/TP_HER2_05True.rda"))

system.time(TP_FGF2_05True <- findCDbyCI(CancerDS =  primaryTumor,ptimeCov = "FGF2"
                                         , PPIquantile = 0.5,findEvent = T))
save(TP_FGF2_05True,file=paste0(dir,"/Experiments/TP_FGF2_05True.rda"))

system.time(TP_VIM_05True <- findCDbyCI(CancerDS =  primaryTumor,ptimeCov = "VIM"
                                        , PPIquantile = 0.5,findEvent = T))
save(TP_VIM_05True,file=paste0(dir,"/Experiments/TP_VIM_05True.rda"))


# # #----0.0
#  system.time(TP_HER2_0True <- findCDbyCI(CancerDS =  primaryTumor,ptimeCov = "HER2"
#                                           , PPIquantile = 0,findEvent = T))
#  save(TP_HER2_0True,file=paste0(dir,"/Experiments/TP_HER2_0True.rda"))
# 
#  system.time(TP_FGF2_0True <- findCDbyCI(CancerDS =  primaryTumor,ptimeCov = "FGF2"
#                                           , PPIquantile = 0,findEvent = T))
#  save(TP_FGF2_0True,file=paste0(dir,"/Experiments/TP_FGF2_0True.rda"))
#  
#  system.time(TP_VIM_0True <- findCDbyCI(CancerDS =  primaryTumor,ptimeCov = "VIM"
#                                          , PPIquantile = 0,findEvent = T))
#  save(TP_VIM_0True,file=paste0(dir,"/Experiments/TP_VIM_0True.rda"))
#  

############ Hoang Dataset ######
#----0.7
system.time(Hoang_HER2_07True <- findCDbyCI(CancerDS =  PseudoHoangData,ptimeCov = "HER2"
                                         , PPIquantile = 0.7,findEvent = T))

save(Hoang_HER2_07True,file=paste0(dir,"/Experiments/Hoang_HER2_07True.rda"))

system.time(Hoang_FGF2_07True <- findCDbyCI(CancerDS =  PseudoHoangData,ptimeCov = "FGF2"
                                        , PPIquantile = 0.7,findEvent = T))
save(Hoang_FGF2_07True,file=paste0(dir,"/Experiments/Hoang_FGF2_07True.rda"))

system.time(Hoang_VIM_07True <- findCDbyCI(CancerDS =  PseudoHoangData,ptimeCov = "VIM"
                                       , PPIquantile = 0.7,findEvent = T))
save(Hoang_VIM_07True,file=paste0(dir,"/Experiments/Hoang_VIM_07True.rda"))

#----0.6
system.time(Hoang_HER2_06True <- findCDbyCI(CancerDS =  PseudoHoangData,ptimeCov = "HER2"
                                            , PPIquantile = 0.6,findEvent = T))

save(Hoang_HER2_06True,file=paste0(dir,"/Experiments/Hoang_HER2_06True.rda"))

system.time(Hoang_FGF2_06True <- findCDbyCI(CancerDS =  PseudoHoangData,ptimeCov = "FGF2"
                                            , PPIquantile = 0.6,findEvent = T))
save(Hoang_FGF2_06True,file=paste0(dir,"/Experiments/Hoang_FGF2_06True.rda"))

system.time(Hoang_VIM_06True <- findCDbyCI(CancerDS =  PseudoHoangData,ptimeCov = "VIM"
                                           , PPIquantile = 0.6,findEvent = T))
save(Hoang_VIM_06True,file=paste0(dir,"/Experiments/Hoang_VIM_06True.rda"))


#----0.5
 system.time(Hoang_HER2_05True <- findCDbyCI(CancerDS =  PseudoHoangData,ptimeCov = "HER2"
                                          , PPIquantile = 0.5,findEvent = T))
 save(Hoang_HER2_05True,file=paste0(dir,"/Experiments/Hoang_HER2_05True.rda"))

 system.time(Hoang_FGF2_05True <- findCDbyCI(CancerDS =  PseudoHoangData,ptimeCov = "FGF2"
                                          , PPIquantile = 0.5,findEvent = T))
 save(Hoang_FGF2_05True,file=paste0(dir,"/Experiments/Hoang_FGF2_05True.rda"))

 system.time(Hoang_VIM_05True <- findCDbyCI(CancerDS =  PseudoHoangData,ptimeCov = "VIM"
                                         , PPIquantile = 0.5,findEvent = T))
 save(Hoang_VIM_05True,file=paste0(dir,"/Experiments/Hoang_VIM_05True.rda"))

  # 
# #----0.0
#  system.time(Hoang_HER2_0True <- findCDbyCI(CancerDS =  PseudoHoangData,ptimeCov = "HER2"
#                                          , PPIquantile = 0,findEvent = T))
#  save(Hoang_HER2_0True,file=paste0(dir,"/Experiments/Hoang_HER2_0True.rda"))
#  
#  system.time(Hoang_FGF2_0True <- findCDbyCI(CancerDS =  PseudoHoangData,ptimeCov = "FGF2"
#                                          , PPIquantile = 0,findEvent = T))
#  save(Hoang_FGF2_0True,file=paste0(dir,"/Experiments/Hoang_FGF2_0True.rda"))
#  
#  system.time(Hoang_VIM_0True <- findCDbyCI(CancerDS =  PseudoHoangData,ptimeCov = "VIM"
#                                         , PPIquantile = 0,findEvent = T))
#  save(Hoang_VIM_0True,file=paste0(dir,"/Experiments/Hoang_VIM_0True.rda"))
# 


########### Single cell ##########


#---- 0.7
system.time(SC_HER2_07True <- findCDbyCI(CancerDS =  SingleCell ,ptimeCov = "HER2"
                                   ,PPIquantile = 0.7, findEvent = T, Step = 1))
save(SC_HER2_07True,file=paste0(dir,"/Experiments/SC_HER2_07True.rda"))



system.time(SC_VIM_07True <- findCDbyCI(CancerDS =  SingleCell ,ptimeCov = "VIM"
                                        ,PPIquantile = 0.7, findEvent = T))
save(SC_VIM_07True,file=paste0(dir,"/Experiments/SC_VIM_07True.rda"))

#---- 0.6
system.time(SC_HER2_06True <- findCDbyCI(CancerDS =  SingleCell ,ptimeCov = "HER2"
                                         ,PPIquantile = 0.6, findEvent = T, Step = 1))
save(SC_HER2_06True,file=paste0(dir,"/Experiments/SC_HER2_06True.rda"))



system.time(SC_VIM_06True <- findCDbyCI(CancerDS =  SingleCell ,ptimeCov = "VIM"
                                        ,PPIquantile = 0.6, findEvent = T))
save(SC_VIM_06True,file=paste0(dir,"/Experiments/SC_VIM_06True.rda"))



#---- 0.5
system.time(SC_HER2_05True <- findCDbyCI(CancerDS =  SingleCell ,ptimeCov = "HER2"
                                         ,PPIquantile = 0.5, findEvent = T, Step = 1))
save(SC_HER2_05True,file=paste0(dir,"/Experiments/SC_HER2_05True.rda"))



system.time(SC_VIM_05True <- findCDbyCI(CancerDS =  SingleCell ,ptimeCov = "VIM"
                                        ,PPIquantile = 0.5, findEvent = T))
save(SC_VIM_05True,file=paste0(dir,"/Experiments/SC_VIM_05True.rda"))

##---- 0.0
# system.time(SC_HER2_0True <- findCDbyCI(CancerDS =  SingleCell ,ptimeCov = "HER2"
#                                          ,PPIquantile = 0, findEvent = T))
# save(SC_HER2_0True,file=paste0(dir,"/Experiments/SC_HER2_0True.rda"))
#  
# 
# 
# system.time(SC_VIM_0True <- findCDbyCI(CancerDS =  SingleCell ,ptimeCov = "VIM"
#                                         ,PPIquantile = 0, findEvent = T))
# save(SC_VIM_0True,file=paste0(dir,"/Experiments/SC_VIM_0True.rda"))
# 
# 
# ##########################################################################
# #*************************************************************************
# #*           FIND EVENT = FALSE
# 
# ############ Primary Tumor#####
# #---- 0.7
#  system.time(TP_HER2_07False <- findCDbyCI(CancerDS =  primaryTumor,ptimeCov = "HER2"
#                                           , PPIquantile = 0.7))
#  save(TP_HER2_07False,file=paste0(dir,"/Experiments/TP_HER2_07False.rda"))
# 
#  system.time(TP_FGF2_07False <- findCDbyCI(CancerDS =  primaryTumor,ptimeCov = "FGF2"
#                                            , PPIquantile = 0.7))
#  save(TP_FGF2_07False,file=paste0(dir,"/Experiments/TP_FGF2_07False.rda"))
# 
#  system.time(TP_VIM_07False <- findCDbyCI(CancerDS =  primaryTumor,ptimeCov = "VIM"
#                                            , PPIquantile = 0.7))
#  save(TP_VIM_07False,file=paste0(dir,"/Experiments/TP_VIM_07False.rda"))
# 
#  
# #---- 0.5
#  system.time(TP_HER2_05False <- findCDbyCI(CancerDS =  primaryTumor,ptimeCov = "HER2"
#                                            , PPIquantile = 0.5))
#  save(TP_HER2_05False,file=paste0(dir,"/Experiments/TP_HER2_05False.rda"))
# 
#  system.time(TP_FGF2_05False <- findCDbyCI(CancerDS =  primaryTumor,ptimeCov = "FGF2"
#                                          , PPIquantile = 0.5)) 
#  save(TP_FGF2_05False,file=paste0(dir,"/Experiments/TP_FGF2_05False.rda"))
# 
#  system.time(TP_VIM_05False <- findCDbyCI(CancerDS =  primaryTumor,ptimeCov = "VIM"
#                                           , PPIquantile = 0.5))
#  save(TP_VIM_05False,file=paste0(dir,"/Experiments/TP_VIM_05False.rda"))
# 
# #---- 0.0
#  
#  system.time(TP_HER2_0False <- findCDbyCI(CancerDS =  primaryTumor,ptimeCov = "HER2"
#                                            , PPIquantile = 0))
#  save(TP_HER2_0False,file=paste0(dir,"/Experiments/TP_HER2_0False.rda"))
# 
#  system.time(TP_FGF2_0False <- findCDbyCI(CancerDS =  primaryTumor,ptimeCov = "FGF2"
#                                            , PPIquantile = 0))
#  save(TP_FGF2_0False,file=paste0(dir,"/Experiments/TP_FGF2_0False.rda"))
# 
#  system.time(TP_VIM_0False <- findCDbyCI(CancerDS =  primaryTumor,ptimeCov = "VIM"
#                                           , PPIquantile = 0))
#  save(TP_VIM_0False,file=paste0(dir,"/Experiments/TP_VIM_0False.rda"))
#  
#  
# ############ Hoang's Dataset#####
#  #---- 0.7
#  system.time(Hoang_HER2_07False <- findCDbyCI(CancerDS =  PseudoHoangData,ptimeCov = "HER2"
#                                            , PPIquantile = 0.7))
#  save(Hoang_HER2_07False,file=paste0(dir,"/Experiments/Hoang_HER2_07False.rda"))
#  
#  system.time(Hoang_FGF2_07False <- findCDbyCI(CancerDS =  PseudoHoangData,ptimeCov = "FGF2"
#                                            , PPIquantile = 0.7))
#  save(Hoang_FGF2_07False,file=paste0(dir,"/Experiments/Hoang_FGF2_07False.rda"))
#  
#  system.time(Hoang_VIM_07False <- findCDbyCI(CancerDS =  PseudoHoangData,ptimeCov = "VIM"
#                                           , PPIquantile = 0.7))
#  save(Hoang_VIM_07False,file=paste0(dir,"/Experiments/Hoang_VIM_07False.rda"))
# 
#   
# #---- 0.5
#  system.time(Hoang_HER2_05False <- findCDbyCI(CancerDS =  PseudoHoangData,ptimeCov = "HER2"
#                                           , PPIquantile = 0.5))
#  save(Hoang_HER2_05False,file=paste0(dir,"/Experiments/Hoang_HER2_05False.rda"))
# 
#  system.time(Hoang_FGF2_05False <- findCDbyCI(CancerDS =  PseudoHoangData,ptimeCov = "FGF2"
#                                           , PPIquantile = 0.5))
#  save(Hoang_FGF2_05False,file=paste0(dir,"/Experiments/Hoang_FGF2_05False.rda"))
# 
#  system.time(Hoang_VIM_05False <- findCDbyCI(CancerDS =  PseudoHoangData,ptimeCov = "VIM"
#                                          , PPIquantile = 0.5))
#  save(Hoang_VIM_05False,file=paste0(dir,"/Experiments/Hoang_VIM_05False.rda"))
# 
#  
# #---- 0.0
#  system.time(Hoang_HER2_0False <- findCDbyCI(CancerDS =  PseudoHoangData,ptimeCov = "HER2"
#                                           , PPIquantile = 0))
#  save(Hoang_HER2_0False,file=paste0(dir,"/Experiments/Hoang_HER2_0False.rda"))
# 
#  system.time(Hoang_FGF2_0False <- findCDbyCI(CancerDS =  PseudoHoangData,ptimeCov = "FGF2"
#                                           , PPIquantile = 0))
#  save(Hoang_FGF2_0False,file=paste0(dir,"/Experiments/Hoang_FGF2_0False.rda"))
# 
#  system.time(Hoang_VIM_0False <- findCDbyCI(CancerDS =  PseudoHoangData,ptimeCov = "VIM"
#                                          , PPIquantile = 0))
#  save(Hoang_VIM_0False,file=paste0(dir,"/Experiments/Hoang_VIM_0False.rda"))
# 
#  ############ Single Cell Dataset#####
# #---- 0.7
# system.time(SC_HER2_07False <- findCDbyCI(CancerDS =  SingleCell ,ptimeCov = "HER2"
#                                         ,PPIquantile = 0.7))
# save(SC_HER2_07False,file=paste0(dir,"/Experiments/SC_HER2_07False.rda"))
# 
# 
# 
#  system.time(SC_VIM_07False <- findCDbyCI(CancerDS =  SingleCell ,ptimeCov = "VIM"
#                                           ,PPIquantile = 0.7))
#  save(SC_VIM_07False,file=paste0(dir,"/Experiments/SC_VIM_07False.rda"))
# 
#  
# #---- 0.5
#  system.time(SC_HER2_05False <- findCDbyCI(CancerDS =  SingleCell ,ptimeCov = "HER2"
#                                            ,PPIquantile = 0.5))
#  save(SC_HER2_05False,file=paste0(dir,"/Experiments/SC_HER2_05False.rda"))
# 
# 
# 
#  system.time(SC_VIM_05False <- findCDbyCI(CancerDS =  SingleCell ,ptimeCov = "VIM"
#                                          ,PPIquantile = 0.5))
#  save(SC_VIM_05False,file=paste0(dir,"/Experiments/SC_VIM_05False.rda"))
# 
# #---- 0.0
#  system.time(SC_HER2_0False <- findCDbyCI(CancerDS =  SingleCell ,ptimeCov = "HER2"
#                                            ,PPIquantile = 0))
#  save(SC_HER2_0False,file=paste0(dir,"/Experiments/SC_HER2_0False.rda"))
#  
# 
#  system.time(SC_VIM_0False <- findCDbyCI(CancerDS =  SingleCell ,ptimeCov = "VIM"
#                                           ,PPIquantile = 0))
#  save(SC_VIM_0False,file=paste0(dir,"/Experiments/SC_VIM_0False.rda"))
#  
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# ###########   
# 
#  
#  
#  ##########################################################################
# ############ ER status   ######
# 
# DS <- PseudoHoangData%>%
#   filter(er_status_by_ihc%in%c("Negative","Positive"))
# 
#  system.time(Hoang_ERstatus_0 <- findCDbyCI(CancerDS =  DS,ptimeCov = "er_status_by_ihc"
#                                       ,PPIquantile = 0))
#  save(Hoang_ERstatus_0,file=paste0(dir,"/Experiments/Hoang_ERstatus_0.rda"))
# 
# 
#  system.time(Hoang_ERstatus_05 <- findCDbyCI(CancerDS =  DS,ptimeCov = "er_status_by_ihc"
#                                        ,PPIquantile = 0.5))
#  save(Hoang_ERstatus_05,file=paste0(dir,"/Experiments/Hoang_ERstatus_05.rda"))
# 
# 
#  system.time(Hoang_ERstatus_07 <- findCDbyCI(CancerDS =  DS,ptimeCov = "er_status_by_ihc"
#                                        ,PPIquantile = 0.7))
#  save(Hoang_ERstatus_07,file=paste0(dir,"/Experiments/Hoang_ERstatus_07.rda"))
# 
# 
#  DS <- primaryTumor%>%
#    filter(er_status_by_ihc%in%c("Negative","Positive"))
# 
#  system.time(TP_ERstatus_0 <- findCDbyCI(CancerDS =  DS,ptimeCov = "er_status_by_ihc"
#                                          ,PPIquantile = 0))
#  save(TP_ERstatus_0,file=paste0(dir,"/Experiments/TP_ERstatus_0.rda"))
# 
# 
#  system.time(TP_ERstatus_05 <- findCDbyCI(CancerDS =  DS,ptimeCov = "er_status_by_ihc"
#                                           ,PPIquantile = 0.5))
#  save(TP_ERstatus_05,file=paste0(dir,"/Experiments/TP_ERstatus_05.rda"))
# 
# 
#  system.time(TP_ERstatus_07 <- findCDbyCI(CancerDS =  DS,ptimeCov = "er_status_by_ihc"
#                                           ,PPIquantile = 0.7))
#  save(TP_ERstatus_07,file=paste0(dir,"/Experiments/TP_ERstatus_07.rda"))
# 
# 
# # covariant <- GSE75688_series_matrix%>%
# #   select(all_of(TumorSamples$sample))
# # covariant <- covariant[10,]
# # covariant <- t(covariant)
# # 
# # DS <- cbind(covariant[row.names(SingleCell),],SingleCell)
# # system.time(categorical_07_SC <- findCDbyCI(CancerDS =  DS ,ptimeCov = "V1"
# #                                     ,PPIquantile = 0.7))
# # save(categorical_07_SC,file=paste0(dir,"/Experiments/categorical_07_SC.rda"))
# # 
# # 
# # system.time(categorical_05_SC <- findCDbyCI(CancerDS =  DS ,ptimeCov = "V1"
# #                                            ,PPIquantile = 0.5))
# # save(categorical_05_SC,file=paste0(dir,"/Experiments/categorical_05_SC.rda"))
# # 
# # #system.time(categorical_0_SC <- findCDbyCI(CancerDS =  DS ,ptimeCov = "V1"
# # #                                            ,PPIquantile = 0))
# # 
# # 
