rm(list = ls())


library(readxl)
library(AMCBGeneUtils)
library(tidyverse)

###############-------load experiments-------#############
path <- "C:/Users/cifam001/OneDrive - University of South Australia/PhD UNISA/CancerDrivers/Causal Impact/R/Experiments"
path <- paste0(path,"/RelEffect(Abs)")
load("C:/Users/cifam001/OneDrive - University of South Australia/PhD UNISA/DATA/CGC.driverNames.rdata")

file_list <- list.files(path=path,all.files = F, recursive = F,include.dirs = F)
file_list <- setdiff(file_list,list.dirs(path=path,full.names = F))
index<- str_detect(file_list
                   , pattern = "SC.+0[6|7]True.+")

file_list <- file_list[index]


for (i in file_list) {
  load(paste0(path,"/",i))
}
experiments <- str_remove(file_list, ".rda")

dataVIM <- SC_HER2_06True$res$CDinfer%>%
  arrange(desc(rank),desc(RelEffect))

intersect(dataVIM$Ensembl.ID[1:50],CGC.driverNames$Ensembl.ID)

###############-------load pseudoGT (literature)-------#############
pseudoGT.dir <- "C:/Users/cifam001/OneDrive - University of South Australia/PhD UNISA/DATA/CDpseudoGT"

CDpapers <- c("Lawrence2014","Martincorena2017"
              ,"Bailey2018","Priestley2019","Dietlein2020"
              ,"Rheinbay2020")

#----Nature_Lawrence 505 pages495-501(2014)
Lawrence2014 <- read_excel(paste0(pseudoGT.dir
                                        , "/Lawrence2014.xlsx")
                                  , sheet = 1, range = "A2:I262")
Lawrence2014 <- Lawrence2014%>%
  dplyr::select(-3)

Lawrence2014 <- cbind(changeGeneId(Lawrence2014$gene)
                             ,Lawrence2014[,2:8])

#----CELL_Martincorena  171 (2017).xlsx
Martincorena2017  <- read_excel(paste0(pseudoGT.dir
                                        , "/Martincorena2017.xlsx")
                                 , sheet = 1)
Martincorena2017 <- changeGeneId(Martincorena2017$Gene)


#------Cell_Bailey_Vol 173 Issue 2 Pg 371-385.e18(2018)
Bailey2018 <-read_excel(paste0(pseudoGT.dir
                                    , "/Bailey2018.xlsx")
                             , sheet = 2, range = "A4:L743")

Bailey2018 <- cbind(changeGeneId(Bailey2018$Gene)
                        ,Bailey2018[,2:12])
BRCA <- list()
BRCA$Bailey2018 <- Bailey2018%>%
  filter(Cancer == "BRCA")

#-----Nature_Priestley Vol 575_7(2019)
Priestley2019 <- read_excel(paste0(pseudoGT.dir
                                          , "/Priestley2019.xlsx")
                                   , sheet = 1, range = "A1:L20071")
Priestley2019 <- cbind(changeGeneId(Priestley2019$gene)
                              ,Priestley2019[,-3])
BRCA$Priestley2019 <- Priestley2019%>%
  filter(cancerType == "Breast")

#------Nature Genetics_Dietlein_Vol 52_pg208-218 (2020)
Dietlein2020 <-read_excel(paste0(pseudoGT.dir
                                    , "/Dietlein2020.xlsx")
                             , sheet = 4, range = "A1:H828")

Dietlein2020 <- cbind(changeGeneId(Dietlein2020$Gene)
                         ,Dietlein2020[,2:8])
BRCA$Dietlein2020 <- Dietlein2020%>%
  filter(`Cancer Type` == c("Breast"))


#------Nature_Rheinbay Vol 578, pages102-111 (2020)
Rheinbay2020 <- read_excel(paste0(pseudoGT.dir
                                         , "/Rheinbay2020.xlsx")
                                  , sheet = 4)
Rheinbay2020 <-  cbind(changeGeneId(Rheinbay2020$gene), Rheinbay2020[,-4])

BRCA$Rheinbay2020 <- Rheinbay2020%>%
  filter(tissue %in%c("Breast-AdenoCa","meta_Breast"))

###############-------Create stratification table-------#############
genes <- rbind(SC_HER2_06True$res$CDinfer
               ,SC_VIM_06True$res$CDinfer)

genes <- unique(genes$Ensembl.ID)
genes <- changeGeneId(genes , from="Ensembl.ID")

DriverCatalogues <- vector(mode="list", length = length(CDpapers))
names(DriverCatalogues) <- CDpapers
aux <- character(0)
for (i in CDpapers) {
  DriverCatalogues[[i]] <- unique(get(i)$HGNC.symbol)
  aux <-union(aux,DriverCatalogues[[i]])
  print(length(intersect(genes$HGNC.symbol,DriverCatalogues[[i]]))/length(DriverCatalogues[[i]]))
}

consistency <- data.frame(HGNC.symbol=aux)
for (i in CDpapers) {
  consistency <- consistency%>%
    mutate("{i}" := aux %in% DriverCatalogues[[i]])
}
consistency <- consistency%>%
  mutate(n=rowSums(across(Lawrence2014:Rheinbay2020)))%>%
  mutate(inCGC = HGNC.symbol%in%CGC.driverNames$HGNC.symbol)

AdditionalCGC <- data.frame(Study=c("SC_combined",CDpapers), n=0)
row.names(AdditionalCGC) <- c("SC_combined",CDpapers)

a <- intersect(CGC.driverNames$HGNC.symbol, genes$HGNC.symbol)
b <- intersect(CGC.driverNames$HGNC.symbol, consistency$HGNC.symbol)
myNovelCGC <- setdiff(a,b)
AdditionalCGC["SC_combined","n"] <- length(myNovelCGC)

for (i in CDpapers) {
  aux <- consistency%>%
    filter(inCGC==T, n==1)%>%
    dplyr::select(all_of(i))
  AdditionalCGC[i,"n"] <- sum(aux)
}

######---------- GGPLOT---------------------
AdditionalCGC$Study <- factor(AdditionalCGC$Study
                              , levels = c(CDpapers,"SC_combined"))
p <- AdditionalCGC %>%
  ggplot(aes(x=Study, y=n, fill=Study)) +
  geom_bar(stat="identity")+
  geom_text(aes(label=n)
            , vjust=-1, color="black")+
  scale_fill_brewer(palette="RdYlBu")+
  theme_minimal()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

p + labs(x="", y = "Unique CGC genes")

SC_HER2_06 <- cbind(changeGeneId(SC_HER2_06True$res$CDinfer$Ensembl.ID
                                , from = "Ensembl.ID")[,-1]
                   ,SC_HER2_06True$res$CDinfer[,-1])

write.csv(SC_HER2_06, file = "SC_HER2_06.csv")

SC_HER2_07 <- cbind(changeGeneId(SC_HER2_07True$res$CDinfer$Ensembl.ID
                                 , from = "Ensembl.ID")[,-1]
                    ,SC_HER2_07True$res$CDinfer[,-1])

write.csv(SC_HER2_07, file = "SC_HER2_07.csv")


SC_VIM_06 <- cbind(changeGeneId(SC_VIM_06True$res$CDinfer$Ensembl.ID
                                 , from = "Ensembl.ID")[,-1]
                    ,SC_VIM_06True$res$CDinfer[,-1])
SC_VIM_07 <- cbind(changeGeneId(SC_VIM_07True$res$CDinfer$Ensembl.ID
                                , from = "Ensembl.ID")[,-1]
                   ,SC_VIM_07True$res$CDinfer[,-1])

write.csv(SC_VIM_06, file = "SC_VIM_06.csv")

SC_VIM_HER.union <- rbind(SC_VIM_06,SC_HER2_06)%>%
  distinct(Ensembl.ID, .keep_all = T)%>%
  arrange(desc(abs(RelEffect)))

SC_VIM_HER.union <- rbind(SC_VIM_06,SC_HER2_06)%>%
  distinct(Ensembl.ID, .keep_all = T)%>%
  arrange(desc(abs(RelEffect)))


SC_VIM_HER.overlap <- SC_VIM_HER.union%>%
  filter(Ensembl.ID%in%intersect(SC_VIM_06$Ensembl.ID
                                 ,SC_HER2_06$Ensembl.ID))%>%
  dplyr::select(-rank)%>%
  mutate(inCGC = Ensembl.ID%in%CGC.driverNames$Ensembl.ID)%>%
  arrange(desc(abs(RelEffect)))

path <- "C:/Users/cifam001/OneDrive - University of South Australia/PhD UNISA/CancerDrivers/Causal Impact/"
write.csv(SC_VIM_HER.overlap, file = paste0(path,"SC_VIM_HER.overlap.csv"))

SC_VIM_HER.table <- matrix(data= "",nrow = 604, ncol = 3)
colnames(SC_VIM_HER.table) <- c("VIMtime","HER2time","CGC")

aux <- SC_VIM_06$HGNC.symbol
SC_VIM_HER.table[1:length(aux),1] <- aux
aux <- SC_HER2_06$HGNC.symbol
SC_VIM_HER.table[1:length(aux),2] <- aux
aux <- intersect(SC_VIM_HER.table[,1:2],CGC.driverNames$HGNC.symbol)
SC_VIM_HER.table[1:length(aux),3] <- aux

write.csv(SC_VIM_HER.table, file = paste0(path,"SC_VIM_HER.table.csv"))

