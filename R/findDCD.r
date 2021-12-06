#' @title findDCD
#'
#' @description findDCD identifies genes driving driving one (or more)
#' significant biological processes along cancer progression based on the
#' hypothesis that the causal relationship between a cancer driver
#' gene and cancer development induces a significant deviation (also referred
#' as causal impact) of a core process from normal to carcinogenic.
#'
#' @usage function(GeneExpression, z=NULL, pathCovariate =NULL
#' , findEvent = T, Step=1, chunk_size= 100
#' , PPItop = 0.3, alpha=0.05, CIniter=200
#' , returnModel=F, elbo_tol=1e-3, project = NULL)\cr
#'
#' @param GeneExpression A \code{matrix} containing mRNA gene expression.\cr
#' Columns represent mRNAs, rows represent samples. Column names can be any of
#' the following: Hugo Gene Symbol (HGNC.symbol), Hugo Gene ID (HGNC.ID),
#' EnsemblID, NBCI.
#' @param z A \code{numeric vector} containing a pseudotime
#' score for each sample. If \code{NULL} (default) a pseudotime score is calculated
#' by using \code{Phenopath} package and the \code{pathCovariate}.
#' @param findEvent If \code{TRUE} (default) samples are ordered in pseudotime order
#' and deviations from normal to cancerogenic are assessed by using
#' \code{CausalImpact} from the package \code{CausalImpact}. The sample with the
#' largest (significant) \{CausalImpact} is labelled as the "event". If \code{FALSE},
#' the sample where the change of sign (from negative to positive) occurs is
#' laballed as the "event".
#' @param Step An \code{integer} indicating the distance between samples to be assesed
#' when \code{findEvent = TRUE}. \code{Step} = 1 (default) means all samples are
#' considered during the \code{findEvent} process.
#' @param chunk_size An \code{integer} defining the number of genes to be analysed
#' at a time. \code{chunk_size = 100} (default) indicates that groups of 100 genes
#' will be analised at a time.
#' @param PPItop A \{numeric} value between 0 and 1 indicating the percentage of
#' PPI genes in the dataset to be selected as putative drivers. PPI genes with the
#' most interactions are selected.
#' @param alpha Significance level for the statistical test.
#' \code{alpha=0.05} by default.
#'#' @param CIniter number of iterations (200 by default) for
#' \code{CausalImpact} modeling.
#' @param returnModel If \code{TRUE}, it includes the complete \code{CausalImpact} model in
#' the outcome of \code{findDCD}. If \code{FALSE} (default), only the most relevant
#' parameters of the \code{CausalImpact} model are returned.
#' @param elbo_tol A \code{numeric} value for the tolerance during the Pseudotime
#' scoring by \code{Phenopath}. The lower the tolerance the stricter the scoring.
#' \code{elbo_tol = 1e-3} by default.
#' @param project An optional parameter with a TCGA project name (e.g. BRCA).
#' If provided, a dummy rank for the inferred dynamic cancer driver is calculated
#' based on the frequency of mutations of those genes in the TCGA project dataset.
#'
#'
#' @author Andres Mauricio Cifuentes_Bernal, Vu VH Pham, Xiaomei Li, Lin Liu, JiuyongLi and Thuc Duy Le
#' @export
#' @seealso \link[DynamicCancerDriver]{findCovariate},
#' \link[DynamicCancerDriver]{parallelCI}
#'
#' @return A \code{list} consisting of the following elements:
#'   \item{\code{res}}{A \code{list} with the results of the DynamicCancerDriver
#'   inference process. Results are listed as follows:
#'   \enumerate{
#'           \item{\code{FS:}}{A \code{vector} containing the names of the putative
#'           cancer drivers}
#'           \item{\code{CausalImpact:}}{Causal impact models of the putative drivers}
#'           \item{\code{CDinfer:}}{Inferred Dynamic Cancer Drivers}
#'           \item{\code{summary:}}{A table with a summary of the results}
#'           }
#'    For each target gene with at least one parent. The index of the parents.}
#'   \item{\code{eventAt}}{A \code{integer} containing the index (after pseudotime
#'   order) of the sample labelled as the "event".}
#'   \item{\code{z}}{Pseudotime score}
#'
#' @examples \dontrun{
#'    data("GSE75688_TPM_tumor", package = "DynamicCancerDriver")
#'
#' ----- Find Dynamic Cancer Drivers, PPI top 40% -----
#' DCD.HER2time_SC <- findDCD(GeneExpression = GSE75688_TPM_tumor
#'                            , pathCovariate = "HER2"
#'                            , PPItop = 0.3
#'                            , findEvent = TRUE)
#' }
#'
#' @references
#'
#'


findDCD <- function(GeneExpression, z=NULL, pathCovariate =NULL
                    , findEvent = T, Step=1, chunk_size= 100
                    , PPItop = 0.3, alpha=0.05, CIniter=200
                    , returnModel=F, elbo_tol=1e-3, project = NULL){

  #---- Function for validating pseudotime vector (if provided) ----
  z_validate <- function(z){
    valid <- FALSE
    if(!is.null(z)){
      z <- as.matrix(z,ncol=1)
      if(length(z)%%nrow(GeneExpression)==0){
        valid <- all(sapply(z, is.numeric))
      }
    }else{
      valid <- TRUE
    }

    if(!valid){
      message("z needs to be NULL or a numeric vector with length == nrow(GeneExpression)")
      message("calculating z using pathCovariate")
      z <- NULL
    }
    return(z)
  }

  #---- pathCovariate_validate: Function to verify if the pathCovariate is valid  ----
  pathCovariate_validate <- function(covariate){
    covariate <- as.matrix(covariate)
    n <- nrow(covariate)
    isNumerical <- all(sapply(covariate, is.numeric))
    if(!isNumerical){
      intersection <- intersect(covariate,colnames(GeneExpression))
      if(length(intersection)<1){
        aux <- covariate
        isEnsembl <- str_starts(aux, "ENSG0", negate = FALSE)
        if(!all(isEnsembl)){
          aux[!isEnsembl] <- changeGeneId(covariate[!isEnsembl]
                                          ,from = "HGNC.symbol"
                                          ,to = "Ensembl.ID")[[2]]
        }
        intersection <- intersect(aux,colnames(GeneExpression))
        if(length(intersection)<1){
          if(n!=nrow(GeneExpression))
            stop("invalid covariate")
        }
      }
      covariate <- GeneExpression[,intersection, drop=F]
      row.names(covariate) <- row.names(GeneExpression)
    }else{
      if(n != nrow(GeneExpression))
        stop("invalid covariate")
    }
    return(covariate)
  }
  #---- If path covariate is categorical, event when pseudotime change its sign ----
  eventAtZero <- function(z){
    z <- as.matrix(z,ncol=1)
    eventAt <- which(z[order(data.matrix(z)),1] >= 0)[1]
    return(list(eventAt=eventAt, z=z))
  }

  #---- function to find the event (only for numerical path covariates)
  DCD.findEvent <- function(z, GeneExpression, pathCovariate.name, Step=Step){
    sControl <- findCovariate(GeneExpression = GeneExpression
                                     , FS = pathCovariate.name)

    if(!require(parallel))
      install.packages("parallel")
    if(!require(CausalImpact))
      install.packages("CausalImpact")


    library(parallel)
    library(CausalImpact)


    DS.order <- GeneExpression%>%
      dplyr::select(all_of(sControl[,1]),all_of(sControl[,2]))%>%
      mutate(z=z, .before = 1)%>%
      arrange(z)
    DS.order <- data.matrix(DS.order)

    n<-nrow(DS.order)
    toTest <- seq.int(5+Step,n-5,by=Step)

    CausalImp<-vector(mode = "list",length = length(toTest))
    names(CausalImp) <- toTest

    copies_of_r <- detectCores(logical = F)-1

    chunk <- split(toTest, ceiling(toTest/chunk_size))
    cl <- makeCluster(copies_of_r)
    clusterExport(cl, c("DS.order","CausalImpact", "zoo", "sControl")
                  ,envir = environment())

    for (i in 1:length(chunk)) {
      message(paste0("------ (DCD.findEvent) chunk ",toString(i), " of "
                     , toString(length(chunk))
                     ," ------"))
      CausalImp[as.character(chunk[[i]])] <-
        parSapply(cl, chunk[[i]]
                  , function(a){

                    set.seed(1)
                    res <- CausalImpact(data = zoo(DS.order[,2:3]
                                                   ,order.by = DS.order[,1])
                                        , pre.period = DS.order[c(1,a-1),1, drop=T]
                                        , post.period =  DS.order[c(a,nrow(GeneExpression)),1, drop=T]
                                        , model.args = list(niter = CIniter) )

                    res <- cbind(c(a,a),res[["summary"]])
                    colnames(res)[1] <- "sample"
                    return(res)
                  }, simplify = F)
      gc()
    }# end for
    stopCluster(cl)

    res <- which(sapply(CausalImp
                        , function(x){all(x["p"]<alpha)}))

    res <- sapply(CausalImp[res],function(x){abs(x[1,"RelEffect"])})

    res <- names(res)[which(res==max(res))]
    return(list(eventAt=as.numeric(res[1]), z=z))
  }

  #---- function to find pseudotime score (based on Phenopath)
  findZ <- function(GeneExpression, FS, pathCovariate, elbo_tol = elbo_tol){
    if(!require(phenopath)){
      if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")

      BiocManager::install("phenopath")
    }
    library(phenopath)
    exprs_obj <- GeneExpression%>%
      dplyr::select(all_of(FS))%>%
      transmute_all(as.numeric)

    exprs_obj <- log(exprs_obj + 1)

    Ptime<-phenopath(exprs_obj = data.matrix(exprs_obj)
                     , x = data.matrix(pathCovariate), elbo_tol = elbo_tol)

    #plot_elbo(Ptime)
    z <- as.matrix(trajectory(Ptime),ncol=1)
    row.names(z)<-row.names(GeneExpression)
    colnames(z)<-"z"
    return(z)
  }

  #############--------  MAIN FUNCTION CODE  --------#############
  GeneExpression <-dplyr::as_tibble(GeneExpression, rownames=NA)  #to keep rownames
  PPIrank <- CreatePPIRank(colnames(GeneExpression))
  temp <- PPIrank%>%
    filter(total>=quantile(total,1-PPItop))

  FS <- temp[[1]]
  rm(temp)

  case <- is.null(z)*4+is.null(pathCovariate)*2+isTRUE(findEvent)+1

  switch(case
         #case 1
         ,{z <- z_validate(z)
         event <- eventAtZero(z)}
         #case 2
         ,{z <- z_validate(z)
         pathCovariate <- pathCovariate_validate(pathCovariate)
         pathCovariate.isNum <-  all(sapply(pathCovariate, is.numeric))
         if(pathCovariate.isNum)
           event <- DCD.findEvent(z = z
                                  ,GeneExpression = GeneExpression
                                  ,pathCovariate.name = colnames(pathCovariate)
                                  ,Step =Step)
         else
           event <- eventAtZero(z)
         }
         #case 3
         , event <- eventAtZero(z)
         #case 4
         , event <- eventAtZero(z)
         #case 5
         ,{pathCovariate <- pathCovariate_validate(pathCovariate)
         z <- findZ(GeneExpression = GeneExpression, FS =  FS
                    , pathCovariate=pathCovariate,elbo_tol = elbo_tol)
         event <- eventAtZero(z)}
         #case 6
         ,{pathCovariate <- pathCovariate_validate(pathCovariate)
         z <- findZ(GeneExpression = GeneExpression, FS =  FS
                    , pathCovariate=pathCovariate, elbo_tol = elbo_tol)
         pathCovariate.isNum <-  all(sapply(pathCovariate, is.numeric))
         if(pathCovariate.isNum)
           event <- DCD.findEvent(z = z
                                  ,GeneExpression = GeneExpression
                                  ,pathCovariate.name = colnames(pathCovariate)
                                  ,Step =Step)
         else
           event <- eventAtZero(z)}
         #default
         ,message("Error: z and pathCovariate missing")
  )

  sControl <- findCovariate(GeneExpression=GeneExpression, FS=FS, PPIrank)
  CausalImp <- parallelCI(GeneExpression = GeneExpression
                          ,sControl = sControl
                          ,z = z
                          ,chunk_size = 50
                          ,eventAt=event$eventAt
                          ,CIniter = CIniter
                          ,returnModel = returnModel)
  if(returnModel){
    CD <- names(which(sapply(CausalImp
                             , function(x){all(x[["summary"]]["p"]<alpha)})))
  }else{
    CD <- names(which(sapply(CausalImp
                             , function(x){all(x["p"]<alpha)})))
  }

  if(!is.null(project)){
    if("patient"%in%colnames(GeneExpression)){
      CDrank <- rankByMut(genesIds = CD
                          ,project = project
                          ,patient_IDs = as.matrix(GeneExpression[,"patient"]))

    }else{
      CDrank <- rankByMut(genesIds = CD
                          ,project = project)
    }
    CDinfer <-  CDrank%>%
      mutate(Ensembl.ID = row.names(CDrank),.before=1)
  }else{
    CDinfer <- data.frame(Ensembl.ID = CD)
  }



  temp <- sapply(CausalImp[CD]
                 , function(x){
                   return(x[1,c("RelEffect","AbsEffect","p")])
                 },simplify = FALSE)

  temp <- do.call(what = rbind,temp)
  temp <- temp%>%
    rownames_to_column(var = "Ensembl.ID")
  CDinfer <- left_join(x = CDinfer, y = temp,by="Ensembl.ID")

  if(!is.null(project)){
    CDinfer <- CDinfer%>%
      arrange(desc(rank),desc(RelEffect))
  }else{
    CDinfer <- CDinfer%>%
      arrange(desc(RelEffect))
  }

  CDinfer <- cbind(AMCBGeneUtils::changeGeneId(CDinfer[,1], from = "Ensembl.ID")[2:4]
                   ,CDinfer[,-1])
  rename(CDinfer, p = p.val)



  data(CGC.driverNames)
  nCD <- length(CD)
  nCDconf <- length(intersect(CGC.driverNames$Ensembl.ID
                              ,CD))
  nFSconf <- length(intersect(CGC.driverNames$Ensembl.ID
                              ,FS))

  summary <- matrix(
    c("Covariant(s) ptime: " ,colnames(pathCovariate)
      ,"Top PPI" ,PPItop
      ,"Find Event", paste0(findEvent,"(event at:",event$eventAt,")")
      ,"nFS" ,length(FS)
      ,"nFSconf", nFSconf
      ,"inferred" ,nCD
      ,"Confirmed" ,nCDconf
      ,"percentage",100*nCDconf/nCD
    ), ncol =2, byrow = T)


  res <- list()
  res$FS <- FS
  res$CausalImpact <- CausalImp[CD]
  res$CDinfer <- CDinfer
  res$summary <- summary
  outcome <- list()
  outcome$res <- res
  outcome$eventAt <- event$eventAt
  outcome$z <- event$z
  return(outcome)


}
