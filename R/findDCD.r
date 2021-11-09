findDCD <- function(GeneExpression, z=NULL, pathCovariate =NULL
                    , findEvent = T, Step=1, chunk_size= 100
                    , PPItop = 0.3, alpha=0.05, CIniter=200
                    , returnModel=F, elbo_tol=1e-3){

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
    sControl <- findSyntheticControl(GeneExpression = GeneExpression
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

  ##########
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
  #use_r("findDCD")
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

  sControl <- findSyntheticControl(GeneExpression=GeneExpression, FS=FS, PPIrank)
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

  if("patient"%in%colnames(GeneExpression)){
    CDrank <- rankByMut(genesIds = CD
                        ,project = "BRCA"
                        , patient_IDs = as.matrix(GeneExpression[,"patient"]))

  }else{
    CDrank <- rankByMut(genesIds = CD
                        ,project = "BRCA")
  }

  CDinfer <-  CDrank%>%
    mutate(Ensembl.ID = row.names(CDrank),.before=1)

  temp <- sapply(CausalImp[CD]
                 , function(x){
                   return(x[1,c("RelEffect","AbsEffect","p")])
                 },simplify = FALSE)

  temp <- do.call(what = rbind,temp)
  temp <- temp%>%
    rownames_to_column(var = "Ensembl.ID")
  CDinfer <- left_join(x = CDinfer, y = temp,by="Ensembl.ID")

  CDinfer <- CDinfer%>%
    arrange(desc(rank),desc(RelEffect))

  data(CGC.driverNames)
  nCD <- length(CD)
  nCDconf <- length(intersect(CGC.driverNames$Ensembl.ID
                              ,CD))
  nFSconf <- length(intersect(CGC.driverNames$Ensembl.ID
                              ,FS))
  phyper <- phyper(nCDconf-1, nFSconf
                   , length(FS)-nFSconf
                   , nCD, lower.tail = F, log.p = FALSE)

  summary <- matrix(
    c("Covariant(s) ptime: " ,colnames(pathCovariate)
      ,"Top PPI" ,PPItop
      ,"Find Event", paste0(findEvent,"(event at:",event$eventAt,")")
      ,"nFS" ,length(FS)
      ,"nFSconf", nFSconf
      ,"inferred" ,nCD
      ,"Confirmed" ,nCDconf
      ,"phyper", phyper
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
