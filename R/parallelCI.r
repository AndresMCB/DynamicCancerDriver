#' @title parallelCI
#'
#' @description parallelCI uses \code{parallel} package for implementing
#' a parallelised calculation of the \code{CausalImpact} on gene expression.
#' It is assumed that the \code{GeneExpression} matrix provided contains
#' pseudotime ordered gene expression. \code{z} is the pseudotime score used for
#' ordering the gene expression, and \code{eventAt} indicates the sample at
#' the most significant change occurs.
#'
#'
#'
#' @usage function(GeneExpression,sControl,z,eventAt
#' , chunk_size = 50, returnModel = F)\cr
#' @inheritParams findDCD
#'
#' @param sControl A 2 column matrix containing a gene ID (1st column) and 1 non-PPI
#'  gene to be uses as covariate for the \code{CausalImpact} model. For a correct
#'  functioning, the pseudotime ordered data of all element in \code{sControl} need
#'  to be included as columns in \code{GeneExpression} matrix.
#' @param z A numeric vector containing the pseudotime score used for ordering
#' the samples in \code{GeneExpression}. For a correct functioning, \code{z} needs to
#' follow ascending order and this order must agree with the order of the samples (rows)
#' of the \code{GeneExpression} matrix.
#' @param eventAt An integer with the index (in pseudotime order) of the sample where
#' the most significant change is inferred to happen.
#' @param chunk_size An integer indicating the number of genes to be passed to
#' each worker during the parallel calculation. (50 by default)
#' @param returnModel A boolean. If \code{TRUE}, the full causal impact model
#' (as calculated by \code{CausalImpact} package) is returned. if \code{FALSE}
#' (default), only the main parameters of the \code{CausalImpact} are returned,
#'
#'
#' @author Andres Mauricio Cifuentes_Bernal, Vu VH Pham, Xiaomei Li, Lin Liu, JiuyongLi and Thuc Duy Le
#' @export
#' @seealso \link[DynamicCancerDriver]{findCovariate},
#' \link[DynamicCancerDriver]{parallelCI}
#'
#' @return A \code{list} where each element is the full \code{CausalImpact} model
#' (if \code{returnModel = TRUE}) or the simplified \code{CausalImpact} model
#' (if \code{returnModel = FALSE}) of one gene in the 1st column of \code{sControl}.
#'
#'
#'
#' @examples \dontrun{
#'    data("GSE75688_TPM_tumor", package = "DynamicCancerDriver")
#'    FS <- colnames(GSE75688_TPM_tumor)[1:100]
#'    GE <- GSE75688_TPM_tumor[,1:500]
#'    sControl <- findCovariate(GeneExpression = GSE75688_TPM_tumor[,1:500]
#'    , FS = FS)
#'
#'    #toy example, using "VIM" as path covariate
#'    z <- findZ(GeneExpression = GE
#'              , FS, pathCovariate = GSE75688_TPM_tumor[,"ENSG00000026025"]
#'              , elbo_tol = 1e-3)
#'    GE <- GE[order(z),,drop=F]
#'    z <- z[order(z),1, drop=F]
#'
#'    parCI <-parallelCI(GE,sControl,z,eventAt=7)
#'
#' }
#'
#' @references
#'
#'


parallelCI <- function(GeneExpression,sControl,z,eventAt
                       , chunk_size = 50
                       , returnModel = F){
  if(!require(parallel))
    install.packages("parallel")
  if(!require(CausalImpact))
    install.packages("CausalImpact")

  library(parallel)
  GeneExpression <-dplyr::as_tibble(GeneExpression, rownames=NA)  #to keep rownames

  CIniter <- 200 # parameter for CausalImpact model calculation

  DS.order <- GeneExpression%>%
    dplyr::select(all_of(sControl[,1]),all_of(sControl[,2]))%>%
    mutate(z=z, .before = 1)%>%
    arrange(z)%>%
    transmute_all(as.numeric)

  DS.order <- as.matrix(DS.order)
  pre.period <- DS.order[c(1,eventAt-1), 1, drop=T]
  post.period <- DS.order[c(eventAt,nrow(GeneExpression)), 1, drop=T]



  CausalImp<-vector(mode = "list",length = nrow(sControl))
  names(CausalImp) <- sControl[,1]

  copies_of_r <- detectCores(logical = F)-1
  n<-nrow(sControl)


  chunk <- split(1:n, ceiling(1:n/chunk_size))
  cl <- makeCluster(copies_of_r)
  clusterExport(cl, c("DS.order","CausalImpact", "zoo", "sControl",
                      "returnModel")
                ,envir = environment())

  for (i in 1:length(chunk)) {
    message(paste0("------ chunk ",toString(i), " of "
                   , toString(length(chunk))
                   ," ------"))
    CausalImp[chunk[[i]]] <-
      parSapply(cl, chunk[[i]]
                , function(a){

                  set.seed(1)
                  res <- CausalImpact(data = zoo(DS.order[,sControl[a,]]
                                                 ,order.by = DS.order[,1])
                                      , pre.period = pre.period
                                      , post.period = post.period
                                      , model.args = list(niter = CIniter) )
                  if(!returnModel)
                    res <- res[["summary"]]
                  return(res)
                }, simplify = F)
    gc()
  }# end for
  stopCluster(cl)

  failed <- which(sapply(CausalImp, is.null))
  if(length(failed)>0)
    CausalImp <- CausalImp[-failed]

  return(CausalImp)
}
