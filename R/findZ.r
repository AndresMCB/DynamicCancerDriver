
#' @title findZ
#'
#' @description \code{findZ} calculates a pseudotime score for each sample
#' by using a \code{pathCovariate} that reasonably encodes the progression
#' of one core biological process during cancer progression (in the sense described by
#' \href{https://www.nature.com/articles/s41467-018-04696-6}{Campbell 2018}).\cr
#' Pseudotime score calculation relies in the procedures implemented in the
#' \href{https://www.bioconductor.org/packages/release/bioc/html/phenopath.html}{phenopath} package.
#'
#' @usage function(GeneExpression, FS, pathCovariate, elbo_tol = 1e-3)\cr
#'
#' @inheritParams findDCD
#' @param FS A \code{character} vector containing the names of the features
#'  to be used for the calculation of the pseudotime score.
#' @param pathCovariate A \code{named vector} containing the data of a path
#' covariate.
#'
#' @author Andres Mauricio Cifuentes_Bernal, Vu VH Pham, Xiaomei Li, Lin Liu, JiuyongLi and Thuc Duy Le
#' @export
#' @seealso \link[DynamicCancerDriver]{findDCD}
#'
#' @return A \code{dataframe} with the following two variables:
#'   \enumerate{
#'           \item{\code{Feature:}}{The vector \code{FS} of features}
#'           \item{\code{scontrol:}}{For each \code{Feature}, the name of the
#'           non-PPI gene with the largest Pearson correlation.}
#'           }
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
#'              }
#'
#'
#' @references

#---- function to find pseudotime score (based on Phenopath)
findZ <- function(GeneExpression, FS, pathCovariate, elbo_tol = 1e-3){
  GeneExpression <-dplyr::as_tibble(GeneExpression, rownames=NA)  #to keep rownames
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
