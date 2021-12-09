#' @title findCovariate
#'
#' @description For each feature in \code{FS}, \code{findCovariate}
#'  finds the non-PPI gene with the largest Pearson correlation with the feature.
#'  For a proper functioning, the data of all features in \code{FS} must be included as columns in
#'  \code{GeneExpression}.
#'
#'
#' @usage function(GeneExpression, FS, PPIrank=NULL)
#'
#' @inheritParams findDCD
#' @param FS A \code{character} vector containing the names of the features
#'  for which a covariate is to be found.
#' @param PPIrank (optional) A \code{dataframe} obtained from
#'  \code{\link{CreatePPIRank}}
#'
#'
#' @author Andres Mauricio Cifuentes_Bernal, Vu VH Pham, Xiaomei Li, Lin Liu, JiuyongLi and Thuc Duy Le
#' @export
#' @seealso \link[DynamicCancerDriver]{findDCD}
#'
#' @return A \code{dataframe} with the following two variables:
#'   \enumerate{
#'           \item{\emph{Feature: }}{The vector \code{FS} of features. \cr
#'           The name of this variable can be 1 of \emph{Ensembl.ID, HGNC.ID, NCBI.ID}
#'           or \emph{HGNC.symbol}}.
#'           \item{\code{scontrol: }}{For each \emph{Feature}, the name of the non-PPI gene with the
#'           largest Pearson correlation}
#'           }
#'
#' @examples \dontrun{
#'    data("GSE75688_TPM_tumor", package = "DynamicCancerDriver")
#'    FS <- colnames(GSE75688_TPM_tumor)[1:100]
#'
#'    sControl <- findCovariate(GeneExpression = GSE75688_TPM_tumor[,1:500]
#'    , FS = FS)
#'    }
#'
#'
#' @references

findCovariate <- function(GeneExpression, FS, PPIrank=NULL){
  GeneExpression <-dplyr::as_tibble(GeneExpression, rownames=NA)

  if(is.null(PPIrank))
    PPIrank <- CreatePPIRank(colnames(GeneExpression))

  PPIinDataset <- intersect(PPIrank[[1]]
                            , colnames(GeneExpression))

  GeneExpression.comp <- GeneExpression%>%
    dplyr::select(!all_of(PPIinDataset))%>%
    transmute_all(as.numeric)

  DS <- GeneExpression%>%
    dplyr::select(all_of(FS))%>%
    transmute_all(as.numeric)


  system.time(Pearson <- cor(x = as.matrix(GeneExpression.comp)
                             , y = as.matrix(DS)))

  # Change to 0 values of self correlation
  temp <- row.names(Pearson)%in%colnames(Pearson)
  index <- row.names(Pearson)[temp]
  for (i in index) {
    Pearson[i,i] <- 0
  }

  sControl <- as.matrix(apply(Pearson, MARGIN = 2, which.max))
  sControl <- cbind(row.names(sControl),row.names(Pearson)[sControl])
  colnames(sControl) <- c("Feature","scontrol")

  return(sControl)
}
