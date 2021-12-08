#' @title findCovariate
#'
#' @description For each feature in \code{FS} \code{findCovariate}
#' finds the non-PPI gene with the largest pearson correlation with the feature.
#' The data of all features in \code{FS}, need to be included in the
#' \code{GeneExpression} parameter.
#'
#'
#' @usage function(GeneExpression, FS, PPIrank=NULL)\cr
#'
#' @inheritParams findDCD
#' @param FS A \code{character} vector containing the names of the features
#'  for which a covariate is to be found.
#' @param PPIrank (optional) A \code{dataframe} obtained from \code{\link{CreatePPIRank}}
#'
#'
#' @author Andres Mauricio Cifuentes_Bernal, Vu VH Pham, Xiaomei Li, Lin Liu, JiuyongLi and Thuc Duy Le
#' @export
#' @seealso \link[DynamicCancerDriver]{findCovariate},
#' \link[DynamicCancerDriver]{parallelCI}
#'
#' @return A \code{dataframe} with the following two variables:
#'   \enumerate{
#'           \item{\code{Feature:}}{The vector \code{FS} of features}
#'           \item{\code{scontrol:}}{The name of the non-PPI gene with the
#'           largest Pearson correlation with the respective feature}
#'           }
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
