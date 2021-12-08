#' @title CreatePPIRank
#'
#' @description CreatePPIRank
#'
#' @usage function(geneIDs = NULL, PPImatrix = NULL)\cr
#'
#' @param geneIDs
#' @param PPImatrix
#'
#'
#' @author Andres Mauricio Cifuentes_Bernal, Vu VH Pham, Xiaomei Li, Lin Liu, JiuyongLi and Thuc Duy Le
#' @export
#' @seealso \link[DynamicCancerDriver]{findCovariate},
#' \link[DynamicCancerDriver]{parallelCI}
#'
#' @return A \code{dataframe} consisting of the following elements:
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
#'    findCovariate(GeneExpression = GSE75688_TPM_tumor[,1:500]
#'                 ,FS =colnames(GSE75688_TPM_tumor[,1:100]))
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



CreatePPIRank <- function(geneIDs = NULL, PPImatrix = NULL){

  if(!require(AMCBGeneUtils))
    devtools::install_github(repo = "AndresMCB/AMCBGeneUtils")

  library(AMCBGeneUtils)
  library(tidyverse)

  IdSource <- GeneIdSource(geneIDs)

  if(is.null(PPImatrix)){
    data("PPI")

    PPI <- PPI %>%
      transmute(input = changeGeneId(.[[2]], from = "NCBI.ID", to = IdSource)[[2]]
                , output = changeGeneId(.[[4]], from = "NCBI.ID", to = IdSource)[[2]])
    col_input <- 1
    col_output <- 2
  }else
  {
    aux <- dim(PPImatrix)
    if(dim(aux)==2 && aux[2]>=2){
      if(aux[2]>2){
        message("PPI matrix has more than 2 columns.")
        message("Using the 1st column with the word \"input\" and 1st column with the word +
                \"output.\"")
      }
      pattern <- regex("input", ignore_case = T)
      col_input <- which(str_detect(colnames(PPImatrix)
                                    ,pattern =pattern))[1]

      if(is.na(col_input))
        error("Not column name including \"input\" found.")

      pattern <- regex("output", ignore_case = T)
      col_output <- which(str_detect(colnames(PPImatrix)
                                     ,pattern =pattern))[1]
      if(is.na(col_output))
        error("Not column name including \"output\" found.")
    }else
    {
      error("PPImatrix has not the right dimensions")
    }

    IdSource <- GeneIdSource(PPI[,col_input])
    if(IdSource!=GeneIdSource(PPI[,col_output]))
      error("input and ouput nodes in the PPI are not +
          in the same nomenclature.")

    PPI <- as.data.frame(PPImatrix[,c(col_input, col_output),drop=F])

    if(IdSource=="HGNC.symbol"){
      PPI[,1]<-updateGeneSymbol(PPI[,1])[["HGNC.symbol"]]
      PPI[,2]<-updateGeneSymbol(PPI[,2])[["HGNC.symbol"]]
    }
  }

  geneIDs <- changeGeneId(geneIDs, from = IdSource)

  # data-specific PPI
  PPI_net <- PPI%>%
    filter(.[[col_input]] %in% geneIDs[[IdSource]])%>%
    filter(.[[col_output]] %in% geneIDs[[IdSource]])


  PPIrank <- unique(c(PPI_net[[col_input]]
                      ,PPI_net[[col_output]]))
  PPIrank <- as.data.frame(na.omit(PPIrank))
  colnames(PPIrank) <- IdSource

  In <- PPI_net%>%
    count(.[[col_input]],name = "n.in")

  colnames(In)[1] <- IdSource

  Out <-PPI_net%>%
    count(.[[col_output]],name = "n.out")
  colnames(Out)[1] <- IdSource

  PPIrank <- left_join(x = PPIrank,y = In, by = IdSource)
  PPIrank <- left_join(x = PPIrank,y = Out, by = IdSource)
  PPIrank[is.na(PPIrank)] <- 0
  PPIrank <- PPIrank%>%
    mutate(total=n.in+n.out)

  return(PPIrank)

}
