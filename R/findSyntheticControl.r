findSyntheticControl <- function(GeneExpression, FS, PPIrank=NULL){
  GeneExpression <-dplyr::as_tibble(GeneExpression, rownames=NA)

  if(is.null(PPIrank))
    PPIrank <- CreatePPIRank(colnames(GeneExpression))

  PPIinDataset <- intersect(PPIrank[[1]]
                            , colnames(GeneExpression))
  #isEnsembl <- which(str_starts(colnames(GeneExpression), "ENSG0", negate = FALSE))
  GeneExpression.comp <- GeneExpression%>%
    #dplyr::select(all_of(isEnsembl))%>%
    dplyr::select(!all_of(PPIinDataset))%>%
    transmute_all(as.numeric)

  DS <- GeneExpression%>%
    dplyr::select(all_of(FS))%>%
    transmute_all(as.numeric)


  system.time(Pearson <- cor(x = as.matrix(GeneExpression.comp)
                             , y = as.matrix(DS)))
  sControl <- as.matrix(apply(Pearson, MARGIN = 2, which.max))
  sControl <- cbind(row.names(sControl),row.names(Pearson)[sControl])
  colnames(sControl) <- c("gene","scontrol")

  return(sControl)
}
