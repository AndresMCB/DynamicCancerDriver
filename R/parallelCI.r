parallelCI <- function(GeneExpression,sControl,z,eventAt
                       , CIniter=200, chunk_size = 50
                       , returnModel = F){
  if(!require(parallel))
    install.packages("parallel")

  library(parallel)


  DS.order <- GeneExpression%>%
    dplyr::select(all_of(sControl[,1]),all_of(sControl[,2]))%>%
    mutate(z=z, .before = 1)%>%
    arrange(z)%>%
    transmute_all(as.numeric)

  DS.order <- as.matrix(DS.order)
  pre.period <- DS.order[c(1,eventAt-1), 1, drop=T]
  post.period <- DS.order[c(eventAt,nrow(GeneExpression)), 1, drop=T]



  CausalImp<-vector(mode = "list",length = nrow(sControl))
  names(CausalImp) <- sControl[,"gene"]

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
