rdCV_RSCA <- function(DATA, Jk, R, LassoSequence, GLassoSequence, n_rep, n_seg, MaxIter, NRSTARTS, nfolds, ncores){
  
  DATA <- data.matrix(DATA)
  nsub <- dim(DATA)[1]
  
  #1 repetition
  randindex <- runif(nsub, 0, 1)
  perc_test <- 1/n_seg
  
  #cl <- makeCluster(ncores) 
  #registerDoSNOW(cl)
  
  e_hat <- list()
  
  #result <- foreach(i = 1:n_seg) %dorng% {
  for(i in 1:n_seg){
    testset_index <- ((i - 1) * perc_test < randindex) & (randindex < i * perc_test)
    testset <- DATA[testset_index, ]
    calibset <- DATA[!testset_index, ]
    
    results_innerloop <- cv_sparseSCA(calibset, Jk, R, MaxIter, NRSTARTS, LassoSequence, GLassoSequence, nfolds, method = "component")
        
    estimatedP <- results_innerloop$PoutFinal
    A <- t(estimatedP) %*% t(testset)    
    SVD_DATA <- svd(A, R, R)
    estimatedT <- SVD_DATA$v %*% t(SVD_DATA$u)
    
    #e_hat <- testset - estimatedT %*% t(estimatedP)  
    e_hat[[i]] <- testset - estimatedT %*% t(estimatedP)  
    #return(e_hat)
  }
  #stopCluster(cl)
  #return(result)
  
  return(e_hat)
}