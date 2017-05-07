rdCV_RSCA <- function(DATA, Jk, R, LassoSequence, GLassoSequence, n_rep, n_seg, MaxIter, NRSTARTS, nfolds){
  
  DATA <- data.matrix(DATA)
  nsub <- dim(DATA)[1]
  
    
  E_hat <- list()
  OptimumLasso <- matrix(NA, n_rep, n_seg)
  OptimumGLasso <- matrix(NA, n_rep, n_seg)
  
  r = 1
  while(r <= n_rep){
    
    randindex <- runif(nsub, 0, 1)
    perc_test <- 1/n_seg

    e_hat <- list()
    
    for(i in 1:n_seg){
      testset_index <- ((i - 1) * perc_test < randindex) & (randindex < i * perc_test)
      testset <- DATA[testset_index, ]
      calibset <- DATA[!testset_index, ]
      
      results_innerloop <- cv_sparseSCA(calibset, Jk, R, MaxIter, NRSTARTS, LassoSequence, GLassoSequence, nfolds, method = "component")
      OptimumLasso[r, i] <- results_innerloop$OptimumLasso
      OptimumGLasso[r, i] <- results_innerloop$OptimumGLasso
      estimatedP <- results_innerloop$PoutFinal
      A <- t(estimatedP) %*% t(testset)    
      SVD_DATA <- svd(A, R, R)
      estimatedT <- SVD_DATA$v %*% t(SVD_DATA$u)
    
      e_hat[[i]] <- testset - estimatedT %*% t(estimatedP)  

    }
    E_hat[[r]] <- e_hat
    r <- r + 1
  }
    
  results <- list()
  results$Lasso <- OptimumLasso
  results$GLasso <- OptimumGLasso
  results$e_hat <- E_hat
  return(results)
    
  
  
}