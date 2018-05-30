MS1_repDCV <- function(DATA, Jk, R, LassoSequence, GLassoSequence, n_rep, n_seg, MaxIter, NRSTARTS, nfolds){
  
####### parameters ###############
# DATA: concatenated data matrix
# Jk: number of variables (a vector, the elements of which correspond to the datablocks)
# R: number of components
# LassoSequence: A sequence of Lasso tuning parameters
# GLassoSequence: A sequence of Group Lasso tuning parameters
# n_rep: number of repetitions
# n_seg: number of segments
# MaxIter: maximum number of iterations
# NRSTARTS: number of multi-starts
# nfolds: number of folds for CV
  
  DATA <- data.matrix(DATA)
  nsub <- dim(DATA)[1]
  
  E_hat <- list()
  OptimumLasso <- matrix(NA, n_rep, n_seg)
  OptimumGLasso <- matrix(NA, n_rep, n_seg)
  
  perc_test <- 1/n_seg
  
  r = 1
  while(r <= n_rep){
    
    cat(sprintf("\n Repetition: %s", r))
    
    Flag_person <- 0  # this, together with the following while(min_person=0), is to ensure that
                     # the number of subject in testset > number of components. !!! important
    randindex <- runif(nsub, 0, 1)
    testset_indexsize <- array(NA, n_seg) #note that each element in the array will be compared to R.
    while(Flag_person == 0){
      for(i in 1:n_seg){
        #randindex <- runif(nsub, 0, 1)
        testset_index <- ((i - 1) * perc_test < randindex) & (randindex < i * perc_test)
        testset_indexsize[i] <- sum(testset_index)
      }
      if(sum(testset_indexsize<=R)>0){
          print("The algorithm is now resampling the persons (without replacement).")
          randindex <- runif(nsub, 0, 1)  #The number of persons should be > than R. if not, resample.
      }else{
          Flag_person <- 1
      }
    }
    
    e_hat <- array()
    
    for(i in 1:n_seg){
      
      cat(sprintf("\n Outer loop: %s", i))
      
      testset_index <- ((i - 1) * perc_test < randindex) & (randindex < i * perc_test)
      testset <- DATA[testset_index, ]
      calibset <- DATA[!testset_index, ]
      
      results_innerloop <- RegularizedSCA::cv_sparseSCA(calibset, Jk, R, MaxIter, NRSTARTS, LassoSequence, GLassoSequence, nfolds, method = "component")
      OptimumLasso[r, i] <- results_innerloop$RecommendedLambda[1]
      OptimumGLasso[r, i] <- results_innerloop$RecommendedLambda[2]
      estimatedP <- results_innerloop$P_hat
      A <- t(estimatedP) %*% t(testset) 
      
      SVD_DATA <- svd(A, R, R)
      estimatedT <- SVD_DATA$v %*% t(SVD_DATA$u)
      
      e_hat[i] <- sum((testset - estimatedT %*% t(estimatedP))^2)  
      
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