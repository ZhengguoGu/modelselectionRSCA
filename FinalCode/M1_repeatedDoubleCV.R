####### parameters #####################################################################
# DATA: concatenated data matrix
# Jk: number of variables (a vector, the elements of which correspond to the datablocks)
# R: number of components
# LassoSequence: A sequence of Lasso tuning parameters
# GLassoSequence: A sequence of Group Lasso tuning parameters
# n_rep: number of repetitions
# n_seg: number of segments
# N_cores: number of cores (for parallel computing)
#######################################################################################

M1_repeatedDoubleCV <- function(DATA, R, Jk, N_cores, LassoSequence,  GLassoSequence, n_rep, n_seg, NRSTARTS, nfolds, MaxIter){
  #library(RegularizedSCA)
  #library(foreach)
  #library(snow)
  #library(doSNOW)
  #library(doRNG)

  if(missing(N_cores)){
    N_cores = 1
  }
  
  if(missing(n_rep)){
    n_rep <- 10
  }
  
  if(missing(n_seg)){
    n_seg <- 3
  }
 
  
  if(missing(LassoSequence)){
    LassoSequence = seq(0.001, RegularizedSCA::maxLGlasso(DATA, Jk, R)$Lasso, length.out = 10)
  }
  
  if(missing(GLassoSequence)){
    GLassoSequence = seq(0.001, RegularizedSCA::maxLGlasso(DATA, Jk, R)$Glasso, length.out = 10)
  }
  
  if(missing(NRSTARTS)){
    NRSTARTS = 1
  }
  
  
  DATA <- data.matrix(DATA)
  nsub <- dim(DATA)[1]

 
################################ Repeated Double CV ###########################################


perc_test <- 1/n_seg

cl <- snow::makeCluster(N_cores)
doSNOW::registerDoSNOW(cl)

#note that set.seed() and %dorng% ensure that parallel computing generates reproducable results.
sim_result <- foreach::foreach(r = 1:n_rep, .combine='cbind') %dorng% {
    
    Flag_person <- 0  # this, together with the following while(min_person=0), is to ensure that
    # the number of subject in testset > number of components. !!! important
    randindex <- stats::runif(nsub, 0, 1)
    testset_indexsize <- array(NA, n_seg) #note that each element in the array will be compared to R.
    while(Flag_person == 0){
      for(i in 1:n_seg){
        #randindex <- runif(nsub, 0, 1)
        testset_index <- ((i - 1) * perc_test < randindex) & (randindex < i * perc_test)
        testset_indexsize[i] <- sum(testset_index)
      }
      if(sum(testset_indexsize<=R)>0){
        randindex <- stats::runif(nsub, 0, 1)  #The number of persons should be > than R. if not, resample.
      }else{
        Flag_person <- 1
      }
    }
  
    #e_hat <- array()
    OptimumLasso <- array()
    OptimumGLasso <- array()
    for(i in 1:n_seg){
    
      testset_index <- ((i - 1) * perc_test < randindex) & (randindex < i * perc_test)
      #testset <- DATA[testset_index, ]
      calibset <- DATA[!testset_index, ]
    
      results_innerloop <- RegularizedSCA::cv_sparseSCA(calibset, Jk, R, MaxIter, NRSTARTS, LassoSequence, GLassoSequence, nfolds, method = "component")
      OptimumLasso[i] <- results_innerloop$RecommendedLambda[1]
      OptimumGLasso[i] <- results_innerloop$RecommendedLambda[2]
      #estimatedP <- results_innerloop$P_hat
      #A <- t(estimatedP) %*% t(testset) 
      #SVD_DATA <- svd(A, R, R)
      #estimatedT <- SVD_DATA$v %*% t(SVD_DATA$u)
      #e_hat[i] <- sum((testset - estimatedT %*% t(estimatedP))^2)  COMMENTS: According to Filzmoser, Liebmann, & Varmuza (2009), no need to 
      #compute MSE. One may simply check the frequency table of the optimal lasso and glasso tuning parameter values. This is why I omit the MSE steps. 
    
    }
    
    final_sim <- cbind(OptimumLasso, OptimumGLasso)
  
    #final_sim[[3]] <- e_hat
    return(final_sim)
}

snow::stopCluster(cl)

#results <- list()
#results$Lasso <- OptimumLasso
#results$GLasso <- OptimumGLasso
#results$e_hat <- E_hat

return_tables <- list(Lasso = table(sim_result[, colnames(sim_result) == "OptimumLasso"]), 
                      GroupLasso = table(sim_result[, colnames(sim_result) == "OptimumGLasso"]))
return(return_tables)
}
