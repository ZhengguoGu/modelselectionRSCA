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

set.seed(112)
library(RegularizedSCA)
library(foreach)
library(doSNOW)
library(doRNG)

N_cores <- 2

DATA <- data.matrix(Data_final)
nsub <- dim(DATA)[1]

LassoSequence = seq(.00001, RegularizedSCA::maxLGlasso(DATA, Jk, R)$Lasso, length.out = 5)
GLassoSequence = seq(.00001, RegularizedSCA::maxLGlasso(DATA, Jk, R)$Glasso, length.out = 5)

n_rep <- 2
n_seg <- 3
################################ Repeated Double CV ###########################################


perc_test <- 1/n_seg

cl <- makeCluster(N_cores)
registerDoSNOW(cl)

#note that set.seed() and %dorng% ensure that parallel computing generates reproducable results.
sim_result <- foreach(r = 1:n_rep, .combine='cbind') %dorng% {
    
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
        randindex <- runif(nsub, 0, 1)  #The number of persons should be > than R. if not, resample.
      }else{
        Flag_person <- 1
      }
    }
  
    #e_hat <- array()
    OptimumLasso <- array()
    OptimumGLasso <- array()
    for(i in 1:n_seg){
    
      testset_index <- ((i - 1) * perc_test < randindex) & (randindex < i * perc_test)
      testset <- DATA[testset_index, ]
      calibset <- DATA[!testset_index, ]
    
      results_innerloop <- RegularizedSCA::cv_sparseSCA(calibset, Jk, R, MaxIter = 400, NRSTARTS = 2, LassoSequence, GLassoSequence, nfolds = 10, method = "component")
      OptimumLasso[i] <- results_innerloop$RecommendedLambda[1]
      OptimumGLasso[i] <- results_innerloop$RecommendedLambda[2]
      estimatedP <- results_innerloop$P_hat
      A <- t(estimatedP) %*% t(testset) 
      
      SVD_DATA <- svd(A, R, R)
      estimatedT <- SVD_DATA$v %*% t(SVD_DATA$u)
    
      #e_hat[i] <- sum((testset - estimatedT %*% t(estimatedP))^2)  
    
    }
    
    final_sim <- list(Lasso = OptimumLasso, Glasso = OptimumGLasso)
  
    #final_sim[[3]] <- e_hat
    return(final_sim)
}

stopCluster(cl)

#results <- list()
#results$Lasso <- OptimumLasso
#results$GLasso <- OptimumGLasso
#results$e_hat <- E_hat

table(OptimumLasso)
table(OptimumGLasso)