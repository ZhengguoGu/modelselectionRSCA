#Bolasso with CV

# Input parameters:
# DATA: Concatenated data matrix (!!! Not standardized! the matrix will be standardized after each resampling)
# Jk: A vector. Each element of this vector is the number of columns of a data block
# N_boots: number of bootstrap relicates
# R: The number of components (R>=2)
# LassoSequence: A vector of lasso tuning parameter values in accending order
# GlassoSequence: A vector of Group Lasso tuning parameter values in accending order


Bolasso_CV <- function(DATA, Jk, R, N_boots, LassoSequence, GLassoSequence){
  
  library(foreach)
  library(doSNOW)
  library(doRNG)
  library(RegularizedSCA)
  
  DATA <- data.matrix(Data_final)  #DATA should be pre-processed at this stage
  
  n_persons <- nrow(DATA)
  person_index <- sample(1:n_persons, n_persons, replace = TRUE)
  Data_sample <- DATA[person_index, ]
  result <- RegularizedSCA::cv_sparseSCA(Data_sample, Jk, R, MaxIter = 400, NRSTARTS = 2, LassoSequence, GLassoSequence, nfolds = 10, method = "component")  
  T_target <- result$T_hat                #We fix the estimated T matrix from the first resampled data.
  #All the estimated T matrix are to be compared to this estimated T.       
  #(This is due to permutation freedom)
  P_indexset <- result$P_hat
  P_indexset[which(P_indexset!=0)] <- 1  #non-zero loadings are marked as 1
  
  #resampling
  cl <- makeCluster(2)
  registerDoSNOW(cl)
  #note that set.seed() and %dorng% ensure that parallel computing generates reproducable results.
  sim_result <- foreach(i = 1:(N_boots-1), .combine='+') %dorng% {
    
    person_index <- sample(1:n_persons, n_persons, replace = TRUE)
    Data_sample <- DATA[person_index, ]
    result <- RegularizedSCA::cv_sparseSCA(Data_sample, Jk, R, MaxIter = 400, NRSTARTS = 2, LassoSequence, GLassoSequence, nfolds = 10, method = "component")
    T_result <- result$T_hat
    perm <- RegularizedSCA::TuckerCoef(T_target, T_result)$perm
    P_result <- result$P_hat[, perm]
    P_result[which(P_result!=0)] <- 1
    
    return(P_result)
    
  }
  stopCluster(cl)
  
  P_indexset <- sim_result + P_indexset
  
  P_indexset[which(P_indexset) != N_boots] <- 0  #Variables that have not been selected N_boots times are left to be zero
  
  # Reestimate P and T, with 20 starts
  Pout3d <- list()
  Tout3d <- list()
  LOSS <- array()
  LOSSvec <- list()
  
  for (n in 1:20) { 
    VarSelectResult <- StrucSCA_withIndex(DATA, Jk, R, P_indexset = P_indexset, MaxIter=400)
    Pout3d[[n]] <- VarSelectResult$Pmatrix
    Tout3d[[n]] <- VarSelectResult$Tmatrix
    LOSS[n] <- VarSelectResult$Loss
    LOSSvec[[n]] <- VarSelectResult$Lossvec
  }
  k <- which(LOSS == min(LOSS))
  if (length(k) > 1) {
    pos <- sample(1:length(k), 1)
    k <- k[pos]
  }
  
  PoutBest <- Pout3d[[k]]  #this is the final, estimated P
  ToutBest <- Tout3d[[k]]  #this is the final, estimated T
  
  return_list <- list(P_hat = PoutBest, T_hat = ToutBest)
  return(return_list)
 
}

  
  