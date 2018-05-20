#Bolasso with CV

# Input parameters:
# DATA: Concatenated data matrix (!!! Not standardized! the matrix will be standardized after each resampling)
# Jk: A vector. Each element of this vector is the number of columns of a data block
# N_boots: number of bootstrap relicates
# R: The number of components (R>=2)
# LassoSequence: A vector of lasso tuning parameter values in accending order
# GlassoSequence: A vector of Group Lasso tuning parameter values in accending order
# w: whether the data should be weighted by its number of variables; If yes, then w=TRUE. This is used because after each resampling, the data needs to be standardized by pre_process(), a function from the RegularizedSCA package. 
# NRSTARTS: Number of multistarts
BolassoCV <- function(DATA, Jk, R, LassoSequence, GlassoSequence, W, NRSTARTS){
  
  DATA <- data.matrix(DATA)
  
  n_persons <- nrow(DATA)
  
  person_index <- sample(1:n_persons, n_persons, replace = TRUE)
  Data_sample <- DATA[person_index, ]
  Data_sample <- RegularizedSCA::pre_process(Data_sample, weight = w)
  result <- cv_sparseSCA(Data_sample, Jk, R, LassoSequence, GLassoSequence)  
  T_target <- result$T_hat            #We fix the estimated T matrix from the first resampled data.
                                      #All the estimated T matrix are to be compared to this estimated T.       
                                      #(This is due to permutation freedom)
  P_indexset <- result$P_hat
  P_indexset[which(P_indexset!=0)] <- 1  #non-zero loadings are marked as 1
  
  #resampling
  i <- 1
  while(i < m){
    person_index <- sample(1:n_persons, n_persons, replace = TRUE)
    Data_sample <- DATA[person_index, ]
    Data_sample <- RegularizedSCA::pre_process(Data_sample, weight = w)
    
    result <- cv_sparseSCA(Data_sample, Jk, R, LassoSequence, GLassoSequence)
    
    T_result <- result$T_hat
    perm <- RegularizedSCA::TuckerCoef(T_target, T_result)$perm
    P_result <- result$P_hat[, perm]
    
    P_result[which(P_result!=0)] <- 1
    
    P_indexset <- P_indexset + P_result
    
  }
  
  # Reestimate P and T
  Pout3d <- list()
  Tout3d <- list()
  LOSS <- array()
  LOSSvec <- list()
  
  for (n in 1:NRSTARTS) {
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
  
  PoutBest <- Pout3d[[k]]
  ToutBest <- Tout3d[[k]]
  
  return_varselect <- list()
  return_varselect$Pmatrix <- PoutBest
  return_varselect$Tmatrix <- ToutBest
  
  return(return_varselect)
}