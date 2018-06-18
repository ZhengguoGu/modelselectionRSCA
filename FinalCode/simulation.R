library(RegularizedSCA)
library(foreach)
library(snow)
library(doSNOW)
library(doRNG)

##############################################################################################
# StrucSCA_withIndex() estimates T and P, given the pre-defined structure of P 
##############################################################################################
StrucSCA_withIndex <- function (DATA, Jk, R, P_indexset, MaxIter) {
  DATA <- data.matrix(DATA)
  
  I_Data <- dim(DATA)[1]
  sumJk <- dim(DATA)[2]
  eps <- 10^(-12)
  if (missing(MaxIter)) {
    MaxIter <- 400
  }
  P <- matrix(stats::rnorm(sumJk * R), nrow = sumJk, ncol = R)
  P[P_indexset == 0] <- 0
  Pt <- t(P)
  
  sqP <- P^2
  residual <- sum(DATA^2)
  Lossc <- residual 
  
  conv <- 0
  iter <- 1
  Lossvec <- array()
  while (conv == 0) {
    
    SVD_DATA <- svd(DATA, R, R)
    Tmat <- SVD_DATA$u
    
    residual <- sum((DATA - Tmat %*% Pt)^2)
    Lossu <- residual 
    
    P <- t(DATA) %*% Tmat
    P[P_indexset == 0] <- 0
    Pt <- t(P)
    
    sqP <- P^2
    residual <- sum((DATA - Tmat %*% Pt)^2)
    Lossu2 <- residual
    if (abs(Lossc - Lossu) < 10^(-9)) {
      Loss <- Lossu
      residual <- residual
      P[abs(P) <= 2 * eps] <- 0
      conv <- 1
    }
    else if (iter > MaxIter) {
      Loss <- Lossu
      residual <- residual
      P[abs(P) <= 2 * eps] <- 0
      conv <- 1
    }
    Lossvec[iter] <- Lossu
    iter <- iter + 1
    Lossc <- Lossu2
  }
  return_varselect <- list()
  return_varselect$Pmatrix <- P
  return_varselect$Tmatrix <- Tmat
  return_varselect$Loss <- Loss
  return_varselect$Lossvec <- Lossvec
  return(return_varselect)
}

####################################################################################################
# Calculate the number of variables correctly selected
###################################################################################################

num_correct <- function (TargetP, EstimatedP){
  
  total_vnumber <- dim(TargetP)[1] * dim(TargetP)[2]
  
  TargetP[which(TargetP != 0)] <- 1
  sum_select <- sum(TargetP)
  sum_zero <- total_vnumber - sum_select
  
  EstimatedP[which(EstimatedP != 0)] <- 1
  
  total_correct <- sum(TargetP == EstimatedP) # this is the total number of variables correctedly selected and zeros correctly retained
  
  prop_correct <- total_correct/total_vnumber
  
  return(prop_correct)
  
}


####################################################################################################################################
####################################################################################################################################

##############################################################################################################
####
#### Simulations: manually change the parameters: PropNoise and Perc0
####
##############################################################################################################



I <- 28
J1 <- 144
J2 <- 44
Jk <- c(J1, J2)
R <- 3
NRSTARTS <- 5
n_rep = 20
n_seg = 3
N_boots = 20

### 1. benchmark CV
set.seed(1)
n_dataset <- 1
N_dataset = 20
RESULT_BenchmarCV <- matrix(NA, N_dataset, 2)
ESTIMATED_P <- list()
ESTIMATED_T <- list()
while(n_dataset <= N_dataset){
  
  filename <- paste("Data_", n_dataset, ".RData", sep = "")
  load(filename)
  
  Lassosequence <- seq(0.0000001, RegularizedSCA::maxLGlasso(my_data_list$data, Jk, R)$Lasso, length.out = 50)
  GLassosequence <- seq(0.0000001, RegularizedSCA::maxLGlasso(my_data_list$data, Jk, R)$Glasso, length.out = 50)
 
  result_sim1_BM <- RegularizedSCA::cv_sparseSCA(my_data_list$data, Jk, R, MaxIter = 400, NRSTARTS, Lassosequence, GLassosequence, nfolds = 7, method = "component") 
  tuckerresult <- RegularizedSCA::TuckerCoef(my_data_list$T_mat, result_sim1_BM$T_hat)
  RESULT_BenchmarCV[n_dataset, 1] <- tuckerresult$tucker_value
  RESULT_BenchmarCV[n_dataset, 2] <- num_correct(my_data_list$P_mat, result_sim1_BM$P_hat[, tuckerresult$perm])  
  
  ESTIMATED_P[[n_dataset]] <- result_sim1_BM$P_hat
  ESTIMATED_T[[n_dataset]] <- T_hat = result_sim1_BM$T_hat
  n_dataset <- n_dataset + 1
}


save(RESULT_BenchmarCV, ESTIMATED_P, ESTIMATED_T, file = "BenchmarkCV.RData")





### 2. repeated Double CV
n_dataset <- 1
N_dataset = 20
RESULT_rdCV <- matrix(NA, N_dataset, 2)
ESTIMATED_PrdCV <- list()
ESTIMATED_TrdCV <- list()

set.seed(1)
while(n_dataset <= N_dataset){
  
  filename <- paste("Data_", n_dataset, ".RData", sep = "")
  load(filename)
  
  Lassosequence <- seq(0.0000001, RegularizedSCA::maxLGlasso(my_data_list$data, Jk, R)$Lasso, length.out = 50)
  GLassosequence <- seq(0.0000001, RegularizedSCA::maxLGlasso(my_data_list$data, Jk, R)$Glasso, length.out = 50)
  
  result_sim1_RDCV <- M1_repeatedDoubleCV(my_data_list$data,  R, Jk, N_cores = 20, LassoSequence = Lassosequence, GLassoSequence = GLassosequence, n_rep , n_seg, NRSTARTS)
  
  temp_lasso <- as.data.frame(result_sim1_RDCV$Lasso)
  temp_lasso$Var1 <- sort(as.numeric(levels(temp_lasso$Var1)))
  LASSO <- max(temp_lasso[temp_lasso[,2] == max(temp_lasso[,2]),1])  #the first max ensures that the largest Lasso value is chosen, in case more than one lasso value is recommended by M1_repeatedDoubleCV
  temp_glasso <- as.data.frame(result_sim1_RDCV$GroupLasso)
  temp_glasso$Var1 <- sort(as.numeric(levels(temp_glasso$Var1)))
  GLASSO <- max(temp_glasso[temp_glasso[,2] == max(temp_glasso[,2]),1])
  
  
  final_RDCV <- RegularizedSCA::sparseSCA(my_data_list$data, Jk, R, LASSO = LASSO, GROUPLASSO = GLASSO, MaxIter = 400,
                                          NRSTARTS = 20, method = "component")
  
  tuckerresult_RDCV <- RegularizedSCA::TuckerCoef(my_data_list$T_mat, final_RDCV$Tmatrix)
  RESULT_rdCV[n_dataset, 1] <- tuckerresult_RDCV$tucker_value 
  RESULT_rdCV[n_dataset, 2] <- num_correct(my_data_list$P_mat, final_RDCV$Pmatrix[, tuckerresult_RDCV$perm])  
  
  ESTIMATED_PrdCV[[n_dataset]] <- final_RDCV$Pmatrix
  ESTIMATED_TrdCV[[n_dataset]] <- final_RDCV$Tmatrix
  n_dataset <- n_dataset + 1
  
}

save(RESULT_rdCV, ESTIMATED_PrdCV, ESTIMATED_TrdCV, file = "RepeatedDCV.RData")



### 3. BIC and IS
n_dataset <- 1
N_dataset = 20
RESULT_BIC <- matrix(NA, N_dataset, 2)
RESULT_IS <- matrix(NA, N_dataset, 2)
ESTIMATED_Pbic <- list()
ESTIMATED_Tbic <- list()
ESTIMATED_PIS <- list()
ESTIMATED_TIS <- list()

set.seed(1)
while(n_dataset <= N_dataset){
  filename <- paste("Data_", n_dataset, ".RData", sep = "")
  load(filename)
  
  Lassosequence <- seq(0.0000001, RegularizedSCA::maxLGlasso(my_data_list$data, Jk, R)$Lasso, length.out = 50)
  GLassosequence <- seq(0.0000001, RegularizedSCA::maxLGlasso(my_data_list$data, Jk, R)$Glasso, length.out = 50)
  
  result_sim_BICIS <- M2_BIC_IS(my_data_list$data, Jk, R, LassoSequence = Lassosequence, GLassoSequence = GLassosequence, NRSTARTS)
  

  Croux_index <- which(result_sim_BICIS$Croux == min(result_sim_BICIS$Croux), arr.ind = T)
  Lasso_croux <- max(Lassosequence[Croux_index[1]])  #max() is used in case multiple lasso values are chosen. 
  GLasso_croux <- max(GLassosequence[Croux_index[2]]) 
  final_croux <- RegularizedSCA::sparseSCA(my_data_list$data, Jk, R, LASSO = Lasso_croux, GROUPLASSO = GLasso_croux, MaxIter = 400, NRSTARTS = 20, method = "component")
  ESTIMATED_Pbic[[n_dataset]] <- final_croux$Pmatrix
  ESTIMATED_Tbic[[n_dataset]] <- final_croux$Tmatrix
  tuckerresult_croux <- RegularizedSCA::TuckerCoef(my_data_list$T_mat, final_croux$Tmatrix)    
  RESULT_BIC[n_dataset, 1] <- tuckerresult_croux$tucker_value 
  RESULT_BIC[n_dataset, 2] <- num_correct(my_data_list$P_mat, final_croux$Pmatrix[, tuckerresult_croux$perm])  
  
  
  IS_index <- which(result_sim_BICIS$IS == max(result_sim_BICIS$IS), arr.ind = T)
  Lasso_IS <- max(Lassosequence[IS_index[1]])
  Glasso_IS <- max(GLassosequence[IS_index[2]])
  final_IS <- RegularizedSCA::sparseSCA(my_data_list$data, Jk, R, LASSO = Lasso_IS, GROUPLASSO = Glasso_IS, MaxIter = 400, NRSTARTS = 20, method = "component")
  ESTIMATED_PIS[[n_dataset]] <- final_IS$Pmatrix
  ESTIMATED_TIS[[n_dataset]] <- final_IS$Tmatrix
  tuckerresult_IS <- RegularizedSCA::TuckerCoef(my_data_list$T_mat, final_IS$Tmatrix)    
  RESULT_IS[n_dataset, 1] <- tuckerresult_IS$tucker_value
  RESULT_IS[n_dataset, 2] <- num_correct(my_data_list$P_mat, final_IS$Pmatrix[, tuckerresult_IS$perm])  
  
  n_dataset <- n_dataset + 1
}

save(RESULT_BIC, RESULT_IS, ESTIMATED_Pbic, ESTIMATED_Tbic, ESTIMATED_PIS, ESTIMATED_TIS,  file = "BIC_IS.RData")



### 4. Bolasso
n_dataset <- 1
N_dataset = 20
RESULT_BoLasso <- matrix(NA, N_dataset, 2)
ESTIMATED_Pbolasso <- list()
ESTIMATED_Tbolasso <- list()

set.seed(1)
while(n_dataset <= N_dataset){
  
  filename <- paste("Data_", n_dataset, ".RData", sep = "")
  load(filename)
  
  Lassosequence <- seq(0.0000001, RegularizedSCA::maxLGlasso(my_data_list$data, Jk, R)$Lasso, length.out = 50)
  GLassosequence <- seq(0.0000001, RegularizedSCA::maxLGlasso(my_data_list$data, Jk, R)$Glasso, length.out = 50)
  
  result_sim1_Bolasso <- Bolasso_CV(my_data_list$data, Jk, R, N_boots, LassoSequence = Lassosequence, GLassoSequence = GLassosequence, NRSTARTS)

  tuckerresult_Bolasso <- RegularizedSCA::TuckerCoef(my_data_list$T_mat, result_sim1_Bolasso$T_hat)    
  RESULT_BoLasso[n_dataset, 1] <- tuckerresult_Bolasso$tucker_value 
  RESULT_BoLasso[n_dataset, 2] <- num_correct(my_data_list$P_mat, result_sim1_Bolasso$P_hat[, tuckerresult_Bolasso$perm])
  ESTIMATED_Pbolasso[[n_dataset]] <- result_sim1_Bolasso$P_hat
  ESTIMATED_Tbolasso[[n_dataset]] <- result_sim1_Bolasso$T_hat
  
  n_dataset <- n_dataset + 1
}
