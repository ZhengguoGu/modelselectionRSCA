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


set.seed(1)
I <- 28
J1 <- 144
J2 <- 44
Jk <- c(J1, J2)
R <- 5
NRSTARTS <- 5
n_rep = 20
n_seg = 3
N_boots = 20

### 1. benchmark CV

n_dataset <- 1
N_dataset = 20
RESULT_BenchmarCV <- matrix(N_dataset, 2)
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
  
  ESTIMATED_P[[n_dataset]] <- list(result_sim1_BM$P_hat)
  ESTIMATED_T[[n_dataset]] <- list(T_hat = result_sim1_BM$T_hat)
  n_dataset <- n_dataset + 1
}










### 2. repeated Double CV
set.seed(1)
result_sim1_RDCV <- M1_repeatedDoubleCV(my_data,  R, Jk, N_cores = 12, LassoSequence = Lassosequence, GLassoSequence = GLassosequence, n_rep , n_seg, NRSTARTS)

result_sim1_RDCV$Lasso  # 0.0747062083396429  the most sparse one was chozen
result_sim1_RDCV$GroupLasso  #0.0640494874513447

final_RDCV <- RegularizedSCA::sparseSCA(my_data, Jk, R, LASSO = 0.0747062083396429, GROUPLASSO = 0.0640494874513447, MaxIter = 400,
                                        NRSTARTS = 20, method = "component")

final_RDCV$Pmatrix
tuckerresult_RDCV <- RegularizedSCA::TuckerCoef(my_data_Ttrue, final_RDCV$Tmatrix)
tuckerresult_RDCV$tucker_value # 0.9997192
num_correct(my_data_Ptrue, final_RDCV$Pmatrix[, tuckerresult_RDCV$perm])  # 0.9819149


### 3. BIC and IS
set.seed(1)
result_sim_BICIS <- M2_BIC_IS(my_data, Jk, R, LassoSequence = Lassosequence, GLassoSequence = GLassosequence, NRSTARTS)

Croux_index <- which(result_sim_BICIS$Croux == min(result_sim_BICIS$Croux), arr.ind = T)
Lasso_croux <- Lassosequence[Croux_index[1]] #1.568828
GLasso_croux <- GLassosequence[Croux_index[2]] #1e-07
final_croux <- RegularizedSCA::sparseSCA(my_data, Jk, R, LASSO = Lasso_croux, GROUPLASSO = GLasso_croux, MaxIter = 400, NRSTARTS = 20, method = "component")
final_croux$Pmatrix
tuckerresult_croux <- RegularizedSCA::TuckerCoef(my_data_Ttrue, final_croux$Tmatrix)    
tuckerresult_croux$tucker_value # #0.9872166
num_correct(my_data_Ptrue, final_croux$Pmatrix[, tuckerresult_croux$perm])  # 0.8489362


IS_index <- which(result_sim_BICIS$IS == max(result_sim_BICIS$IS), arr.ind = T)
Lasso_IS <- Lassosequence[IS_index[1]]  # 0.7470612
Glasso_IS <- GLassosequence[IS_index[2]]  #1e-07
final_IS <- RegularizedSCA::sparseSCA(my_data, Jk, R, LASSO = Lasso_IS, GROUPLASSO = Glasso_IS, MaxIter = 400, NRSTARTS = 20, method = "component")
final_IS$Pmatrix
tuckerresult_IS <- RegularizedSCA::TuckerCoef(my_data_Ttrue, final_IS$Tmatrix)    
tuckerresult_IS$tucker_value # 0.9949985
num_correct(my_data_Ptrue, final_IS$Pmatrix[, tuckerresult_IS$perm])  #  0.9159574


### 4. Bolasso
set.seed(1)
result_sim1_Bolasso <- Bolasso_CV(my_data, Jk, R, N_boots, LassoSequence = Lassosequence, GLassoSequence = GLassosequence, NRSTARTS)


tuckerresult_Bolasso <- RegularizedSCA::TuckerCoef(my_data_Ttrue, result_sim1_Bolasso$T_hat)    
tuckerresult_Bolasso$tucker_value  # 0.8920457
num_correct(my_data_Ptrue, result_sim1_Bolasso$P_hat[, tuckerresult_Bolasso$perm])
