####### parameters ########################################################################
# DATA: concatenated data matrix
# Jk: number of variables (a vector, the elements of which correspond to the datablocks)
# R: number of components
# LassoSequence: A sequence of Lasso tuning parameters
# GLassoSequence: A sequence of Group Lasso tuning parameters
###########################################################################################

load("Functions.R")  #load functions
load("data.RData")

M2_BIC_IS <- function(DATA, Jk, R, LassoSequence, GLassoSequence, NRSTARTS){
  
  DATA <- data.matrix(Data_final)
  n_sub <- dim(DATA)[1]
  
  if(missing(LassoSequence)){
    LassoSequence = seq(0.001, RegularizedSCA::maxLGlasso(DATA, Jk, R)$Lasso, length.out = 20)
  }
  
  if(missing(GLassoSequence)){
    GLassoSequence = seq(0.001, RegularizedSCA::maxLGlasso(DATA, Jk, R)$Glasso, length.out = 20)
  }
  
  if(missing(NRSTARTS)){
    NRSTARTS = 5
  }
  
  VarSelect0 <- RegularizedSCA::sparseSCA(DATA, Jk, R, LASSO = 0, GROUPLASSO = 0, MaxIter = 400, NRSTARTS = NRSTARTS, method = "component")
  P_hat0 <- VarSelect0$Pmatrix
  T_hat0 <- VarSelect0$Tmatrix
  
  V_0 <- sum(DATA - T_hat0%*%t(P_hat0))^2  # this is for BIC_Croux and BIC_GUO
  error_var <- V_0 / n_sub  #this is for BIC_Guo
  
  V_oo <- sum(DATA)^2  # this is for Index of sparseness (IS)
  V_s <- sum(T_hat0%*%t(P_hat0))^2  # this is for IS
  
  BIC_Croux <- matrix(NA, length(LassoSequence), length(GLassoSequence))
  BIC_Guo <- matrix(NA, length(LassoSequence), length(GLassoSequence))
  IS <- matrix(NA, length(LassoSequence), length(GLassoSequence))
  for(i in 1:length(LassoSequence)){
    for(j in 1:length(GLassoSequence)){
      
      VarSelect <- RegularizedSCA::sparseSCA(DATA, Jk, R, LASSO = LassoSequence[i], GROUPLASSO = GLassoSequence[j], MaxIter = 400, NRSTARTS = 5, method = "component")
      P_hat <- VarSelect$Pmatrix
      T_hat <- VarSelect$Tmatrix
      V_tilde <- sum(DATA - T_hat %*% t(P_hat))^2
      
      BIC_Croux[i, j] <- V_tilde / V_0 + sum(P_hat != 0) * log(n_sub) / n_sub  # this is the BIC index prop0sed by Croux et al.
      BIC_Guo[i, j] <- V_tilde / error_var + sum(P_hat != 0) * log(n_sub)      # this is the BIC index proposed by Guo et al.
      
      V_a <- sum(T_hat %*% t(P_hat))^2  # this is for IS
      IS[i, j] <- V_a * V_s / V_oo^2 * sum(P_hat == 0) /(sum(Jk) * R)          # this is index of sparseness
    }
  }
  
  BIC_list <- list(Croux = BIC_Croux, GUO = BIC_Guo, IS = IS)
  return(BIC_list)
}



