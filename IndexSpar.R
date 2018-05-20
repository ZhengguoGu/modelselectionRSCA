#Index of Sparseness (Gajjar, Kulahci, and Palazoglu, 2017; Tendafilov, 2014)

# Input parameters:
# DATA: Concatenated data matrix
# Jk: A vector containing number of variables in the concatinated data matrix
# R: Number of components (R>=2)
# LassoSequence: A vector of lasso tuning parameter values in accending order
# GlassoSequence: A vector of Group Lasso tuning parameter values in accending order

IndexSpar <- function(DATA, Jk, R, LassoSequence, GlassoSequence){
  
  DATA <- data.matrix(DATA)
  
  V_o <- matrix_var(DATA)
  
  result <- sparseSCA(DATA, Jk, R, Lasso = 0, GROUPLASSO = 0, NRSTARTS=1)
  T_hat <- result$Tmatrix
  P_hat <- result$Pmatrix 
  V_s <- matrix_var(T_hat%*%t(P_hat))
  
  V_so2JkR <- V_s/(V_o^2 * sum(Jk) * R)
    
  IS <- matrix(NA, length(LassoSequence), length(GlassoSequence))
  
  for(j in 1:length(GlassoSequence)){
    for(i in 1:length(LassoSequence)){
      
      result <- sparseSCA(DATA, Jk, R, LassoSequence[i], GlassoSequence[j], NRSTARTS=1)
      T_hat <- result$Tmatrix
      P_hat <- result$Pmatrix
      
      V_a <- matrix_var(T_hat%*%t(P_hat))
      
      N_0loading <- sum(P_hat==0)
      
      IS[i,j] <- V_a * N_0loading/V_so2JkR 
    }
  }
  
  
  return(IS)
  
}