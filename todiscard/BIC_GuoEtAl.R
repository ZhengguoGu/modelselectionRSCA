# Another BIC: Guo, James, Levina, Michailidis, & Zhu (2010)

# Input parameters:
# DATA: Concatenated data matrix
# T_hat: Estimated T matrix (given Lambda_L, and Lambda_G)
# P_hat: Estimated P matrix (given Lambda_L, and Lambda_G)
# T_hat0: Estimated T matrix (given Lambda_L=0, and Lambda_G=0)
# P_hat0: Estimated P matrix (given Lambda_L=0, and Lambda_G=0)

BIC_GuoEtAl <- function(DATA, T_hat, P_hat, T_hat0, P_hat0){
  
  I_DATA <- dim(DATA)[1]
  error_var <- sum((DATA - T_hat0 %*% t(P_hat0))^2)/I_DATA
  
  df_RSCA <- sum(P_hat ==0)
  
  BIC_result <- sum((DATA - T_hat %*% t(P_hat))^2)/error_var + log(I_DATA) * df_RSCA
  
  return(BIC_result)
}