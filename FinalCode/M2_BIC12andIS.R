####### parameters ###############
# DATA: concatenated data matrix
# Jk: number of variables (a vector, the elements of which correspond to the datablocks)
# R: number of components
# LassoSequence: A sequence of Lasso tuning parameters
# GLassoSequence: A sequence of Group Lasso tuning parameters
# MaxIter: maximum number of iterations
# NRSTARTS: number of multi-starts

load() functions
DATA <- data.matrix(DATA)
n_sub <- dim(DATA)[1]
Jk =
R = 
LassoSequence = 
GLassoSequence = 

MaxIter = 400
NRSTARTS = 2



VarSelect0 <- RegularizedSCA::sparseSCA(DATA, Jk, R, LASSO = 0, GROUPLASSO = 0, MaxIter, NRSTARTS, method = "component")
P_hat0 <- VarSelect0$Pmatrix
T_hat0 <- VarSelect0$Tmatrix

V_0 <- matrix_var(DATA - T_hat0%*%t(P_hat0))  # this is for BIC_Croux and BIC_GUO
error_var <- V_0 / n  #this is for BIC_Guo

V_oo <- matrix_var(DATA)  # this is for Index of sparseness (IS)
V_s <- matrix_var(T_hat0%*%t(P_hat0))  # this is for IS

BIC_Croux <- matrix(NA, length(LassoSequence), length(GLassoSequence))
BIC_Guo <- matrix(NA, length(LassoSequence), length(GLassoSequence))
IS <- matrix(NA, length(LassoSequence), length(GLassoSequence))
for(i in 1:length(LassoSequence)){
  for(j in 1:length(GLassoSequence))
    
    VarSelect <- RegularizedSCA::sparseSCA(DATA, Jk, R, LASSO = LassoSequence[i], GROUPLASSO = GLassoSequence[j], MaxIter, NRSTARTS, method = "component")
    P_hat <- VarSelect$Pmatrix
    T_hat <- VarSelect$Tmatrix
    V_tilde <- matrix_var(DATA - T_hat %*% t(P_hat))
    
    BIC_Croux[i, j] <- V_tilde / V_0 + sum(P_hat != 0) * log(n_sub) / n_sub
    BIC_GUO[i, j] <- V_tilde / error_var + sum(P_hat != 0) * log(n_sub)
     
    V_a <- matrix_var(T_hat %*% t(P_hat))  # this is for IS
    IS[i, j] <- V_a * V_s / V_oo^2 * sum(P_hat == 0) /(sum(Jk) * R)
  
}





