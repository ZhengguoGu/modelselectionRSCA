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
V_0 <- matrix_var(DATA - T_hat0%*%t(P_hat0))

BIC <- matrix(NA, length(LassoSequence), length(GLassoSequence))
for(i in 1:length(LassoSequence)){
  for(j in 1:length(GLassoSequence))
    
    VarSelect <- RegularizedSCA::sparseSCA(DATA, Jk, R, LASSO = LassoSequence[i], GROUPLASSO = GLassoSequence[j], MaxIter, NRSTARTS, method = "component")
    P_hat <- VarSelect$Pmatrix
    T_hat <- VarSelect$Tmatrix
    V_tilde <- matrix_var(DATA - T_hat %*% t(P_hat))
    
    BIC[i, j] <- V_tilde/V_0 + sum(P_hat != 0) * log(n_sub) / n_sub
}





