#######################################################
####  Illustrative applications 
#######################################################
library(devtools)
install_github("ZhengguoGu/RegularizedSCA")
library(RegularizedSCA)
library(foreach)
library(snow)
library(doSNOW)
library(doRNG)

########## 1. re-analysis of the parent-child relationship survey data
#1) load data
load("family_data.RData")
set.seed(115)


data<- cbind(pre_process(family_data[[1]]), pre_process(family_data[[2]]), pre_process(family_data[[3]]))
num_var <- cbind(dim(family_data[[1]])[2], dim(family_data[[2]])[2], dim(family_data[[3]])[2])
R <- 5 # known based on previous research

Lassosequence <- seq(0.0000001, maxLGlasso(data, num_var, R)$Lasso, length.out = 50)
GLassosequence <- seq(0.0000001, maxLGlasso(data, num_var, R)$Glasso, length.out = 50)

#2) load function M1_repeatedDoubleCV.R

ptm <- proc.time()
result_fam_RDCV <- M1_repeatedDoubleCV(data,  R, num_var, N_cores = 5, LassoSequence = Lassosequence, GLassoSequence = GLassosequence, n_rep=20, n_seg=3, NRSTARTS = 5)

temp_lasso <- as.data.frame(result_fam_RDCV$Lasso)
temp_lasso$Var1 <- sort(as.numeric(levels(temp_lasso$Var1)))
LASSO <- max(temp_lasso[temp_lasso[,2] == max(temp_lasso[,2]),1])  #the first max ensures that the largest Lasso value is chosen, in case more than one lasso value is recommended by M1_repeatedDoubleCV
temp_glasso <- as.data.frame(result_fam_RDCV$GroupLasso)
temp_glasso$Var1 <- sort(as.numeric(levels(temp_glasso$Var1)))
GLASSO <- max(temp_glasso[temp_glasso[,2] == max(temp_glasso[,2]),1])

final_RDCV <- RegularizedSCA::sparseSCA(data, num_var, R, LASSO = LASSO, GROUPLASSO = GLASSO, MaxIter = 400,
                                        NRSTARTS = 20, method = "component")

savetime_family_RdCV <- proc.time() - ptm
final_RDCV$Pmatrix
save(final_RDCV, savetime_family_RdCV, file="family_RdCV.RData")

#3) load function M2_BIC12andIS.R 
ptm <- proc.time()
result_fam_BICIS <- M2_BIC_IS(data, num_var, R, LassoSequence = Lassosequence, GLassoSequence = GLassosequence, NRSTARTS=5)

IS_index <- which(result_fam_BICIS$IS == max(result_fam_BICIS$IS), arr.ind = T)
Lasso_IS <- max(Lassosequence[IS_index[1]])
Glasso_IS <- max(GLassosequence[IS_index[2]])
final_IS <- RegularizedSCA::sparseSCA(data, num_var, R, LASSO = Lasso_IS, GROUPLASSO = Glasso_IS, MaxIter = 400, NRSTARTS = 20, method = "component")
savetime_family_IS <- proc.time() - ptm
save(final_IS, savetime_family_IS, file="family_IS.RData")
