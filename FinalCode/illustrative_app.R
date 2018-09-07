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
load(family_data.RData)
set.seed(115)

#2) load function M1_repeatedDoubleCV.R and pre-process data
data<- cbind(pre_process(family_data[[1]]), pre_process(family_data[[2]]), pre_process(family_data[[3]]))
num_var <- cbind(dim(family_data[[1]])[2], dim(family_data[[2]])[2], dim(family_data[[3]])[2])
R <- 5 # known based on previous research

Lassosequence <- seq(0.0000001, maxLGlasso(data, num_var, R)$Lasso, length.out = 50)
GLassosequence <- seq(0.0000001, maxLGlasso(data, num_var, R)$Glasso, length.out = 50)

ptm <- proc.time()
result_fam_RDCV <- M1_repeatedDoubleCV(data,  R, num_var, N_cores = 10, LassoSequence = Lassosequence, GLassoSequence = GLassosequence, n_rep=20, n_seg=3, NRSTARTS = 5)
proc.time() - ptm
