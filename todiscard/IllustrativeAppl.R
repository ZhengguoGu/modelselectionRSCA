###### Illustrative applications ##################################
library(devtools)
install_github("ZhengguoGu/RegularizedSCA")
library(RegularizedSCA)

###1. family data
# 1) load data

# 2) analyze data
set.seed(1)

data<- cbind(pre_process(family_data[[1]]), pre_process(family_data[[2]]), pre_process(family_data[[3]]))
num_var <- cbind(dim(family_data[[1]])[2], dim(family_data[[2]])[2], dim(family_data[[3]])[2])

R <- 5
Lassosequence <- seq(0.0000001, RegularizedSCA::maxLGlasso(data, num_var, R)$Lasso, length.out = 50)
GLassosequence <- seq(0.0000001, RegularizedSCA::maxLGlasso(data, num_var, R)$Glasso, length.out = 50)

result_sim_BICIS <- M2_BIC_IS(data, num_var, R, LassoSequence = Lassosequence, GLassoSequence = GLassosequence, NRSTARTS=10)

IS_index <- which(result_sim_BICIS$IS == max(result_sim_BICIS$IS), arr.ind = T)
Lasso_IS <- max(Lassosequence[IS_index[1]])
Glasso_IS <- max(GLassosequence[IS_index[2]])
final_IS <- RegularizedSCA::sparseSCA(data, num_var, R, LASSO = Lasso_IS, GROUPLASSO = Glasso_IS, MaxIter = 400, NRSTARTS = 20, method = "component")
familyP <- final_IS$Pmatrix
familyT <-final_IS$Tmatrix

###2. metabolomics data
#1) load data
metab_data <- read.table("C:\\Users\\Zhengguo\\Documents\\modelselectionRSCA\\IllustrativeApplications\\metabolomics.dat", sep = ",")
metab_label <- c(read.csv("C:\\Users\\Zhengguo\\Documents\\modelselectionRSCA\\IllustrativeApplications\\metabolomicsLabel.csv", 
                                 stringsAsFactors = F, header = F))$`V1`
colnames(metab_data) <- metab_label  #note that the data have been pre-processed.

#2) analyze data
set.seed(1)
num_var <- c(144,44)
R <- 5

Lassosequence <- seq(0.0000001, RegularizedSCA::maxLGlasso(metab_data, num_var, R)$Lasso, length.out = 50)
GLassosequence <- seq(0.0000001, RegularizedSCA::maxLGlasso(metab_data, num_var, R)$Glasso, length.out = 50)

metab_result <- M2_BIC_IS(metab_data, num_var, R, LassoSequence = Lassosequence, GLassoSequence = GLassosequence, NRSTARTS=10)
IS_index <- which(metab_result$IS == max(metab_result$IS), arr.ind = T)
Lasso_IS <- max(Lassosequence[IS_index[1]])
Glasso_IS <- max(GLassosequence[IS_index[2]])
final_IS <- RegularizedSCA::sparseSCA(metab_data, num_var, R, LASSO = Lasso_IS, GROUPLASSO = Glasso_IS, MaxIter = 400, NRSTARTS = 20, method = "component")
metaP <- final_IS$Pmatrix
metaT <-final_IS$Tmatrix


###3. hering data
library(RegularizedSCA)
set.seed(1)
ChemPhy <- pre_process(Herring$Herring_ChemPhy)
Sensory <- pre_process(Herring$Herring_Sensory)
herring_data <- cbind(ChemPhy, Sensory)
num_var <- cbind(dim(ChemPhy)[2], dim(Sensory)[2])
R <- 4

Lassosequence <- seq(0.0000001, RegularizedSCA::maxLGlasso(herring_data, num_var, R)$Lasso, length.out = 50)
GLassosequence <- seq(0.0000001, RegularizedSCA::maxLGlasso(herring_data, num_var, R)$Glasso, length.out = 50)

metab_result <- M2_BIC_IS(herring_data, num_var, R, LassoSequence = Lassosequence, GLassoSequence = GLassosequence, NRSTARTS=10)
IS_index <- which(metab_result$IS == max(metab_result$IS), arr.ind = T)
Lasso_IS <- max(Lassosequence[IS_index[1]])
Glasso_IS <- max(GLassosequence[IS_index[2]])
final_IS <- RegularizedSCA::sparseSCA(herring_data, num_var, R, LASSO = Lasso_IS, GROUPLASSO = Glasso_IS, MaxIter = 400, NRSTARTS = 20, method = "component")
herringP <- final_IS$Pmatrix
herringT <-final_IS$Tmatrix

###4. GliomaData


