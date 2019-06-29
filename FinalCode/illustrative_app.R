#######################################################
####  Illustrative applications 
#######################################################
library(devtools)
install_github("ZhengguoGu/RegularizedSCA")  #load the lastest R package from Github
library(RegularizedSCA)


######### 1. analysis of the herring data  #####

#1) load the package and data, and pre-process the data
library(RegularizedSCA)
names(Herring)    #the herring data is included in the package RegularizedSCA

ChemPhy <- pre_process(Herring$Herring_ChemPhy)
Sensory <- pre_process(Herring$Herring_Sensory)
herring_data <- cbind(ChemPhy, Sensory)
num_var <- cbind(dim(ChemPhy)[2], dim(Sensory)[2])

R <- 4 # known based on previous research Gu and Van Deun 2018

Lassosequence <- seq(0.0000001, maxLGlasso(herring_data, num_var, R)$Lasso, length.out = 50)
GLassosequence <- seq(0.0000001, maxLGlasso(herring_data, num_var, R)$Glasso, length.out = 50)

#2) load function M2_BIC12andIS.R 
source("M2_BIC12andIS.R")

set.seed(115)
ptm <- proc.time()
result_her_BICIS <- M2_BIC_IS(herring_data, num_var, R, LassoSequence = Lassosequence, GLassoSequence = GLassosequence, NRSTARTS=5)

IS_index <- which(result_her_BICIS$IS == max(result_her_BICIS$IS), arr.ind = T)
Lasso_IS <- max(Lassosequence[IS_index[1]])
Glasso_IS <- max(GLassosequence[IS_index[2]])
final_her_IS <- RegularizedSCA::sparseSCA(herring_data, num_var, R, LASSO = Lasso_IS, GROUPLASSO = Glasso_IS, MaxIter = 400, NRSTARTS = 20, method = "component")
savetime_her_IS <- proc.time() - ptm
save(final_her_IS, savetime_her_IS, file="her_IS.RData")


#4) generate a table 
load("her_IS.RData")

set.seed(115)

final_her_IS <- undoShrinkage(herring_data, R, 
                              final_her_IS$Pmatrix)  
final_her_IS$Pmatrix
final_her_IS$Tmatrix
write.table(final_her_IS$Pmatrix, "final_herring_P.csv", sep = ",")
write.table(final_her_IS$Tmatrix, "final_herring_T.csv", sep = ",")

######## 2. analysis of metabolomics data  ########################
# 1) load data 
MyData <- read.csv(file="metabolomics.csv", header=TRUE, sep=",")

meta1 <- pre_process(MyData[, 1:144], weight = TRUE)
meta2 <- pre_process(MyData[, 145:188], weight = TRUE)
metabolomics_data <- cbind(meta1, meta2)
num_var <- c(144,4)

R <- 5 # known based on previous research Gu and Van Deun 2016

Lassosequence <- seq(0.0000001, maxLGlasso(metabolomics_data , num_var, R)$Lasso, length.out = 50)
GLassosequence <- seq(0.0000001, maxLGlasso(metabolomics_data , num_var, R)$Glasso, length.out = 50)

#2) load function M2_BIC12andIS.R 
source("M2_BIC12andIS.R")

set.seed(115)
ptm <- proc.time()
result_meta_BICIS <- M2_BIC_IS(metabolomics_data, num_var, R, LassoSequence = Lassosequence, GLassoSequence = GLassosequence, NRSTARTS=5)

IS_index <- which(result_meta_BICIS$IS == max(result_meta_BICIS$IS), arr.ind = T)
Lasso_IS <- max(Lassosequence[IS_index[1]])
Glasso_IS <- max(GLassosequence[IS_index[2]])
final_meta_IS <- RegularizedSCA::sparseSCA(metabolomics_data, num_var, R, LASSO = Lasso_IS, GROUPLASSO = Glasso_IS, MaxIter = 400, NRSTARTS = 20, method = "component")
savetime_meta_IS <- proc.time() - ptm
final_meta_IS$Pmatrix
save(final_meta_IS, savetime_meta_IS, file="meta_IS.RData")

set.seed(115)
final_meta_IS <- undoShrinkage(metabolomics_data, R, 
                               final_meta_IS$Pmatrix) 
final_meta_IS$Pmatrix

#4)  We draw a heatmap for IS
Pmat <- final_meta_IS$Pmatrix
keepname <- rownames(Pmat)
short_name <- array()  #some of the variable names are too long, we can shorten the the names (in case we want to) 
for(i in 1:length(keepname)){
  short_name[i] <- substring(keepname[i], first = 1, last = 27)
}

colnames(Pmat) <- c('C1', 'C2', 'C3', 'C4', 'C5')

library(ggplot2)
names <- short_name
component <- colnames(Pmat)
PmatVec <- c(Pmat)
names <- rep(names, 5)
component <- rep(component, each = 188)

# note that part of the ggplot code below is from https://learnr.wordpress.com/2010/01/26/ggplot2-quick-heatmap-plotting/
# which is a website for drawing heatmap using ggplot2. 
Pmat_dataframe <- data.frame(Loadings = PmatVec, Variables = factor(names, ordered = T, levels = short_name), Components = component)

p <- ggplot(Pmat_dataframe, aes(x = Components, y = Variables) )+
  geom_tile(aes(fill = Loadings), colour = "white") +
  scale_fill_gradient2(low="green", mid = "black", high = "red") 

p + theme_grey(base_size = 15) + labs(x = "", y = "") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0))

########## 3. re-analysis of the parent-child relationship survey data  #########################
#1) load data
load("family_data.RData")

data<- cbind(pre_process(family_data[[1]]), pre_process(family_data[[2]]), pre_process(family_data[[3]]))
num_var <- cbind(dim(family_data[[1]])[2], dim(family_data[[2]])[2], dim(family_data[[3]])[2])
R <- 5 # known based on previous research

Lassosequence <- seq(0.0000001, maxLGlasso(data, num_var, R)$Lasso, length.out = 50)
GLassosequence <- seq(0.0000001, maxLGlasso(data, num_var, R)$Glasso, length.out = 50)


#2) load function M2_BIC12andIS.R 
source("M2_BIC12andIS.R")

set.seed(115)
ptm <- proc.time()
result_fam_BICIS <- M2_BIC_IS(data, num_var, R, LassoSequence = Lassosequence, GLassoSequence = GLassosequence, NRSTARTS=5)

IS_index <- which(result_fam_BICIS$IS == max(result_fam_BICIS$IS), arr.ind = T)
Lasso_IS <- max(Lassosequence[IS_index[1]])
Glasso_IS <- max(GLassosequence[IS_index[2]])
final_IS <- RegularizedSCA::sparseSCA(data, num_var, R, LASSO = Lasso_IS, GROUPLASSO = Glasso_IS, MaxIter = 400, NRSTARTS = 20, method = "component")
savetime_family_IS <- proc.time() - ptm
save(final_IS, savetime_family_IS, file="family_IS.RData")

#4) Undo the shrinkage and generate a table 
# In Table 2, the component loading matrix obtained from Gu and Van Deun 2018, the authors undo the shrinkage, Hence, we undo the shrinkage here. 

load("familytarget.RData")  # this is the T matrix of the Family data from Gu and Van deun 2018.

perm2 <- RegularizedSCA::TuckerCoef(family_target, final_IS$Tmatrix)$perm
final_IS_result <- final_IS$Pmatrix[, perm2]   #final P matrix for IS


set.seed(115)
final_fam_IS <- undoShrinkage(data, R = 5, 
                              final_IS_result)  #position of components changed so as to be compared to the results by RdCV
final_fam_IS$Pmatrix
write.table(final_fam_IS$Pmatrix, "final_fam.csv", sep = ",")

#!note, the final_fam_IS$Pmatrix is different from the the reported matrix in paper in terms of signs. because as SCA is invariant with respect to the sign, we manually changed the signs
# so that it is easier to compared to Figure 2. The interpretation does not change as a result of change of signs. 





