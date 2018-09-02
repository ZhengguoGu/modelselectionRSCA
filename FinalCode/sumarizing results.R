###############################################################
###################   sumarizing results           ############
###############################################################

################ PART 1: plots ################################

library(ggplot2)


# Sim_1 0.5% noise and 30% zero
  #(note: load data by hand)

tucker_result <- cbind(RESULT_BenchmarCV[,1],
                       RESULT_rdCV[, 1],
                       RESULT_BIC[,1],
                       RESULT_IS[,1],
                       RESULT_BoLasso[,1], 
                       RESULT_StabS[,1])
colnames(tucker_result) <- c("CV", "RdCV", "BIC", "IS", "BoLasso", "Stab. selection")

boxplot(tucker_result, ylab = "Tucker congruence", main = "0.5% noise and 30% zeros", ylim = c(0,1), 
        cex.axis = 1.3, cex.lab = 1.25, cex.main = 1.5)

PL <- cbind(RESULT_BenchmarCV[,2],
            RESULT_rdCV[, 2],
            RESULT_BIC[,2],
            RESULT_IS[,2],
            RESULT_BoLasso[,2],
            RESULT_StabS[,2])
colnames(PL) <- c("CV", "RdCV", "BIC", "IS", "BoLasso", "Stab. selection")
boxplot(PL, ylab = "Proportion of loadings correctly selected", main = "0.5% noise and 30% zeros", ylim = c(0,1), 
        cex.axis = 1.3, cex.lab = 1.25, cex.main = 1.5)


# Sim_2 0.5% noise and 50% zero
#(note: load data by hand)

tucker_result <- cbind(RESULT_BenchmarCV[,1],
                       RESULT_rdCV[, 1],
                       RESULT_BIC[,1],
                       RESULT_IS[,1],
                       RESULT_BoLasso[,1],
                       RESULT_StabS[,1])
colnames(tucker_result) <- c("CV", "RdCV", "BIC", "IS", "BoLasso", "Stab. selection")

boxplot(tucker_result, ylab = "Tucker congruence", main = "0.5% noise and 50% zeros", ylim = c(0,1), 
        cex.axis = 1.3, cex.lab = 1.25, cex.main = 1.5)

PL <- cbind(RESULT_BenchmarCV[,2],
            RESULT_rdCV[, 2],
            RESULT_BIC[,2],
            RESULT_IS[,2],
            RESULT_BoLasso[,2],
            RESULT_StabS[,2])
colnames(PL) <- c("CV", "RdCV", "BIC", "IS", "BoLasso", "Stab. selection")
boxplot(PL, ylab = "Proportion of loadings correctly selected", main = "0.5% noise and 50% zeros", ylim = c(0,1), 
        cex.axis = 1.3, cex.lab = 1.25, cex.main = 1.5)

# Sim_3 5% noise and 30% zero
#(note: load data by hand)

tucker_result <- cbind(RESULT_BenchmarCV[,1],
                       RESULT_rdCV[, 1],
                       RESULT_BIC[,1],
                       RESULT_IS[,1],
                       RESULT_BoLasso[,1],
                       RESULT_StabS[,1])
colnames(tucker_result) <- c("CV", "RdCV", "BIC", "IS", "BoLasso", "Stab. selection")

boxplot(tucker_result, ylab = "Tucker congruence", main = "5% noise and 30% zeros", ylim = c(0,1), 
        cex.axis = 1.3, cex.lab = 1.25, cex.main = 1.5)

PL <- cbind(RESULT_BenchmarCV[,2],
            RESULT_rdCV[, 2],
            RESULT_BIC[,2],
            RESULT_IS[,2],
            RESULT_BoLasso[,2],
            RESULT_StabS[,2])
colnames(PL) <- c("CV", "RdCV", "BIC", "IS", "BoLasso", "Stab. selection")
boxplot(PL, ylab = "Proportion of loadings correctly selected", main = "5% noise and 30% zeros", ylim = c(0,1),
        cex.axis = 1.3, cex.lab = 1.25, cex.main = 1.5)


# Sim_4 5% noise and 50% zero
#(note: load data by hand)

tucker_result <- cbind(RESULT_BenchmarCV[,1],
                       RESULT_rdCV[, 1],
                       RESULT_BIC[,1],
                       RESULT_IS[,1],
                       RESULT_BoLasso[,1],
                       RESULT_StabS[,1])
colnames(tucker_result) <- c("CV", "RdCV", "BIC", "IS", "BoLasso", "Stab. selection")

boxplot(tucker_result, ylab = "Tucker congruence", main = "5% noise and 50% zeros", ylim = c(0,1), 
        cex.axis = 1.3, cex.lab = 1.25, cex.main = 1.5)

PL <- cbind(RESULT_BenchmarCV[,2],
            RESULT_rdCV[, 2],
            RESULT_BIC[,2],
            RESULT_IS[,2],
            RESULT_BoLasso[,2],
            RESULT_StabS[,2])
colnames(PL) <- c("CV", "RdCV", "BIC", "IS", "BoLasso", "Stab. selection")
boxplot(PL, ylab = "Proportion of loadings correctly selected", main = "5% noise and 50% zeros", ylim = c(0,1), 
        cex.axis = 1.3, cex.lab = 1.25, cex.main = 1.5)



# Sim_5 30% noise and 30% zero
#(note: load data by hand)

tucker_result <- cbind(RESULT_BenchmarCV[,1],
                       RESULT_rdCV[, 1],
                       RESULT_BIC[,1],
                       RESULT_IS[,1],
                       RESULT_BoLasso[,1],
                       RESULT_StabS[,1])
colnames(tucker_result) <- c("CV", "RdCV", "BIC", "IS", "BoLasso", "Stab. selection")

boxplot(tucker_result, ylab = "Tucker congruence", main = "30% noise and 30% zeros", ylim = c(0,1),
        cex.axis = 1.3, cex.lab = 1.25, cex.main = 1.5)

PL <- cbind(RESULT_BenchmarCV[,2],
            RESULT_rdCV[, 2],
            RESULT_BIC[,2],
            RESULT_IS[,2],
            RESULT_BoLasso[,2],
            RESULT_StabS[,2])
colnames(PL) <- c("CV", "RdCV", "BIC", "IS", "BoLasso", "Stab. selection")
boxplot(PL, ylab = "Proportion of loadings correctly selected", main = "30% noise and 30% zeros", ylim = c(0,1), 
        cex.axis = 1.3, cex.lab = 1.25, cex.main = 1.5)


# Sim_6 30% noise and 50% zero
#(note: load data by hand)

tucker_result <- cbind(RESULT_BenchmarCV[,1],
                       RESULT_rdCV[, 1],
                       RESULT_BIC[,1],
                       RESULT_IS[,1],
                       RESULT_BoLasso[,1],
                       RESULT_StabS[,1])
colnames(tucker_result) <- c("CV", "RdCV", "BIC", "IS", "BoLasso", "Stab. selection")

boxplot(tucker_result, ylab = "Tucker congruence", main = "30% noise and 50% zeros", ylim = c(0,1),
        cex.axis = 1.3, cex.lab = 1.25, cex.main = 1.5)

PL <- cbind(RESULT_BenchmarCV[,2],
            RESULT_rdCV[, 2],
            RESULT_BIC[,2],
            RESULT_IS[,2],
            RESULT_BoLasso[,2],
            RESULT_StabS[,2])
colnames(PL) <- c("CV", "RdCV", "BIC", "IS", "BoLasso", "Stab. selection")
boxplot(PL, ylab = "Proportion of loadings correctly selected", main = "30% noise and 50% zeros", ylim = c(0,1), 
        cex.axis = 1.3, cex.lab = 1.25, cex.main = 1.5)


##########################################################################################################################


##################### PART 2: Table: number of loadings correctly identified  ##########################

load("BenchmarkCV.RData")
load("RepeatedDCV.RData")
load("BIC_IS.RData")
load("BOLASSO.RData")
load("Stability.RData")
n_dataset <- 1
N_dataset = 20

while(n_dataset <= N_dataset){
  
  filename <- paste("Data_", n_dataset, ".RData", sep = "")
  load(filename)

  result_sim1_BM <- RegularizedSCA::cv_sparseSCA(POST_data, Jk, R, MaxIter = 400, NRSTARTS, Lassosequence, GLassosequence, nfolds = 7, method = "component") 
  tuckerresult <- RegularizedSCA::TuckerCoef(my_data_list$T_mat, result_sim1_BM$T_hat)
  RESULT_BenchmarCV[n_dataset, 1] <- tuckerresult$tucker_value
  RESULT_BenchmarCV[n_dataset, 2] <- num_correct(my_data_list$P_mat, result_sim1_BM$P_hat[, tuckerresult$perm])  
  
  ESTIMATED_P[[n_dataset]] <- result_sim1_BM$P_hat
  ESTIMATED_T[[n_dataset]] <- result_sim1_BM$T_hat
  n_dataset <- n_dataset + 1
}
