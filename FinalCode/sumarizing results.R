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
numNo0_correct <- function(MATa, MATb){
  num_corr <- sum((MATa != 0) & (MATb != 0))
  return(num_corr)
}
num0_correct <- function(MATa, MATb){
  num_corr <- sum((MATa == 0) & (MATb == 0))
  return(num_corr)
}

truncate_threshold <- 0.0001  #I noticed that many estimated loadings were very small, but they were not exactly zero. it would be interesting to see
                              # if we force loadings < truncate_threshold to be zeros. 

load("BenchmarkCV.RData")
load("RepeatedDCV.RData")
load("BIC_IS.RData")
load("BOLASSO.RData")
load("Stability.RData")
n_dataset <- 1

numNo0_true <- array()
num0_true <- array()

numNo0_crt_Benchmark <- array()   #number of non-zero loadings 
num0_Benchmark <- array()
numNo0_crt_RdCV <- array() 
num0_RdCV <- array()
numNo0_crt_BIC <- array() 
num0_BIC <- array()
numNo0_crt_IS <- array() 
num0_IS <- array()
numNo0_crt_Bolasso <- array() 
num0_Bolasso <- array()
numNo0_crt_Stab <- array() 
num0_Stab <- array()

 
while(n_dataset <= 20){
  
  filename <- paste("Data_", n_dataset, ".RData", sep = "")
  load(filename)

  numNo0_true[n_dataset] <- sum(my_data_list$P_mat !=0)
  num0_true[n_dataset] <- sum(my_data_list$P_mat ==0)
  
  #BenchmarkCV
  tuckerresult <- RegularizedSCA::TuckerCoef(my_data_list$T_mat, ESTIMATED_T[[n_dataset]]$T_hat)
  numNo0_crt_Benchmark[n_dataset]<- numNo0_correct(ESTIMATED_P[[n_dataset]][[1]][, tuckerresult$perm], my_data_list$P_mat )
  num0_Benchmark[n_dataset] <- num0_correct(ESTIMATED_P[[n_dataset]][[1]][, tuckerresult$perm], my_data_list$P_mat )
  
  #RdCV
  tuckerresult <- RegularizedSCA::TuckerCoef(my_data_list$T_mat, ESTIMATED_TrdCV[[n_dataset]])
  numNo0_crt_RdCV[n_dataset] <- numNo0_correct(ESTIMATED_PrdCV[[n_dataset]][, tuckerresult$perm], my_data_list$P_mat )
  num0_RdCV[n_dataset]<- num0_correct(ESTIMATED_PrdCV[[n_dataset]][, tuckerresult$perm], my_data_list$P_mat )
  
  #BIC
  tuckerresult <- RegularizedSCA::TuckerCoef(my_data_list$T_mat, ESTIMATED_Tbic[[n_dataset]])
  numNo0_crt_BIC[n_dataset] <- numNo0_correct(ESTIMATED_Pbic[[n_dataset]][, tuckerresult$perm], my_data_list$P_mat )
  num0_BICv <- num0_correct(ESTIMATED_Pbic[[n_dataset]][, tuckerresult$perm], my_data_list$P_mat )
  
  #IS
  tuckerresult <- RegularizedSCA::TuckerCoef(my_data_list$T_mat, ESTIMATED_TIS[[n_dataset]])
  numNo0_crt_IS[n_dataset] <- numNo0_correct(ESTIMATED_PIS[[n_dataset]][, tuckerresult$perm], my_data_list$P_mat )
  num0_IS[n_dataset] <- num0_correct(ESTIMATED_PIS[[n_dataset]][, tuckerresult$perm], my_data_list$P_mat )
  
  #Bolasso
  tuckerresult <- RegularizedSCA::TuckerCoef(my_data_list$T_mat, ESTIMATED_Tbolasso[[n_dataset]])
  numNo0_crt_Bolasso[n_dataset] <- numNo0_correct(ESTIMATED_Pbolasso[[n_dataset]][, tuckerresult$perm], my_data_list$P_mat )
  num0_Bolasso[n_dataset] <- num0_correct(ESTIMATED_Pbolasso[[n_dataset]][, tuckerresult$perm], my_data_list$P_mat )
  
  #Stability Selection
  tuckerresult <- RegularizedSCA::TuckerCoef(my_data_list$T_mat, ESTIMATED_TStabS[[n_dataset]])
  numNo0_crt_Stab[n_dataset] <- numNo0_correct(ESTIMATED_PStabS[[n_dataset]][, tuckerresult$perm], my_data_list$P_mat )
  num0_Stab[n_dataset] <- num0_correct(ESTIMATED_PStabS[[n_dataset]][, tuckerresult$perm], my_data_list$P_mat )
  
  n_dataset <- n_dataset + 1
}


