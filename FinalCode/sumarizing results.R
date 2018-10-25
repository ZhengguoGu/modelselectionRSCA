################################################################################
###################   sumarizing results of the simulation study    ############
################################################################################

################ PART 1: plots for figure 2 and 5 ################################

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


##################### PART 2: plots for figures 3 and 4 ##########################

# (authors' comment: first we have to record the number of non-zero loadings that are correctly identified
# and also the number of zero loadings that are correctly identified.)

numNo0_correct <- function(MATa, MATb){
  num_corr <- sum((MATa != 0) & (MATb != 0))
  return(num_corr)
}
num0_correct <- function(MATa, MATb){
  num_corr <- sum((MATa == 0) & (MATb == 0))
  return(num_corr)
}


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
  tuckerresult <- RegularizedSCA::TuckerCoef(my_data_list$T_mat, ESTIMATED_T[[n_dataset]])
  numNo0_crt_Benchmark[n_dataset]<- numNo0_correct(ESTIMATED_P[[n_dataset]][, tuckerresult$perm], my_data_list$P_mat )
  num0_Benchmark[n_dataset] <- num0_correct(ESTIMATED_P[[n_dataset]][, tuckerresult$perm], my_data_list$P_mat )
  
  #RdCV
  tuckerresult <- RegularizedSCA::TuckerCoef(my_data_list$T_mat, ESTIMATED_TrdCV[[n_dataset]])
  numNo0_crt_RdCV[n_dataset] <- numNo0_correct(ESTIMATED_PrdCV[[n_dataset]][, tuckerresult$perm], my_data_list$P_mat )
  num0_RdCV[n_dataset]<- num0_correct(ESTIMATED_PrdCV[[n_dataset]][, tuckerresult$perm], my_data_list$P_mat )
  
  #BIC
  tuckerresult <- RegularizedSCA::TuckerCoef(my_data_list$T_mat, ESTIMATED_Tbic[[n_dataset]])
  numNo0_crt_BIC[n_dataset] <- numNo0_correct(ESTIMATED_Pbic[[n_dataset]][, tuckerresult$perm], my_data_list$P_mat )
  num0_BIC[n_dataset] <- num0_correct(ESTIMATED_Pbic[[n_dataset]][, tuckerresult$perm], my_data_list$P_mat )
  
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

Ratio_numNo0_crt_Benchmark <- numNo0_crt_Benchmark / numNo0_true 
Ratio_num0_Benchmark <- num0_Benchmark / num0_true
Ratio_numNo0_crt_RdCV <- numNo0_crt_RdCV / numNo0_true
Ratio_num0_RdCV <- num0_RdCV / num0_true
Ratio_numNo0_crt_BIC <- numNo0_crt_BIC / numNo0_true
Ratio_num0_BIC <- num0_BIC / num0_true
Ratio_numNo0_crt_IS <- numNo0_crt_IS / numNo0_true
Ratio_num0_IS <- num0_IS / num0_true
Ratio_numNo0_crt_Bolasso <- numNo0_crt_Bolasso / numNo0_true
Ratio_num0_Bolasso <- num0_Bolasso / num0_true
Ratio_numNo0_crt_Stab <- numNo0_crt_Stab / numNo0_true
Ratio_num0_Stab <- num0_Stab / num0_true


library(ggplot2)
loadings_result <- cbind(Ratio_numNo0_crt_Benchmark,
                         Ratio_numNo0_crt_RdCV,
                         Ratio_numNo0_crt_BIC,
                         Ratio_numNo0_crt_IS,
                         Ratio_numNo0_crt_Bolasso, 
                         Ratio_numNo0_crt_Stab)
colnames(loadings_result) <- c("CV", "RdCV", "BIC", "IS", "BoLasso", "Stab. selection")

boxplot(loadings_result, ylab = "Proportion of non-zero loadings correctly selected", main = "30% noise and 50% zeros", ylim = c(0,1), 
        cex.axis = 1.3, cex.lab = 1.25, cex.main = 1.5)

loadings_result <- cbind(Ratio_num0_Benchmark,
                         Ratio_num0_RdCV,
                         Ratio_num0_BIC,
                         Ratio_num0_IS,
                         Ratio_num0_Bolasso, 
                         Ratio_num0_Stab)
colnames(loadings_result) <- c("CV", "RdCV", "BIC", "IS", "BoLasso", "Stab. selection")

boxplot(loadings_result, ylab = "Proportion of zero loadings correctly identified", main = "30% noise and 50% zeros", ylim = c(0,1), 
        cex.axis = 1.3, cex.lab = 1.25, cex.main = 1.5)


