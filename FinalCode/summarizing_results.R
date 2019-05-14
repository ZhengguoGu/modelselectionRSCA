################################################################################
###################   sumarizing results of the simulation study    ############
###################               (for revision)                    ############
################################################################################

library(ggplot2)
library(reshape2)
library(gridExtra)

################ PART 1a: Boxplots (2 data blocks), variable correctly selected and zeros correctly identified ################################

# I. summurizing data; 

# Sim_1 0.5% noise and 30% zero ######
#(note: load data by hand)
tucker_result_Sim1 <- cbind(RESULT_BenchmarCV[,1],
                       RESULT_rdCV[, 1],
                       RESULT_BIC[,1],
                       RESULT_IS[,1],
                       RESULT_BoLasso[,1], 
                       RESULT_StabS[,1], 
                       "0.5% noise, 30% zeros")
colnames(tucker_result_Sim1) <- c("CV", "RdCV", "BIC", "IS", "BoLasso", "Stab. selection", "condition")

PL_Sim1 <- cbind(RESULT_BenchmarCV[,2],
            RESULT_rdCV[, 2],
            RESULT_BIC[,2],
            RESULT_IS[,2],
            RESULT_BoLasso[,2],
            RESULT_StabS[,2], 
            "0.5% noise, 30% zeros")
colnames(PL_Sim1) <- c("CV", "RdCV", "BIC", "IS", "BoLasso", "Stab. selection", "condition")
save(file = "sim1.RData", tucker_result_Sim1, PL_Sim1)
#######

# Sim_2 0.5% noise and 50% zero #####################
#(note: load data by hand)
tucker_result_Sim2 <- cbind(RESULT_BenchmarCV[,1],
                            RESULT_rdCV[, 1],
                            RESULT_BIC[,1],
                            RESULT_IS[,1],
                            RESULT_BoLasso[,1], 
                            RESULT_StabS[,1], 
                            "0.5% noise, 50% zeros")
colnames(tucker_result_Sim2) <- c("CV", "RdCV", "BIC", "IS", "BoLasso", "Stab. selection", "condition")

PL_Sim2 <- cbind(RESULT_BenchmarCV[,2],
                 RESULT_rdCV[, 2],
                 RESULT_BIC[,2],
                 RESULT_IS[,2],
                 RESULT_BoLasso[,2],
                 RESULT_StabS[,2], 
                 "0.5% noise, 50% zeros")
colnames(PL_Sim2) <- c("CV", "RdCV", "BIC", "IS", "BoLasso", "Stab. selection", "condition")
save(file = "sim2.RData", tucker_result_Sim2, PL_Sim2)
###########

# Sim_3 30% noise and 30% zero #####################
#(note: load data by hand)
tucker_result_Sim3 <- cbind(RESULT_BenchmarCV[,1],
                            RESULT_rdCV[, 1],
                            RESULT_BIC[,1],
                            RESULT_IS[,1],
                            RESULT_BoLasso[,1], 
                            RESULT_StabS[,1], 
                            "30% noise, 30% zeros")
colnames(tucker_result_Sim3) <- c("CV", "RdCV", "BIC", "IS", "BoLasso", "Stab. selection", "condition")

PL_Sim3 <- cbind(RESULT_BenchmarCV[,2],
                 RESULT_rdCV[, 2],
                 RESULT_BIC[,2],
                 RESULT_IS[,2],
                 RESULT_BoLasso[,2],
                 RESULT_StabS[,2], 
                 "30% noise, 30% zeros")
colnames(PL_Sim3) <- c("CV", "RdCV", "BIC", "IS", "BoLasso", "Stab. selection", "condition")
save(file = "sim3.RData", tucker_result_Sim3, PL_Sim3)
###########

# Sim_4 30% noise and 50% zero #####################
#(note: load data by hand)
tucker_result_Sim4 <- cbind(RESULT_BenchmarCV[,1],
                            RESULT_rdCV[, 1],
                            RESULT_BIC[,1],
                            RESULT_IS[,1],
                            RESULT_BoLasso[,1], 
                            RESULT_StabS[,1], 
                            "30% noise, 50% zeros")
colnames(tucker_result_Sim4) <- c("CV", "RdCV", "BIC", "IS", "BoLasso", "Stab. selection", "condition")

PL_Sim4 <- cbind(RESULT_BenchmarCV[,2],
                 RESULT_rdCV[, 2],
                 RESULT_BIC[,2],
                 RESULT_IS[,2],
                 RESULT_BoLasso[,2],
                 RESULT_StabS[,2], 
                 "30% noise, 50% zeros")
colnames(PL_Sim4) <- c("CV", "RdCV", "BIC", "IS", "BoLasso", "Stab. selection", "condition")
save(file = "sim4.RData", tucker_result_Sim4, PL_Sim4)
###########

load("sim1.RData")
load("sim2.RData")
load("sim3.RData")
load("sim4.RData")

PL <- rbind(PL_Sim1, PL_Sim2, PL_Sim3, PL_Sim4)
PL_final<- data.frame(apply(PL[, 1:6], 2, as.numeric))
PL_final$condition <- PL[, 7]
colnames(PL_final)[c(5, 6)] <- c("BL", "SS")
dat_temp <- melt(PL_final,id.vars="condition", measure.vars=c("CV", "RdCV", "BIC", "IS", "BL", "SS"))


p <- ggplot(dat_temp, aes(x = variable, y = value)) +
  geom_boxplot()+
  scale_y_continuous(name = "Proportion of loadings correctedly selected", limits = c(0, 1)) +
  scale_x_discrete(name = "Variable selection methods") +
  ggtitle("I=80, J1=40, J2=10") +    #do not forget to manually change this.
  theme_bw() +
  theme(plot.title = element_text(size = 14, family = "Tahoma", face = "bold"),
        text = element_text(size = 14, family = "Tahoma"),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 12))+
  facet_grid(. ~ condition)
p

Tucker_results <- rbind(tucker_result_Sim1, tucker_result_Sim2, tucker_result_Sim3, tucker_result_Sim4)
Tucker_final<- data.frame(apply(Tucker_results[, 1:6], 2, as.numeric))
Tucker_final$condition <- Tucker_results[, 7]
colnames(Tucker_final)[c(5, 6)] <- c("BL", "SS")
dat_temp <- melt(Tucker_final,id.vars="condition", measure.vars=c("CV", "RdCV", "BIC", "IS", "BL", "SS"))

p_tucker <- ggplot(dat_temp, aes(x = variable, y = value)) +
  geom_boxplot()+
  scale_y_continuous(name = "Tucker congruence", limits = c(0, 1)) +
  scale_x_discrete(name = "Variable selection methods") +
  ggtitle("I=80, J1=40, J2=10") +    #do not forget to manually change this.
  theme_bw() +
  theme(plot.title = element_text(size = 14, family = "Tahoma", face = "bold"),
        text = element_text(size = 14, family = "Tahoma"),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 12))+
  facet_grid(. ~ condition)
p_tucker
################################################################################################

################ PART 1b: boxplot (2 data blocks), seperately for variable correctly selected and for zeros corrected identified ####################
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

ratio_nonzero_zero <- function(file_names){

  #file_names: it starts with "I_20_J1_120_J2_30", or something like this. 
  fnames <- paste(file_names, "_benchmark_CV.RData", sep = "")
  load(fnames)
  fnames <- paste(file_names, "_BIC_IS.RData", sep = "")
  load(fnames)
  fnames <- paste(file_names, "_BOLASSO.RData", sep = "")
  load(fnames)
  fnames <- paste(file_names, "_RepeatedDCV.RData", sep = "")
  load(fnames)
  fnames <- paste(file_names, "_Stability.RData", sep = "")
  load(fnames)
  
  ###
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
  
  result <- cbind(Ratio_numNo0_crt_Benchmark, Ratio_numNo0_crt_RdCV, Ratio_numNo0_crt_BIC, Ratio_numNo0_crt_IS, Ratio_numNo0_crt_Bolasso, Ratio_numNo0_crt_Stab,
                  Ratio_num0_Benchmark, Ratio_num0_RdCV, Ratio_num0_BIC, Ratio_num0_IS, Ratio_num0_Bolasso, Ratio_num0_Stab)
  return(result)
}

### folder 2block_I20_J120_30
file_heading <- "I_20_J1_120_J2_30"
# note, change directory to the correct subfolder: for example, for "Sim_1 0.5% noise and 30% zero" (see below), change directory to the 0_005noise_0_3zeros folder
# Sim_1 0.5% noise and 30% zero ######
result_sim1 <- ratio_nonzero_zero(file_heading)
save(result_sim1, file = "seperate_sim1.RData")
# Sim_2 0.5% noise and 50% zero ######
result_sim2 <- ratio_nonzero_zero(file_heading)
save(result_sim2, file = "seperate_sim2.RData")
# Sim_3 30% noise and 30% zero #######
result_sim3 <- ratio_nonzero_zero(file_heading)
save(result_sim3, file = "seperate_sim3.RData")
# Sim_4 30% noise and 50% zero #######
result_sim4 <- ratio_nonzero_zero(file_heading)
save(result_sim4, file = "seperate_sim4.RData")

### folder 2block_I20_J40_10
file_heading <- "I_20_J1_40_J2_10"  
# Sim_1 0.5% noise and 30% zero ######
result_sim1 <- ratio_nonzero_zero(file_heading)
save(result_sim1, file = "seperate_sim1.RData")
# Sim_2 0.5% noise and 50% zero ######
result_sim2 <- ratio_nonzero_zero(file_heading)
save(result_sim2, file = "seperate_sim2.RData")
# Sim_3 30% noise and 30% zero #######
result_sim3 <- ratio_nonzero_zero(file_heading)
save(result_sim3, file = "seperate_sim3.RData")
# Sim_4 30% noise and 50% zero #######
result_sim4 <- ratio_nonzero_zero(file_heading)
save(result_sim4, file = "seperate_sim4.RData")

### folder 2block_I80_J40_10
file_heading <- "I_80_J1_40_J2_10"
# Sim_1 0.5% noise and 30% zero ######
result_sim1 <- ratio_nonzero_zero(file_heading)
save(result_sim1, file = "seperate_sim1.RData")
# Sim_2 0.5% noise and 50% zero ######
result_sim2 <- ratio_nonzero_zero(file_heading)
save(result_sim2, file = "seperate_sim2.RData")
# Sim_3 30% noise and 30% zero #######
result_sim3 <- ratio_nonzero_zero(file_heading)
save(result_sim3, file = "seperate_sim3.RData")
# Sim_4 30% noise and 50% zero #######
result_sim4 <- ratio_nonzero_zero(file_heading)
save(result_sim4, file = "seperate_sim4.RData")
