# N_dataset <- 20 

set.seed(1)
PropNoise = 0.005
Perc0 = 0.3

I <- 28
J1 <- 144
J2 <- 44
Jk <- c(J1, J2)
R <- 3
N_dataset <- 20  #how many datasets to generate
#NRSTARTS <- 5
#n_rep = 20
#n_seg = 3
#N_boots = 20


################################################################################################################################
# A function for generating data
################################################################################################################################
Data_generation <- function(I, J1, J2, R, PropNoise, Perc0, pre_pro, N_dataset){
  
  R = 3 # in this simulation, we use R = 3
  
  if(missing(I)){
    I = 28
  }
  
  if(missing(J1)){
    J1 = 144
  }
  
  if(missing(J2)){
    J2 = 44
  }
  

  if(missing(pre_pro)){
    pre_pro = "Yes"
  }
  
  Jk <- c(J1, J2)
  sumJk <- sum(J1 + J2)
  
  DATA1 <- matrix(rnorm(I*J1, mean = 0, sd = 1), I, J1)
  DATA2 <- matrix(rnorm(I*J2, mean = 0, sd = 1), I, J2)
  DATA <- cbind(DATA1, DATA2)
  
  svddata <- svd(DATA, R, R)
  Ttrue <- svddata$u
  PTrueC <- as.matrix(svddata$v) %*% diag(svddata$d[1:R])   #note that only the first R eigen values are needed.
  
  PTrueCBlock1 <- PTrueC[1:J1,]
  PTrueCBlock2 <- PTrueC[(J1+1):(J1+J2),]
  
  PTrueCBlock1[, 2] <- 0
  PTrueCBlock2[, 3] <- 0
  
  v <- sample(1:J1, size = round(Perc0*J1), replace = F)
  PTrueCBlock1[, 3][v] <- 0
  v <- sample(1:J2, size = round(Perc0*J2), replace = F)
  PTrueCBlock2[, 2][v] <- 0
  PTrueC_bind <- rbind(PTrueCBlock1, PTrueCBlock2)
  v <- sample(1:sumJk, size = round(Perc0*sumJk), replace = F)
  PTrueC_bind[, 1][v] <- 0
  PTrueCnew <- PTrueC_bind
  
  XTrue <- Ttrue %*% t(PTrueCnew)
  SSXtrue <- sum(XTrue ^ 2)
  
  Noise <- matrix(rnorm(I*(J1+J2), mean = 0, sd = 1), I, J1+J2)
  SSNoise <- sum(Noise ^ 2)
  g <- sqrt(PropNoise*SSXtrue/(SSNoise-PropNoise*SSNoise))
  NoiseNew <- g*Noise
  #SSNoiseNew <- sum(NoiseNew ^ 2)
  Xgenerate <- XTrue + NoiseNew
  #SSXgenerate <- sum(Xgenerate ^ 2)
  #NoiseVSgenerate <- SSNoiseNew/SSXgenerate
  
  if(pre_pro == "Yes"){
    
    ##### Data preprocessing ###################
    #library(RegularizedSCA)
    Data_1 <- RegularizedSCA::pre_process(Xgenerate[, 1:J1], weight = T)
    Data_2 <- RegularizedSCA::pre_process(Xgenerate[, (J1+1):(J1+J2)], weight = T)
    Data_final <- cbind(Data_1, Data_2)
    
  }else{
    
    Data_final <- Xgenerate
    
  }
  
  Data_final <- list(data = Xgenerate, T_mat = Ttrue, P_mat = PTrueCnew)
  return(Data_final)
  
}
###########################################################################################################

k <- 1
while(k <= N_dataset){
  
  my_data_list <- Data_generation(I, J1, J2, R, PropNoise, Perc0, pre_pro = "Yes")
  
  filename <- paste("Data_", k, ".RData", sep = "")
  save(my_data_list, PropNoise, Perc0, file = filename)
  #beta_paramter <- paste("beta_", num_test, ".RData", sep = "")
  #save(beta_pre, beta1, beta2, file = beta_paramter)
  k <- k + 1
  
}
