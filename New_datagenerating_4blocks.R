set.seed(1)

Perc0 = .5   #.3 or .5
PropNoise = 0.005 # .005, .3

I <- 20     
J1 <- 120    
J2 <- 30    
J3 <- 40   
J4 <- 10
Jk <- c(J1, J2, J3, J4)
R <- 3
N_dataset <- 20  #how many datasets to generate



################################################################################################################################
# A function for generating data
################################################################################################################################
Data_generation <- function(I, J1, J2, J3, J4, R, PropNoise, Perc0, N_dataset){
  
  Jk <- c(J1, J2, J3, J4)
  sumJk <- sum(J1 + J2 + J3 + J4)
  
  DATA1 <- matrix(rnorm(I*J1, mean = 0, sd = 1), I, J1)
  DATA2 <- matrix(rnorm(I*J2, mean = 0, sd = 1), I, J2)
  DATA3 <- matrix(rnorm(I*J3, mean = 0, sd = 1), I, J3)
  DATA4 <- matrix(rnorm(I*J4, mean = 0, sd = 1), I, J4)
  DATA <- cbind(DATA1, DATA2, DATA3, DATA4)
  
  svddata <- svd(DATA, R, R)
  Ttrue <- svddata$u
  PTrueC <- as.matrix(svddata$v) %*% diag(svddata$d[1:R])   #note that only the first R eigen values are needed.
  
  PTrueCBlock1 <- PTrueC[1:J1,]
  PTrueCBlock2 <- PTrueC[(J1+1):(J1+J2),]
  PTrueCBlock3 <- PTrueC[(J1+J2+1):(J1+J2+J3),]
  PTrueCBlock4 <- PTrueC[(J1+J2+J3+1):(J1+J2+J3+J4),]
  
  PTrueCBlock1[, 3] <- 0
  PTrueCBlock3[, 2] <- 0
  PTrueCBlock4[, c(2,3)] <- 0
  
  v <- sample(1:J1, size = round(Perc0*J1), replace = F)
  PTrueCBlock1[, 1][v] <- 0
  v <- sample(1:J1, size = round(Perc0*J1), replace = F)
  PTrueCBlock1[, 2][v] <- 0
  
  v <- sample(1:J2, size = round(Perc0*J2), replace = F)
  PTrueCBlock2[, 1][v] <- 0
  v <- sample(1:J2, size = round(Perc0*J2), replace = F)
  PTrueCBlock2[, 2][v] <- 0
  v <- sample(1:J2, size = round(Perc0*J2), replace = F)
  PTrueCBlock2[, 3][v] <- 0
  
  v <- sample(1:J3, size = round(Perc0*J3), replace = F)
  PTrueCBlock3[, 1][v] <- 0
  v <- sample(1:J3, size = round(Perc0*J3), replace = F)
  PTrueCBlock3[, 3][v] <- 0
  
  v <- sample(1:J4, size = round(Perc0*J4), replace = F)
  PTrueCBlock4[, 1][v] <- 0
  
  PTrueC_bind <- rbind(PTrueCBlock1, PTrueCBlock2, PTrueCBlock3, PTrueCBlock4)
  PTrueCnew <- PTrueC_bind
  
  XTrue <- Ttrue %*% t(PTrueCnew)
  SSXtrue <- sum(XTrue ^ 2)
  
  Noise <- matrix(rnorm(I*(J1+J2+J3+J4), mean = 0, sd = 1), I, J1+J2+J3+J4)
  SSNoise <- sum(Noise ^ 2)
  g <- sqrt(PropNoise*SSXtrue/(SSNoise-PropNoise*SSNoise))
  NoiseNew <- g*Noise
  #SSNoiseNew <- sum(NoiseNew ^ 2)
  Xgenerate <- XTrue + NoiseNew
  #SSXgenerate <- sum(Xgenerate ^ 2)
  #NoiseVSgenerate <- SSNoiseNew/SSXgenerate
  
  
  
  Data_final <- list(data = Xgenerate, T_mat = Ttrue, P_mat = PTrueCnew)
  return(Data_final)
  
}
###########################################################################################################

k <- 1
while(k <= N_dataset){
  
  my_data_list <- Data_generation(I, J1, J2, J3, J4, R, PropNoise, Perc0)
  
  filename <- paste("Data_", k, ".RData", sep = "")
  save(my_data_list, PropNoise, Perc0, file = filename)
  #beta_paramter <- paste("beta_", num_test, ".RData", sep = "")
  #save(beta_pre, beta1, beta2, file = beta_paramter)
  k <- k + 1
  
}
