# this is to test repeated double cross validation
library(foreach)
library(psychometric)
library(doSNOW)
library(doRNG)


set.seed(112)

I <- 28
J1 <- 144
J2 <-44
Jk <- c(J1, J2)
sumJk <- sum(J1 + J2)
R <- 5


PropNoise <- 0.05
Perc0 <- 0.3

NRSTARTS <- 20
Ndataset <- 50
MAXITER <- 400

DATA1 <- matrix(rnorm(I*J1, mean = 0, sd = 1), I, J1)
DATA2 <- matrix(rnorm(I*J2, mean = 0, sd = 1), I, J2)
DATA <- cbind(DATA1, DATA2)

svddata <- svd(DATA, R, R)
Ttrue <- svddata$u
PTrueC <- as.matrix(svddata$v) %*% diag(svddata$d[1:R])   #note that only the first R eigen values are needed.

PTrueCBlock1 <- PTrueC[1:J1,]
PTrueCBlock2 <- PTrueC[(J1+1):(J1+J2),]

v1 <- c(1, 2, 3)
PTrueCBlock1[, v1] <- 0
v2 <- c(4, 5)
PTrueCBlock2[, v2] <- 0



PTrueCBlock1_vec <- as.vector(PTrueCBlock1[, v2])
v <- sample(1:(J1*2), size = round(Perc0*(J1*2)), replace=F)
PTrueCBlock1_vec[v] <- 0
PTrueCBlock1[, v2] <- matrix(PTrueCBlock1_vec, nrow = J1, ncol = 2)

PTrueCBlock2_vec <- as.vector(PTrueCBlock2[, v1])
v <- sample(1:(J2*3), size = round(Perc0*(J2*3)), replace=F)
PTrueCBlock2_vec[v] <- 0
PTrueCBlock2[, v1] <- matrix(PTrueCBlock2_vec, nrow = J2, ncol = 3)

PTrueCnew <- rbind(PTrueCBlock1, PTrueCBlock2)

XTrue <- Ttrue %*% t(PTrueCnew)
SSXtrue <- sum(XTrue ^ 2)

Noise <- matrix(rnorm(I*(J1+J2), mean = 0, sd = 1), I, J1+J2)
SSNoise <- sum(Noise ^ 2)
g <- sqrt(PropNoise*SSXtrue/(SSNoise-PropNoise*SSNoise))
NoiseNew <- g*Noise
SSNoiseNew <- sum(NoiseNew ^ 2)
Xgenerate <- XTrue + NoiseNew
SSXgenerate <- sum(Xgenerate ^ 2)
NoiseVSgenerate <- SSNoiseNew/SSXgenerate


result <- rdCV_RSCA(DATA=Xgenerate, Jk=Jk, R=R, 
                    LassoSequence = seq(.01, 2, by = 0.1), 
                    GLassoSequence = seq(.5, 2, by = 0.1), 
                    n_rep = 5, n_seg = 3, 
                    MaxIter = MAXITER, NRSTARTS = NRSTARTS, 
                    nfolds = 5)


#####################################

OptimumLasso <- .11
OptimumGLasso <- .5
method <- "component"
Pout3dopt <- list()
Tout3dopt <- list()
LOSSopt <- array()
for (n in 1:NRSTARTS){
  if(method == "datablock"){
    VarSelectResultOpt <- CDfriedmanV1(DATA=Xgenerate, Jk, R, OptimumLasso, OptimumGLasso, MaxIter = MAXITER)
  }else if (method == "component"){
    VarSelectResultOpt <- CDfriedmanV2(DATA=Xgenerate, Jk, R, OptimumLasso, OptimumGLasso, MaxIter = MAXITER)
  }
  Pout3dopt[[n]] <- VarSelectResultOpt$Pmatrix
  Tout3dopt[[n]] <- VarSelectResultOpt$Tmatrix
  LOSSopt[n] <- VarSelectResultOpt$Loss
}
k <- which(LOSSopt == min(LOSSopt))
if (length(k)>1){
  pos <- sample(1:length(k), 1)
  k <- k[pos]
}
PoutFinal <- Pout3dopt[[k]]
ToutFinal <- Tout3dopt[[k]]
TuckerResults <- TuckerCoef(Ttrue, ToutFinal)
TuckerValues <- TuckerResults$tucker_value
PoutBest <- PoutFinal[, TuckerResults$perm]

indSelectedC <- which(PoutBest != 0)
indDropedC <- which(PoutBest == 0)
Proportion <- (sum(PTrueCnew[indSelectedC] != 0) + sum(PTrueCnew[indDropedC] == 0))/(sumJk*R)
