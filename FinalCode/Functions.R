
##########################################################
# martrix_var() calculates the variance in a matrix
##########################################################
matrix_var <- function(X){
  
  X <- data.matrix(X)
  Mvar <- sum(X^2)
  
  return(Mvar)
}


##############################################################################################
# StrucSCA_withIndex() estimates T and P, given the pre-defined structure of P 
##############################################################################################
StrucSCA_withIndex <- function (DATA, Jk, R, P_indexset, MaxIter) {
  DATA <- data.matrix(DATA)
  
  I_Data <- dim(DATA)[1]
  sumJk <- dim(DATA)[2]
  eps <- 10^(-12)
  if (missing(MaxIter)) {
    MaxIter <- 400
  }
  P <- matrix(stats::rnorm(sumJk * R), nrow = sumJk, ncol = R)
  P[P_indexset == 0] <- 0
  Pt <- t(P)
  
  sqP <- P^2
  residual <- sum(DATA^2)
  Lossc <- residual 
  
  conv <- 0
  iter <- 1
  Lossvec <- array()
  while (conv == 0) {
    
    SVD_DATA <- svd(DATA, R, R)
    Tmat <- SVD_DATA$u
    
    residual <- sum((DATA - Tmat %*% Pt)^2)
    Lossu <- residual 
    
    P <- t(DATA) %*% Tmat
    P[P_indexset == 0] <- 0
    Pt <- t(P)
    
    sqP <- P^2
    residual <- sum((DATA - Tmat %*% Pt)^2)
    Lossu2 <- residual
    if (abs(Lossc - Lossu) < 10^(-9)) {
      Loss <- Lossu
      residual <- residual
      P[abs(P) <= 2 * eps] <- 0
      conv <- 1
    }
    else if (iter > MaxIter) {
      Loss <- Lossu
      residual <- residual
      P[abs(P) <= 2 * eps] <- 0
      conv <- 1
    }
    Lossvec[iter] <- Lossu
    iter <- iter + 1
    Lossc <- Lossu2
  }
  return_varselect <- list()
  return_varselect$Pmatrix <- P
  return_varselect$Tmatrix <- Tmat
  return_varselect$Loss <- Loss
  return_varselect$Lossvec <- Lossvec
  return(return_varselect)
}