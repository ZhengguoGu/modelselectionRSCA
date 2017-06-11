######## 1. CDfriedmanV2 ##########################

#Variable selection with Lasso and Group Lasso (component-wise)

CDfriedmanV2 <- function (DATA, Jk, R, LASSO, GROUPLASSO, MaxIter) 
{
  DATA <- data.matrix(DATA)
  I_Data <- dim(DATA)[1]
  sumJk <- dim(DATA)[2]
  eps <- 10^(-12)
  if (missing(MaxIter)) {
    MaxIter <- 400
  }
  P <- matrix(rnorm(sumJk * R), nrow = sumJk, ncol = R)
  Pt <- t(P)
  absP <- abs(P)
  pen_l <- LASSO * sum(absP)
  sqP <- P^2
  L <- 1
  pen_g <- 0
  for (i in 1:length(Jk)) {
    U <- L + Jk[i] - 1
    sqrtsumP <- sqrt(colSums(sqP[L:U, ]))/sqrt(Jk[i])
    pen_g <- pen_g + GROUPLASSO * sum(sqrtsumP) * Jk[i]
    L <- U + 1
  }
  residual <- sum(DATA^2)
  Lossc <- residual + pen_l + pen_g
  conv <- 0
  iter <- 1
  Lossvec <- array()
  while (conv == 0) {
    if (LASSO == 0 & GROUPLASSO == 0) {
      SVD_DATA <- svd(DATA, R, R)
      Tmat <- SVD_DATA$u
    }
    else {
      A <- Pt %*% t(DATA)
      SVD_DATA <- svd(A, R, R)
      Tmat <- SVD_DATA$v %*% t(SVD_DATA$u)
    }
    residual <- sum((DATA - Tmat %*% Pt)^2)
    Lossu <- residual + pen_l + pen_g
    if (LASSO == 0 & GROUPLASSO == 0) {
      P <- t(DATA) %*% Tmat
      Pt <- t(P)
    }
    else {
      L <- 1
      for (i in 1:length(Jk)) {
        U <- L + Jk[i] - 1
        Pt_1 <- Pt[, c(L:U)]
        data <- DATA[, c(L:U)]
        if (sum(abs(Pt_1)) != 0) {
          for (r in 1:R) {
            copy_Pt1 <- Pt_1
            copy_Pt1[r, ] <- 0
            matrix_R <- data - Tmat %*% copy_Pt1
            matrix_R <- t(matrix_R)
            S_2Vec_Lambda <- soft_th(2 * matrix_R %*% 
                                       Tmat[, r], LASSO)
            S_2Vec_Lambda_norm <- sqrt(sum(S_2Vec_Lambda^2))
            s_l2 <- 0.5 - 0.5 * GROUPLASSO * sqrt(Jk[i])/S_2Vec_Lambda_norm
            if (s_l2 <= 0 | S_2Vec_Lambda_norm == 0) {
              Pt[r, c(L:U)] <- 0
            }
            else if (s_l2 > 0) {
              Pt[r, c(L:U)] <- s_l2 * c(S_2Vec_Lambda)
            }
          }
        }
        L <- U + 1
      }
      P <- t(Pt)
    }
    pen_l <- LASSO * sum(abs(P))
    sqP <- P^2
    L <- 1
    pen_g <- 0
    for (i in 1:length(Jk)) {
      U <- L + Jk[i] - 1
      sqrtsumP <- sqrt(colSums(sqP[L:U, ]))/sqrt(Jk[i])
      pen_g <- pen_g + GROUPLASSO * sum(sqrtsumP) * Jk[i]
      L <- U + 1
    }
    residual <- sum((DATA - Tmat %*% Pt)^2)
    Lossu2 <- residual + pen_l + pen_g
    if (abs(Lossc - Lossu) < 10^(-9)) {
      Loss <- Lossu
      residual <- residual
      lassopen <- pen_l
      Glassopen <- pen_g
      P[abs(P) <= 2 * eps] <- 0
      conv <- 1
    }
    else if (iter > MaxIter | LASSO == 0) {
      Loss <- Lossu
      residual <- residual
      lassopen <- pen_l
      Glassopen <- pen_g
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

###############################################

#. multistart procedure.

BIC_SCA <- function (DATA, Jk, R, LASSO, GROUPLASSO, MaxIter, NRSTARTS) 
{
  n_sub <- dim(DATA)[1]
  
  if (missing(NRSTARTS)) {
    NRSTARTS <- 20
  }

  if (R == 1) {
    stop("Parameter R = 1 is not allowed.")
  }
  Pout3d <- list()
  Tout3d <- list()
  LOSS <- array()
  LOSSvec <- list()
  for (n in 1:NRSTARTS) {
          
    VarSelectResult <- CDfriedmanV2(DATA, Jk, R, LASSO, 
                                      GROUPLASSO, MaxIter)
    
    Pout3d[[n]] <- VarSelectResult$Pmatrix
    Tout3d[[n]] <- VarSelectResult$Tmatrix
    LOSS[n] <- VarSelectResult$Loss
    LOSSvec[[n]] <- VarSelectResult$Lossvec
  }
  k <- which(LOSS == min(LOSS))
  if (length(k) > 1) {
    pos <- sample(1:length(k), 1)
    k <- k[pos]
  }
  
  Res_mat <- DATA - Tout3d[[k]] %*% t(Pout3d[[k]])
  variance_sparse <- sum((apply(Res_mat, 2, var)) * (n_sub-1)/n_sub)  #empirical variance = population variance
  
  
  VarSelectResult0 <- CDfriedmanV2(DATA, Jk, R, LASSO=0, 
                                  GROUPLASSO=0, MaxIter)
  Tmat0 <- VarSelectResult0$Tmatrix
  Pmat0 <- VarSelectResult0$Pmatrix
  Res_mat0 <- DATA - Tmat0 %*% t(Pmat0)
  variance_0 <- sum((apply(Res_mat0, 2, var)) * (n_sub-1)/n_sub)  #empirical variance = population variance
  
  bic <- variance_sparse/VarSelectResult0 + sum(Pout3d[[k]]!=0) * ln(n_sub)/n_sub
  
    
  return_varselect <- list()
  return_varselect$Pmatrix <- Pout3d[[k]]
  return_varselect$Tmatrix <- Tout3d[[k]]
  return_varselect$Lossvec <- LOSSvec
  return_varselect$P0 <- Pmat0
  return_varselect$T0 <- Tmat0
  return_varselect$BIC <- bic
  return(return_varselect)
}

