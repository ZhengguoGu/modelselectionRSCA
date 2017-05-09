##################################  Functions for rdCV.R ###############################################

######## 1. CDfriedmanV2 ##########################

#Variable selection with Lasso and Group Lasso (component-wise)
CDfriedmanV2 <- function(DATA, Jk, R, LASSO, GROUPLASSO, MaxIter){
  
  DATA <- data.matrix(DATA)
  I_Data <- dim(DATA)[1]
  sumJk <- dim(DATA)[2]
  eps <- 10^(-12)
  
  if(missing(MaxIter)){
    MaxIter <- 400
  }
  
  #initialize P
  P <- matrix(rnorm(sumJk * R), nrow = sumJk, ncol = R)
  Pt <- t(P)
  
  absP <- abs(P)
  pen_l <- LASSO * sum(absP)
  sqP <- P^2
  L <- 1
  pen_g <- 0
  for (i in 1:length(Jk)){
    
    U <- L + Jk[i] - 1
    sqrtsumP <- sqrt(colSums(sqP[L:U, ]))/sqrt(Jk[i])
    pen_g <- pen_g + GROUPLASSO * sum(sqrtsumP) * Jk[i]
    L <- U + 1
  }
  
  residual <- sum(DATA ^ 2)
  Lossc <- residual + pen_l + pen_g
  
  conv <- 0
  iter <- 1
  Lossvec <- array()
  
  while (conv == 0){
    
    #update Tmat, note that Tmax refers to T matrix
    if (LASSO == 0 & GROUPLASSO == 0){
      SVD_DATA <- svd(DATA, R, R)  #note this is different from the matlab svds function. need to test it!!
      Tmat <- SVD_DATA$u
    }
    else {
      A <- Pt %*% t(DATA)
      SVD_DATA <- svd(A, R, R)
      Tmat <- SVD_DATA$v %*% t(SVD_DATA$u)
    }
    
    residual <- sum((DATA - Tmat %*% Pt)^2)
    Lossu <- residual + pen_l + pen_g
    
    #update P
    if (LASSO == 0 & GROUPLASSO == 0){
      
      P <- t(DATA) %*% Tmat
      Pt <- t(P)
      
    }
    else{
      
      L <- 1
      for (i in 1:length(Jk)){ #iterate over groups
        
        U <- L + Jk[i] - 1
        Pt_1 <- Pt[ ,c(L:U)]
        data <- DATA[ ,c(L:U)]
        
        
        if (sum(abs(Pt_1)) != 0){
          # calculate S(2(t_r \otimes I)' Vec(R), Lambda_L
          
          for(r in 1:R){
            
            copy_Pt1 <- Pt_1
            copy_Pt1[r, ] <- 0
            matrix_R <- data - Tmat %*% copy_Pt1
            matrix_R <- t(matrix_R)
            S_2Vec_Lambda <- soft_th(2 * matrix_R %*% Tmat[,r], LASSO)
            S_2Vec_Lambda_norm <- sqrt(sum(S_2Vec_Lambda^2))
            
            s_l2 <- 0.5 - 0.5 * GROUPLASSO * sqrt(Jk[i])/S_2Vec_Lambda_norm
            
            if(s_l2 <= 0 | S_2Vec_Lambda_norm==0){
              Pt[r ,c(L:U)] <- 0
            } else if(s_l2 > 0){
              Pt[r ,c(L:U)] <- s_l2 * c(S_2Vec_Lambda)
            }
          }
        }
        L <- U + 1
      }
      
      P <- t(Pt)
    }
    
    
    pen_l <- LASSO*sum(abs(P))
    sqP <- P^2
    L <- 1
    pen_g <- 0
    for (i in 1:length(Jk)){
      U <- L + Jk[i] - 1
      sqrtsumP <- sqrt(colSums(sqP[L:U, ]))/sqrt(Jk[i])
      pen_g <- pen_g + GROUPLASSO * sum(sqrtsumP) * Jk[i]
      L <- U + 1
    }
    
    residual <- sum((DATA - Tmat %*% Pt)^2)
    Lossu2 <- residual + pen_l + pen_g
    
    #check convergence
    if (abs(Lossc-Lossu)< 10^(-9)) {
      Loss <- Lossu
      residual <- residual
      lassopen <- pen_l
      Glassopen <- pen_g
      P[abs(P) <= 2 * eps] <- 0
      conv <- 1
    } else if (iter > MaxIter | LASSO == 0){
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


######### 2. cv_sparseSCA ###################################################
#Note: This code is originally from the package RSCA, but has been adjusted for rdCV_RSCA

cv_sparseSCA <- function(DATA, Jk, R, MaxIter, NRSTARTS, LassoSequence, GLassoSequence, nfolds, method){
  
  DATA <- data.matrix(DATA)
  
  PRESS <- matrix(0, length(LassoSequence), length(GLassoSequence))
  se_MSE <- matrix(0, length(LassoSequence), length(GLassoSequence))
  varselected <- matrix(0, length(LassoSequence), length(GLassoSequence))
  percentRemove <- 1/nfolds
  
  ii <- 0
  while(ii != 1){ #this procedure is to make sure that the training sample do not have an entire row/column of NA's
    
    randmat <- matrix(runif(nrow(DATA) * ncol(DATA)), ncol = ncol(DATA))
    jj <- 0
    for (i in 1:nfolds){
      ToRemove <- ((i - 1) * percentRemove < randmat) & (randmat < i * percentRemove) # this idea is from PMA package
      DATArm <- DATA
      DATArm[ToRemove] <- NA
      
      if(ncol(DATArm) %in% rowSums(is.na(DATArm))){
        break
      } else if(nrow(DATArm) %in% colSums(is.na(DATArm))){
        break
      }else{
        jj <- jj+1
      }
    }
    
    if(jj == nfolds){
      ii <- 1
    }
  }
  
  
  for (g in 1: length(GLassoSequence)){
    
    cat(sprintf("\n GLassoSequence: %s", GLassoSequence[g]))
    
    for (l in 1:length(LassoSequence)){
      
      cat(sprintf("\n LassoSequence: %s", LassoSequence[l]))
      
      if(method == "datablock"){
        Forvarselected <- CDfriedmanV1(DATA, Jk, R, LassoSequence[l], GLassoSequence[g], MaxIter)
      }else if (method == "component"){
        Forvarselected <- CDfriedmanV2(DATA, Jk, R, LassoSequence[l], GLassoSequence[g], MaxIter)
      }
      varselected[l,g] <- sum(Forvarselected$Pmatrix != 0)  #how many variables in P have been selected?
      
      error_x <- array()
      for (i in 1:nfolds){
        
        cat(sprintf("\n Inner loop: %s", i))
        
        ToRemove <- ((i - 1) * percentRemove < randmat) & (randmat < i * percentRemove) # this idea is from PMA package
        DATArm <- DATA
        DATArm[ToRemove] <- NA
        
        for(c in 1:ncol(DATA)){
          indexc <- !is.na(DATArm[, c])
          DATArm[, c][!indexc] <- mean(DATArm[, c][indexc]) #missing values are replaced by column means
        }
        
        Pout3d <- list()
        Tout3d <- list()
        LOSS <- array()

        for (n in 1:NRSTARTS){
          if(method == "datablock"){
            VarSelectResult <- CDfriedmanV1(DATArm, Jk, R, LassoSequence[l], GLassoSequence[g], MaxIter)
          }else if (method == "component"){
            VarSelectResult <- CDfriedmanV2(DATArm, Jk, R, LassoSequence[l], GLassoSequence[g], MaxIter)
          }
          Pout3d[[n]] <- VarSelectResult$Pmatrix
          Tout3d[[n]] <- VarSelectResult$Tmatrix
          LOSS[n] <- VarSelectResult$Loss
        }
        k <- which(LOSS == min(LOSS))
        if (length(k)>1){
          pos <- sample(1:length(k), 1)
          k <- k[pos]
        }
        PoutBest <- Pout3d[[k]]
        ToutBest <- Tout3d[[k]]
        
        DATA_hat <- ToutBest%*%t(PoutBest)
        error_x[i] <- sum((DATA[ToRemove] - DATA_hat[ToRemove])^2)
        
      }
      
      
      PRESS[l,g] <- sum(error_x)/nfolds
      se_MSE[l,g] <- sd(error_x)/sqrt(nfolds)
    }
  }
  
  vec_PRESS <- c(PRESS)
  vec_se <- c(se_MSE)
  vec_varsel <- c(varselected)
 

    
  lowestPress <- min(vec_PRESS)
  if(length(which(vec_PRESS == lowestPress))>1){
    #this means there are more than 1 candidate. 
    max_SE <- max(vec_se[which(vec_PRESS == lowestPress)])  #Note, this is our decision rule: we choose the highest se, so that P may be sparser. 
  } else{
    max_SE <- vec_se[which(vec_PRESS == lowestPress)]
  } 
  
  lowestplus1SE <- lowestPress + max_SE
  LGLindex <- tail(which(abs(PRESS - lowestplus1SE) == min(abs(PRESS - lowestplus1SE)), arr.ind = TRUE), n=1) #in case of multiple best Lasso and Glasso combinations, we choose the last pair
                                                                                                                     # which (usually) is the most sparse result, as the last one has the highest Group Lasso
  OptimumLasso <- LassoSequence[LGLindex[1]]
  OptimumGLasso <- GLassoSequence[LGLindex[2]]
  
  #re-run the analysis with optimum Lasso and Group Lasso
  Pout3dopt <- list()
  LOSSopt <- array()
  
  for (n in 1:NRSTARTS){
    if(method == "datablock"){
      VarSelectResultOpt <- CDfriedmanV1(DATA, Jk, R, OptimumLasso, OptimumGLasso, MaxIter)
    }else if (method == "component"){
      VarSelectResultOpt <- CDfriedmanV2(DATA, Jk, R, OptimumLasso, OptimumGLasso, MaxIter)
    }
    Pout3dopt[[n]] <- VarSelectResultOpt$Pmatrix
    LOSSopt[n] <- VarSelectResultOpt$Loss
  }
  k <- which(LOSSopt == min(LOSSopt))
  if (length(k)>1){
    pos <- sample(1:length(k), 1)
    k <- k[pos]
  }
  PoutFinal <- Pout3dopt[[k]]

  return_crossvali <- list()
  return_crossvali$PoutFinal <- PoutFinal
  return_crossvali$OptimumLasso <- OptimumLasso
  return_crossvali$OptimumGLasso <- OptimumGLasso

  return(return_crossvali)
  
}

############ 3. soft thresholding  ######################

## function soft-thresholding

soft_th <- function(X, lambda){
  
  # assume X is a matrix
  result <- X
  index1 <- which(X > lambda)
  result[index1] <- X[index1] - lambda
  index2 <- which(X < -lambda)
  result[index2] <- X[index2] + lambda
  index0 <- which(X <= lambda & X >= -lambda)
  result[index0] <- 0
  
  return(result)
}

############# 4. rdCV_RSCA  ###############################

rdCV_RSCA <- function(DATA, Jk, R, LassoSequence, GLassoSequence, n_rep, n_seg, MaxIter, NRSTARTS, nfolds){
  
  DATA <- data.matrix(DATA)
  nsub <- dim(DATA)[1]
  
  
  E_hat <- list()
  OptimumLasso <- matrix(NA, n_rep, n_seg)
  OptimumGLasso <- matrix(NA, n_rep, n_seg)
  
  r = 1
  while(r <= n_rep){
    
    cat(sprintf("\n Repetition: %s", r))
    
    randindex <- runif(nsub, 0, 1)
    perc_test <- 1/n_seg
    
    e_hat <- list()
    
    for(i in 1:n_seg){
      
      cat(sprintf("\n Outer loop: %s", i))
      
      testset_index <- ((i - 1) * perc_test < randindex) & (randindex < i * perc_test)
      testset <- DATA[testset_index, ]
      calibset <- DATA[!testset_index, ]
      
      results_innerloop <- cv_sparseSCA(calibset, Jk, R, MaxIter, NRSTARTS, LassoSequence, GLassoSequence, nfolds, method = "component")
      OptimumLasso[r, i] <- results_innerloop$OptimumLasso
      OptimumGLasso[r, i] <- results_innerloop$OptimumGLasso
      estimatedP <- results_innerloop$PoutFinal
      A <- t(estimatedP) %*% t(testset)    
      SVD_DATA <- svd(A, R, R)
      estimatedT <- SVD_DATA$v %*% t(SVD_DATA$u)
      
      e_hat[[i]] <- testset - estimatedT %*% t(estimatedP)  
      
    }
    E_hat[[r]] <- e_hat
    r <- r + 1
  }
  
  results <- list()
  results$Lasso <- OptimumLasso
  results$GLasso <- OptimumGLasso
  results$e_hat <- E_hat
  return(results)

}

######### 5. Tucker Coef   #######################

TuckerCoef <- function(MatrixA, MatrixB){
  
  nrow_data <- dim(MatrixA)[1]
  ncol_data <- dim(MatrixA)[2]
  INDIC_Mat <- gtools::permutations(ncol_data, ncol_data)
  ncol_INDIC <- dim(INDIC_Mat)[1]
  TUCK <- array(NA, dim = c(ncol_INDIC, ncol_data))
  tucker_values <- array()
  tuckerr <- array()
  for(i in 1: ncol_INDIC) {
    MatrixB_perm <- MatrixB[, INDIC_Mat[i,]]
    teller <- 1
    
    for (r in 1: ncol_data){
      vec1 <- MatrixA[, r]
      vec2 <- MatrixB_perm[, r]
      cp <- t(vec1) %*% vec2
      var1 <- t(vec1) %*% vec1
      var2 <- t(vec2) %*% vec2
      
      if (var1 > 0 & var2 > 0){
        tuckerr[teller] <- psych::tr(cp)/sqrt(psych::tr(var1)*psych::tr(var2))
        teller <- teller + 1
      } else if (var2 == 0){
        tuckerr[teller] <- 0
        teller <- teller + 1
      }
    }
    
    tucker_values[i] <- mean(abs(tuckerr))
    TUCK[i,] <- tuckerr
  }
  
  k <- which(tucker_values == max(tucker_values))
  k <- k[1]
  
  perm <- INDIC_Mat[k,]
  tucker_value <- max(tucker_values)
  tucker_vector <- TUCK[k, ]
  
  return_tucker <- list()
  return_tucker$perm <- perm
  return_tucker$tucker_value <- tucker_value
  return_tucker$tucker_vector <- tucker_vector
  return(return_tucker)
}



