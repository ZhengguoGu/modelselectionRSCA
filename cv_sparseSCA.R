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
    
    for (l in 1:length(LassoSequence)){
      
      if(method == "datablock"){
        Forvarselected <- CDfriedmanV1(DATA, Jk, R, LassoSequence[l], GLassoSequence[g], MaxIter)
      }else if (method == "component"){
        Forvarselected <- CDfriedmanV2(DATA, Jk, R, LassoSequence[l], GLassoSequence[g], MaxIter)
      }
      varselected[l,g] <- sum(Forvarselected$Pmatrix != 0)  #how many variables in P have been selected?
      
      error_x <- array()
      for (i in 1:nfolds){
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
