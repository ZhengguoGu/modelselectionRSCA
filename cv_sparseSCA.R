cv_sparseSCA <- function(DATA, Jk, R, MaxIter, NRSTARTS, LassoSequence, GLassoSequence, nfolds, method){
  
  DATA <- data.matrix(DATA)
  plotlog <- 0
  
  if(missing(LassoSequence) | missing(GLassoSequence)){
    
    results <- maxLGlasso(DATA, Jk, R)
    GLassomax <- results$Glasso
    Lassomax <- results$Lasso
    
    if(missing(LassoSequence) & missing(GLassoSequence)){
      LassoSequence <- exp(seq(from = log(0.00000001), to = log(Lassomax), length.out = 20))
      GLassoSequence <- seq(from = 0.00000001, to = GLassomax, length.out = 10)  #note that Glasso is not on the log scale, because it is not helpful.
      plotlog <- 1
    } else if(missing(LassoSequence) & (length(GLassoSequence) == 1)){
      LassoSequence <- exp(seq(from = log(0.00000001), to = log(Lassomax), length.out = 50))
      plotlog <- 1
    } else if(missing(GLassoSequence) & (length(LassoSequence) == 1)){
      GLassoSequence <- seq(from = 0.00000001, to = GLassomax, length.out = 50)
      plotlog <- 1
    }
    
  }
  
  if (min(GLassoSequence) < 0) {
    stop("Group Lasso tuning parameter must be non-negative!")
  }
  
  if (min(LassoSequence) < 0) {
    stop("Lasso tuning parameter must be non-negative!")
  }
  if(missing(MaxIter)){
    MaxIter <- 400
  }
  
  if(missing(NRSTARTS)){
    NRSTARTS <- 1
  }
  
  if(missing(nfolds)){
    nfolds <- 10
  }
  if (nfolds < 2){
    stop("Must be at least 2 folds!")
  }
  
  if(missing(method)){
    method <- "component"
  }
  
  
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
      
      cat(sprintf("\nThe cross-validation procedure might take a while to finish. Please be patient."))
      cat(sprintf("\nGroup Lasso: %s", GLassoSequence[g]))
      cat(sprintf("\nLasso: %s", LassoSequence[l]))
      
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
  upper <- vec_PRESS + vec_se
  lower <- vec_PRESS - vec_se 
  #lasso_index0 <- rep(1:length(LassoSequence), length(GLassoSequence))
  #Glasso_index0 <- rep(1:length(GLassoSequence), each=length(LassoSequence))
  
  #lasso_index <- paste("L", lasso_index0)
  #Glasso_index<- factor(paste("G", Glasso_index0), levels=paste("G", 1:length(GLassoSequence)))
  
  
  
  if (length(LassoSequence)>=2 & length(GLassoSequence)>=2){ #### CASE1: multiple lasso and glasso
    
    lasso_index0 <- rep(LassoSequence, length(GLassoSequence))
    Glasso_index0 <- rep(1:length(GLassoSequence), each=length(LassoSequence))
    Glasso_index0 <- factor(paste("G", round(Glasso_index0, digits = 3)), levels=paste("G", 1:length(GLassoSequence)))
    
    lowestPress <- min(vec_PRESS)
    
    if(length(which(vec_PRESS == lowestPress))>1){
      #this happens when min(vec_PRESS) contains multiple values
      lowestplus1SE <- lowestPress + min(vec_se[which(vec_PRESS == lowestPress)]) 
    } else{
      lowestplus1SE <- lowestPress + vec_se[which(vec_PRESS == lowestPress)] #plot 1SE rule region 
    }
    lasso1 <- array()
    lasso2 <- array()
    for(i in 1:length(GLassoSequence)){
      
      pressindex <- which(abs(PRESS[, i]-lowestplus1SE)==min(abs(PRESS[, i]-lowestplus1SE)))  #comment: it seems that we have to have an index number here, instead of inserting the index directly in PRESS[index, i]
      pressindex <- pressindex[length(pressindex)]  #this is because in case of large Glasso values, preindex is a vector, we choose the one with the most sparse results
      lasso1temp <- LassoSequence[pressindex]
      
      if(PRESS[pressindex, i] - lowestplus1SE > 0 ){
        if(LassoSequence[pressindex] == LassoSequence[1]){
          lasso2temp <- lasso1temp  #otherwise lasso2 is out of the boundary
        } else{
          lasso2temp <- LassoSequence[pressindex - 1]
          if(PRESS[pressindex - 1, i]  - lowestplus1SE > 0){
            lasso2temp <- lasso1temp
          }
        }
        #the following condition concerns a rare case 
        
        lasso1[i] <- lasso2temp
        lasso2[i] <- lasso1temp
      }  else if (PRESS[pressindex, i] - lowestplus1SE < 0 ){
        if(LassoSequence[pressindex] == LassoSequence[length(LassoSequence)]){
          lasso2temp <- lasso1temp #otherwise lasso2 is out of the boundary
        } else{
          lasso2temp <- LassoSequence[pressindex + 1]
          #the following condition concerns a rare case 
          if(PRESS[pressindex + 1, i]  - lowestplus1SE < 0) {
            lasso2temp <- lasso1temp
          }
        }
        
        lasso1[i] <- lasso1temp
        lasso2[i] <- lasso2temp
        
      } else { #this is when a PRESS value lies exactly on the 1SE dotted line 
        
        lasso1[i] <- lasso1temp
        lasso2[i] <- lasso1temp
      }
    }
    
    lambdaregion <- cbind(lasso1, lasso2)
    l1matrix <- t(cbind(lasso1, matrix(NA, length(lasso1), length(LassoSequence)-1)))
    l2matrix <- t(cbind(lasso2, matrix(NA, length(lasso2), length(LassoSequence)-1)))
    
    if(plotlog == 1){
      df <- data.frame(GLassoI = Glasso_index0, LassoI = log(lasso_index0), Press = vec_PRESS, Upper = upper, Lower = lower, l1s = c(log(l1matrix)), l2s = c(log(l2matrix)))
      df2 <- data.frame(GLassoI = Glasso_index0, LassoI = log(lasso_index0),  Var = vec_varsel, l1s = c(log(l1matrix)), l2s = c(log(l2matrix)))
      xtag <- "log(Lasso)"
    } else{
      df <- data.frame(GLassoI = Glasso_index0, LassoI = lasso_index0, Press = vec_PRESS, Upper = upper, Lower = lower, l1s = c(l1matrix), l2s = c(l2matrix))
      df2 <- data.frame(GLassoI = Glasso_index0, LassoI = lasso_index0,  Var = vec_varsel, l1s = c(l1matrix), l2s = c(l2matrix))
      xtag <- "Lasso"
    }
    
    
    p1 <- ggplot2::ggplot(df, ggplot2::aes(x=LassoI,y=Press,group=GLassoI)) +
      ggplot2::facet_grid(.~GLassoI)+
      ggplot2::geom_errorbar(ggplot2::aes(ymin=Lower,ymax=Upper, group=GLassoI), width=.1) +
      ggplot2::geom_point(ggplot2::aes(x=LassoI,y=Press,group=GLassoI)) +
      ggplot2::geom_hline(yintercept = upper[which(vec_PRESS == min(vec_PRESS))], linetype = 3) +
      ggplot2::geom_vline(data = subset(df, !is.na(l1s)), ggplot2::aes(xintercept = l1s), linetype = "longdash", col = "red" ) +
      ggplot2::geom_vline(data = subset(df, !is.na(l2s)), ggplot2::aes(xintercept = l2s), linetype = "longdash", col = "red" )      
    
    p1 <- p1 + ggplot2::labs(x = xtag, y="Predicted Mean Squared Errors +/- 1SE")
    
    p2 <- ggplot2::ggplot(df2, ggplot2::aes(x=LassoI,y=vec_varsel,group=GLassoI)) +
      ggplot2::facet_grid(.~GLassoI)+
      ggplot2::geom_point(ggplot2::aes(x=LassoI,y=vec_varsel,group=GLassoI)) +
      ggplot2::geom_vline(data = subset(df, !is.na(l1s)), ggplot2::aes(xintercept = l1s), linetype = "longdash", col = "red" ) +
      ggplot2::geom_vline(data = subset(df, !is.na(l2s)), ggplot2::aes(xintercept = l2s), linetype = "longdash", col = "red" )  
    p2 <- p2 + ggplot2::labs(x = xtag, y="Variables selected in P matrix")
    
    p <- list()
    p[[1]] <- p1
    p[[2]] <- p2
    
  } else if(length(LassoSequence)>=2 & length(GLassoSequence)==1){ #### CASE2: multiple lasso, one glasso
    
    df <- data.frame(LassoI = LassoSequence, Press = vec_PRESS, Upper = upper, Lower = lower)
    
    lowestPress <- min(vec_PRESS)
    lowestplus1SE <- lowestPress + vec_se[which(vec_PRESS == lowestPress)] #plot 1SE rule region 
    lasso1 <- df$LassoI[which(abs(vec_PRESS-lowestplus1SE)==min(abs(vec_PRESS-lowestplus1SE)))]
    if(vec_PRESS[which(abs(vec_PRESS-lowestplus1SE)==min(abs(vec_PRESS-lowestplus1SE)))] - lowestplus1SE > 0 ){
      if(df$LassoI[which(abs(vec_PRESS-lowestplus1SE)==min(abs(vec_PRESS-lowestplus1SE)))] == df$LassoI[1]){
        lasso2 <- lasso1 #otherwise lasso2 is out of the boundary
      } else{
        lasso2 <- df$LassoI[which(abs(vec_PRESS-lowestplus1SE)==min(abs(vec_PRESS-lowestplus1SE))) - 1]
      }
      #the following condition concerns a rare case 
      if((vec_PRESS[which(abs(vec_PRESS-lowestplus1SE)==min(abs(vec_PRESS-lowestplus1SE))) - 1]  - lowestplus1SE > 0) &
         (df$LassoI[which(abs(vec_PRESS-lowestplus1SE)==min(abs(vec_PRESS-lowestplus1SE)))] != df$LassoI[1])){
        lasso2 <- lasso1
        cat("\nWarning! The region for proper tuning parameter values is not available. A possible value is ", lasso1)
        lambdaregion <- glasso1
      } else{
        lambdaregion <- c(lasso2, lasso1)
      }
    } else if (vec_PRESS[which(abs(vec_PRESS-lowestplus1SE)==min(abs(vec_PRESS-lowestplus1SE)))] - lowestplus1SE < 0 ){
      if(df$LassoI[which(abs(vec_PRESS-lowestplus1SE)==min(abs(vec_PRESS-lowestplus1SE)))] == df$LassoI[length(LassoSequence)]){
        lasso2 <- lasso1 #otherwise lasso2 is out of the boundary
      } else{
        lasso2 <- df$LassoI[which(abs(vec_PRESS-lowestplus1SE)==min(abs(vec_PRESS-lowestplus1SE))) + 1]
      }
      
      #the following condition concerns a rare case 
      if((vec_PRESS[which(abs(vec_PRESS-lowestplus1SE)==min(abs(vec_PRESS-lowestplus1SE))) + 1]  - lowestplus1SE < 0) & 
         (df$LassoI[which(abs(vec_PRESS-lowestplus1SE)==min(abs(vec_PRESS-lowestplus1SE)))] != df$LassoI[length(LassoSequence)])){
        lasso2 <- lasso1
        cat("\nWarning! The region for proper tuning parameter values is not available. A possible value is ", lasso1)
        lambdaregion <- lasso1
      } else{
        lambdaregion <- c(lasso1, lasso2)
      }
      
    } else {#this is when a PRESS value lies exactly on the 1SE dotted line
      lasso2 <- lasso1
      lambdaregion <- lasso1
    }
    
    if (plotlog == 1){
      df$LassoI = log(LassoSequence)
      p <- ggplot2::ggplot(df, ggplot2::aes(x=LassoI,y=Press)) +
        ggplot2::geom_errorbar(ggplot2::aes(ymin=Lower,ymax=Upper), width=.1) +
        ggplot2::geom_point(ggplot2::aes(x=LassoI,y=Press))+
        ggplot2::geom_hline(yintercept = upper[which(vec_PRESS == min(vec_PRESS))], linetype = 3)+
        #ggplot2::scale_x_discrete(limits=lasso_index[1:length(LassoSequence)]) +
        ggplot2::geom_vline(xintercept = log(lasso1), 
                            linetype = "longdash", col = "red") +
        ggplot2::geom_vline(xintercept = log(lasso2), 
                            linetype = "longdash", col = "red") 
      p <- p + ggplot2::labs(x = "Lasso (on log scale)", y="Predicted Mean Squared Errors +/- 1SE")
    } else{
      p <- ggplot2::ggplot(df, ggplot2::aes(x=LassoI,y=Press)) +
        ggplot2::geom_errorbar(ggplot2::aes(ymin=Lower,ymax=Upper), width=.1) +
        ggplot2::geom_point(ggplot2::aes(x=LassoI,y=Press))+
        ggplot2::geom_hline(yintercept = upper[which(vec_PRESS == min(vec_PRESS))], linetype = 3)+
        #ggplot2::scale_x_discrete(limits=lasso_index[1:length(LassoSequence)]) +
        ggplot2::geom_vline(xintercept = lasso1, 
                            linetype = "longdash", col = "red") +
        ggplot2::geom_vline(xintercept = lasso2, 
                            linetype = "longdash", col = "red") 
      p <- p + ggplot2::labs(x = "Lasso", y="Predicted Mean Squared Errors +/- 1SE")
    }
  } else if(length(LassoSequence)==1 & length(GLassoSequence)>= 2){####CASE 3: Multiple glasso, one lasso
    
    df <- data.frame(GLassoI = GLassoSequence, Press = vec_PRESS, Upper = upper, Lower = lower)
    
    lowestPress <- min(vec_PRESS)
    lowestplus1SE <- lowestPress + vec_se[which(vec_PRESS == lowestPress)] #plot 1SE rule region 
    glasso1 <- df$GLassoI[which(abs(vec_PRESS-lowestplus1SE)==min(abs(vec_PRESS-lowestplus1SE)))]
    if(vec_PRESS[which(abs(vec_PRESS-lowestplus1SE)==min(abs(vec_PRESS-lowestplus1SE)))] - lowestplus1SE > 0 ){
      if(df$GLassoI[which(abs(vec_PRESS-lowestplus1SE)==min(abs(vec_PRESS-lowestplus1SE)))] == df$GLassoI[1]){
        glasso2 <- glasso1  #otherwise Glasso2 is out of the boundary
      } else{
        glasso2 <- df$GLassoI[which(abs(vec_PRESS-lowestplus1SE)==min(abs(vec_PRESS-lowestplus1SE))) - 1]
      }
      #the following condition concerns a rare case 
      if((vec_PRESS[which(abs(vec_PRESS-lowestplus1SE)==min(abs(vec_PRESS-lowestplus1SE))) - 1]  - lowestplus1SE > 0) &
         (df$GLassoI[which(abs(vec_PRESS-lowestplus1SE)==min(abs(vec_PRESS-lowestplus1SE)))] != df$GLassoI[1])){
        glasso2 <- glasso1
        cat("\nWarning! The region for proper tuning parameter values is not available. A possible value is ", glasso1)
        lambdaregion <- glasso1
      } else{
        lambdaregion <- c(glasso2, glasso1)
      }
      
    } else if (vec_PRESS[which(abs(vec_PRESS-lowestplus1SE)==min(abs(vec_PRESS-lowestplus1SE)))] - lowestplus1SE < 0 ){
      if(df$GLassoI[which(abs(vec_PRESS-lowestplus1SE)==min(abs(vec_PRESS-lowestplus1SE)))] == df$GLassoI[length(GLassoSequence)]){
        glasso2 <- glasso1 #otherwise glasso2 is out of the boundary
      } else{
        glasso2 <- df$GLassoI[which(abs(vec_PRESS-lowestplus1SE)==min(abs(vec_PRESS-lowestplus1SE))) + 1]
      }
      
      #the following condition concerns a rare case 
      if((vec_PRESS[which(abs(vec_PRESS-lowestplus1SE)==min(abs(vec_PRESS-lowestplus1SE))) + 1]  - lowestplus1SE < 0) & 
         (df$GLassoI[which(abs(vec_PRESS-lowestplus1SE)==min(abs(vec_PRESS-lowestplus1SE)))] != df$GLassoI[length(GLassoSequence)])){
        glasso2 <- glasso1
        cat("\nWarning! The region for proper tuning parameter values is not available. A possible value is ", glasso1)
        lambdaregion <- glasso1
      } else{
        lambdaregion <- c(glasso1, glasso2)
      }
      
    } else { #this is when a PRESS value lies exactly on the 1SE dotted line 
      glasso2 <- glasso1
      lambdaregion <- glasso1
    }
    
    if (plotlog == 1){
      df$GLassoI = log(GLassoSequence)
      p <- ggplot2::ggplot(df, ggplot2::aes(x=GLassoI,y=Press)) +
        ggplot2::geom_errorbar(ggplot2::aes(ymin=Lower,ymax=Upper), width=.1) +
        ggplot2::geom_point(ggplot2::aes(x=GLassoI,y=Press))+
        ggplot2::geom_hline(yintercept = upper[which(vec_PRESS == min(vec_PRESS))], linetype = 3) +
        #ggplot2::scale_x_discrete(limits=Glasso_index[1:length(GLassoSequence)])
        ggplot2::geom_vline(xintercept = log(glasso1), 
                            linetype = "longdash", col = "red") +
        ggplot2::geom_vline(xintercept = log(glasso2), 
                            linetype = "longdash", col = "red") 
      p <- p + ggplot2::labs(x = "Group Lasso (on log scale)", y="Predicted Mean Squared Errors +/- 1SE")
    } else{
      p <- ggplot2::ggplot(df, ggplot2::aes(x=GLassoI,y=Press)) +
        ggplot2::geom_errorbar(ggplot2::aes(ymin=Lower,ymax=Upper), width=.1) +
        ggplot2::geom_point(ggplot2::aes(x=GLassoI,y=Press))+
        ggplot2::geom_hline(yintercept = upper[which(vec_PRESS == min(vec_PRESS))], linetype = 3) +
        #ggplot2::scale_x_discrete(limits=Glasso_index[1:length(GLassoSequence)])
        ggplot2::geom_vline(xintercept = glasso1, 
                            linetype = "longdash", col = "red") +
        ggplot2::geom_vline(xintercept = glasso2, 
                            linetype = "longdash", col = "red") 
      p <- p + ggplot2::labs(x = "Group Lasso", y="Predicted Mean Squared Errors +/- 1SE")
    }
    
  }
  
  colnames(lambdaregion) <- c("lower bound", "upper bound")
  return_crossvali <- list()
  return_crossvali$PRESS <- PRESS
  return_crossvali$Press1SE <- lowestplus1SE
  return_crossvali$plot <- p
  return_crossvali$Lasso_values <- LassoSequence
  return_crossvali$Glasso_values <- GLassoSequence
  return_crossvali$Lambdaregion <- lambdaregion
  return(return_crossvali)
  
}
