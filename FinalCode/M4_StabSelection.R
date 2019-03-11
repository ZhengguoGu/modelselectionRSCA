####### parameters ########################################################################
# DATA: concatenated data matrix
# Jk: number of variables (a vector, the elements of which correspond to the datablocks)
# R: number of components
# LassoSequence: A sequence of Lasso tuning parameters
# N_loading: The total number of non-zero loadings in the P matrix
# Thr: threshold for the probability of each loading; the default is .9
# NRSTARTS: Number of random starts
# N_cores: number of CPU logical cores can be used
###########################################################################################

M4_StabSelection <- function(DATA, Jk, R, LassoSequence, N_loading, Thr, NRSTARTS, N_cores){
  
  DATA <- data.matrix(DATA)
  n_persons <- nrow(DATA)
  
  if(missing(LassoSequence)){
    LassoSequence = exp(seq(from = log(0.00000001), to = log(RegularizedSCA::maxLGlasso(DATA, Jk, R)$Lasso), length.out = 100))
  }
   
    if(missing(NRSTARTS)){
    NRSTARTS = 5
  }
  
  if(missing(Thr)){
    Thr = .5
  }
  
  # need a target matrix for T, so that all the estimated T can be compared to it.
  result <- RegularizedSCA::cv_sparseSCA(DATA, Jk, R, MaxIter = 400, NRSTARTS = NRSTARTS, LassoSequence, GLassoSequence=0, nfolds = 7, method = "component")  
  T_target <- result$T_hat                #We fix the estimated T matrix. All the estimated P in the following resampling procedure will be rotated 
                                          #after comparing the estimated T with with T_target.
                                          #Arguably, using CV to generate T_target might not be a good idea. By the time I'm writing this code, I have already
                                          #found out that index of sparseness works better. But CV is the most well known method, so I use CV here.
 
  LassoSequence <- sort(LassoSequence, decreasing = T) # the largest value enters first
  
  P_prob <- list()
  #### this is for P_prob[[1]], which is when used to compare to P_prob[[j]] (j= 2, ...) so as to record the highest probability across all j's (j=1,...)
  cl <- snow::makeCluster(N_cores)
  doSNOW::registerDoSNOW(cl)
  #note that set.seed() and %dorng% ensure that parallel computing generates reproducable results.
  sim_result <- foreach::foreach(i = 1:100, .combine='+') %dorng% {
    
    person_index <- sample(1:n_persons, n_persons/2, replace = F)
    result <- RegularizedSCA::sparseSCA(DATA[person_index, ], Jk, R, LASSO = LassoSequence[1], GROUPLASSO = 0, MaxIter = 400, NRSTARTS = NRSTARTS, method = "component")
    perm <- RegularizedSCA::TuckerCoef(T_target[person_index, ], result$Tmatrix)$perm
    P_result <- result$Pmatrix[, perm]
    P_result[which(P_result!=0)] <- 1
    
    return(P_result)
    
  }
  snow::stopCluster(cl)
  P_prob[[1]] <- data.frame(sim_result)/100
  P_final <- P_prob[[1]]  #P_final will be updated after each comparison (see below)
  j <- 2
  while(j <= length(LassoSequence)){  #note that it is not necessary to run j all the way to the length(LassoSequence), to some point, this algorthm will stop because enough non-zero loadings are gathered.
    cl <- snow::makeCluster(N_cores)
    doSNOW::registerDoSNOW(cl)
    #note that set.seed() and %dorng% ensure that parallel computing generates reproducable results.
    sim_result <- foreach::foreach(i = 1:100, .combine='+') %dorng% {
      
      person_index <- sample(1:n_persons, n_persons/2, replace = F)
      result <- RegularizedSCA::sparseSCA(DATA[person_index, ], Jk, R, LASSO = LassoSequence[j], GROUPLASSO = 0, MaxIter = 400, NRSTARTS = NRSTARTS, method = "component")
      perm <- RegularizedSCA::TuckerCoef(T_target[person_index, ], result$Tmatrix)$perm
      P_result <- result$Pmatrix[, perm]
      P_result[which(P_result!=0)] <- 1
      
      return(P_result)
      
    }
    snow::stopCluster(cl)
    
    P_prob[[j]] <- data.frame(sim_result)/100
    
    index <- which((P_final < P_prob[[j]]), arr.ind = TRUE)  # which((P_final < P_prob[[j]]) == TRUE, arr.ind = TRUE) is also fine
    P_final[index] <- P_prob[[j]][index]  #update the P_final matrix, keep the highest probability
    
    if(sum(P_prob[[j]] >= Thr) > N_loading){
       j <- length(LassoSequence)  # just a way to skip the remaining Lasso values, since the total number of non-zero loadings already > N_loading prespecified by the user
    }
    
    j <- j+1
    print(j)
  }
  
  # Note that in the above code if(sum(P_prob[[j]] >= Thr) > N_loading), it could happen that at j, sum(P_prob[[j]]) >> N_loading. We would want to reduce the size of sum(P_prob[[j]]) close to N_loading
  thr_new <- sort(as.matrix(P_final), decreasing = TRUE)[N_loading] # this is the lowest probability whose corresponding loading should not be zero, ACCORDING TO N_loading
  if(dim(which(P_final== thr_new, arr.ind = TRUE))[1] > 1){
    print("It's likely that the total # of non-0 loadings in the final P_hat greatly exceeds N_loading!")  #this is problematic, in this case, more than one loading corresponds to the lowest probability. 
  }
  
  P_final[which(P_final < thr_new, arr.ind = TRUE)] <- 0
  
  # Reestimate P and T, with 5 starts
  Pout3d <- list()
  Tout3d <- list()
  LOSS <- array()
  LOSSvec <- list()
  
  for (n in 1:NRSTARTS) { 
    VarSelectResult <- StrucSCA_withIndex(DATA, Jk, R, P_indexset = P_final, MaxIter=300)
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
  
  PoutBest <- Pout3d[[k]]  #this is the final, estimated P
  ToutBest <- Tout3d[[k]]  #this is the final, estimated T
  
  return_list <- list(P_hat = PoutBest, T_hat = ToutBest)
  return(return_list)
  
  
}