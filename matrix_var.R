matrix_var <- function(X){
  
  X <- data.matrix(X)
  Mvar <- sum(X^2)
  
  return(Mvar)
}
