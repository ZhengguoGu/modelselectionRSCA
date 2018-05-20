# variance measure of a matrix
# Note: this measure is suggested by Katrijn, need to ask for reference. 

matrix_var <- function(X){
  X <- data.matrix(X)
  n_sub <- nrow(X)
  sum((apply(X, 2, var)) * (n_sub-1)/n_sub)
}