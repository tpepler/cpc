G.algorithm <- function(T.mat, nvec)
{
  # Implementation of the G algorithm as described in Flury 1988 (p 181)
  
  k <- dim(T.mat)[3]	# number of groups
  Q <- matrix(rep(200, 4), nrow = 2, ncol = 2)
  Q.nuut <- diag(2)
  delta <- matrix(NA, nrow = k, ncol = 2)
  
  while(abs(Q[2, 1] - Q.nuut[2, 1]) > 1e-07){
    Q <- Q.nuut
    for(i in 1:k){
      for(j in 1:2){
        delta[i, j] <- t(Q[, j]) %*% T.mat[, , i] %*% Q[, j]
      }
    }
    U <- matrix(rep(0, 4), nrow = 2, ncol = 2)
    for(i in 1:k){
      U <- U + nvec[i] * ((delta[i, 1] - delta[i, 2]) / (delta[i, 1] * delta[i, 2])) * T.mat[, , i]
    }
    Q.nuut <- eigvec(U)
  }
  return(Q.nuut)
}
