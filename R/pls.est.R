pls.est <- function(X, Y)
{
  # Implementation of the partial least squares (PLS) algorithm
  # as described in Elements of Statistical Learning -- Hastie et al. (2009)
  
  # X: independent variables matrix
  # Y: dependent/response variable vector
  
  X.orig <- as.matrix(X)
  X.sdinv <- diag(1 / apply(X.orig, 2, sd))
  X.orig <- t(t(X.orig) - apply(X.orig, 2, mean)) %*% X.sdinv # standardise columns of X.orig to have zero mean and unit variance
  Y <- as.vector(Y)
  
  ## PLS, as from Algorithm 3.3 (p81) in Hastie et al. (2009)
  X <- X.orig
  p <- ncol(X)
  y0 <- mean(Y) * rep(1, times = nrow(X))
  Z <- matrix(0, nrow = nrow(X), ncol = p)
  Ymat <- matrix(NA, nrow = nrow(X), ncol = p)
  for(m in 1:p){
    for(j in 1:p){
      phi <- X[, j] %*% Y
      Z[, m] <- Z[, m] + phi * X[, j]
    }
    theta <- (Z[, m] %*% Y) / (Z[, m] %*% Z[, m])
    if(m == 1){Ymat[, m] <- y0 + theta * Z[, m]}
    else{Ymat[, m] <- Ymat[, m - 1] + theta * Z[, m]}
    for(j in 1:p){
      X[, j] <- X[, j] - ((Z[, m] %*% X[, j]) / (Z[, m] %*% Z[, m])) %*% Z[, m]
    }
  }
  
  B <- solve(t(X.orig) %*% X.orig) %*% t(X.orig) %*% Z
  for(j in 1:p){
    B[, j] <- B[, j] / sqrt(B[, j] %*% B[, j])
  }
  
  return(list(pls.scores = Z, pls.loadings = B))
}
