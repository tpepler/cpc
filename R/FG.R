FG <- function(covmats, nvec, method = c('ML', 'LS'))
{
  # Implementation of the FG algorithm as described in Flury 1988 (p 178),
  # with Maximum Likelihood estimation. Least Squares estimation method as
  # described in Beaghen (1997) Canonical variate analysis and related
  # methods with longitudinal data (PhD dissertation, Appendix 2).
  
  # covmats: array of covariance matrices to be simultaneously diagonalized,
  # created by a command such as covmats <- array(NA, dim = c(p, p, k))
  # nvec: vector of sample sizes for the covariance matrices in covmats
  # method: estimation method, either Maximum Likelihood ('ML', the default)
  # or Least Squares ('LS')
  
  p <- dim(covmats)[2]
  B <- diag(p)
  DIFF <- 100
  k <- dim(covmats)[3]
  while(DIFF > 1e-09){
    B.old <- B
    T.mat <- array(NA, dim = c(2, 2, k))
    for(m in 1:(p - 1)){                                       
      for(j in (m + 1):p){    	# m<-1; j<-2
        vek <- c(m, j)
        for(i in 1:k){
          T.mat[, , i] <- t(B[, vek]) %*% covmats[, , i] %*% B[, vek]
        }
        J <- G.algorithm(T.mat, nvec, method = method[1])
        B[, vek] <- B[, vek] %*% J
      }
    }
    for(i in 1:p){
      for(j in 1:p){
        DIFF <- abs(B[i, j] - B.old[i, j])
      }
    }
  }
  
  # Order the columns of B
  
  diagvals <- 0
  for(i in 1:k){
    diagvals <- diagvals + diag(t(B) %*% covmats[, , i] %*% B)
  }
  
  B <- B[, order(diagvals, decreasing = TRUE)]
  
  diagvalsmat <- matrix(NA, nrow = p, ncol = k)
  for(i in 1:k){
    diagvalsmat[, i] <- diag(t(B) %*% covmats[, , i] %*% B)
  }
  
  return(list(B = B, diagvals = diagvalsmat))
}