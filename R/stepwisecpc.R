stepwisecpc <- function(covmats, nvec) {
  
  # Stepwise CPC as described in the paper by N. Trendafilov (2010)
  
  p <- dim(covmats)[1]
  k <- length(nvec)
  ntot <- sum(nvec)
  
  # Calculate pooled covariance matrix
  
  Spooled <- matrix(0, nrow = p, ncol = p)
  for (j in 1:k){
    Spooled <- Spooled + (nvec[j] - 1)*covmats[, , j]
  }
  Spooled <- Spooled / (ntot - k)
  
  # Initial values for stepwise CPC algorithm
  
  Qmat <- matrix(0, nrow = p, ncol = p)
  Qtilde <- eigen(Spooled, symmetric = TRUE)$vectors
  pimat <- diag(p)
  muvec <- rep(0, k)
  elmax <- 10	# maximum number of iterations
  
  # Calculate stepwise CPCs
  
  for (j in 1:p) {
    xvec <- as.vector(Qtilde[, j])
    xvec <- pimat%*%xvec
    
    for (i in 1:k) {
      muvec[i] <- t(xvec)%*%covmats[, , i]%*%xvec
    }
    
    for (el in 1:elmax) {
      Smat <- matrix(0, nrow = p, ncol = p)
      for (i in 1:k) {
        Smat <- Smat + (nvec[i] - 1)*covmats[, , i] / muvec[i]
      }
      yvec <- pimat%*%Smat%*%xvec
      xvec <- yvec / as.numeric(sqrt(t(yvec)%*%yvec))
      for (i in 1:k) {
        muvec[i] <- t(xvec)%*%covmats[, , i]%*%xvec
      }
    }
    
    Qmat[, j] <- xvec
    pimat <- pimat - xvec%*%t(xvec)
    
  }
  
  eigenvals <- matrix(0, p, k)
  for (i in 1:k) {
    eigenvals[, i] <- t((diag(t(Qmat)%*%covmats[, , i]%*%Qmat)))
  }
  
  results <- list(B = Qmat, eigenvals = eigenvals)
  return(results)
}