B.partial <- function(covmats, nvec, B = cpc::FG(covmats = covmats, nvec = nvec)$B, commonvec.order, q)
{
  # Estimates matrices of common (and non-common) eigenvectors for k groups
  
  # covmats: array of sample covariance matrices for the k groups
  # nvec: vector of sample sizes of the k groups
  # B: matrix of common eigenvectors (estimated under the assumption of full CPC)
  # commonvec.order: order of the common eigenvectors in B (with the q truly common eigenvectors in the first q positions)
  # q: number of eigenvectors common to all k groups
  
  k <- dim(covmats)[3]
  p <- dim(covmats)[1]
  B <- B[, commonvec.order]
  Bmats <- array(NA, dim = c(p, p, k))
  for(i in 1:k){
    B1 <- B[, 1:q, drop = FALSE]
    B2 <- B[, (q + 1):p]
    Q1 <- eigen(t(B2) %*% covmats[, , i] %*% B2)$vectors
    B21 <- B2 %*% Q1
    Bmats[, , i] <- cbind(B1, B21)
  }
  return(Bmats)
}