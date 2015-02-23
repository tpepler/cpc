cpc.test <- function(covmats, nvec, B = cpc::FG(covmats = covmats, nvec = nvec)$B)
{
  # covmats: array of covariance matrices to be tested for CPC vs. complete heterogeneity
  # nvec: vector of sample sizes of the k groups
  # B: modal matrix (orthogonal p x p matrix diagonalising the k covariance matrices simultaneously)
  
  k <- dim(covmats)[3]
  p <- dim(covmats)[1]
  
  covmats.cpc <- array(NA, dim = c(p, p, k))
  chi2total <- 0
  for(i in 1:k){
    lambda <- diag(t(B) %*% covmats[, , i] %*% B)
    covmats.cpc[, , i] <- B %*% diag(lambda) %*% t(B)
    chi2total <- chi2total + (nvec[i] - 1) * log(det(covmats.cpc[, , i]) / det(covmats[, , i]))
  }
  
  df <- k * (0.5 * p * (p - 1) + p) - (0.5 * p * (p - 1) + k * p)
  
  return(list(chi.square = chi2total, df = df, covmats.cpc = covmats.cpc))
}
