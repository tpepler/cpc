bootveccor <- function(origdata, veccormat, nvec, reps = 1000)
{
  ## Function to calculate the bootstrap vector correlations
  
  # ONLY ABLE TO HANDLE 2 GROUPS AT THIS STAGE!!!
  
  # origdata: list of the grouped sample data
  # veccormat: matrix of the p eigenvector combinations with the largest dot products
  # nvec: vector of group sizes
  # reps: number of bootstrap replications
  
  numcomb <- nrow(veccormat)
  p <- ncol(origdata[[1]])  # number of variables
  k <- length(nvec)	# number of groups
  bootreps <- matrix(NA, ncol = numcomb, nrow = reps)
  for(r in 1:reps){
    group.PCA <- array(NA, dim = c(p, p, k))
    for(i in 1:k){
      bootdata <- origdata[[i]][sample(c(1:nvec[i]), size = nvec[i], replace = TRUE),]
      group.PCA[, , i] <- eigen(cov(bootdata))$vectors
    }
    for(ci in 1:numcomb){
      bootreps[r, ci] <- abs(t(group.PCA[, veccormat[ci, 1], 1]) %*% group.PCA[, veccormat[ci, 2], 2])
    }
  }
  return(bootreps)
}