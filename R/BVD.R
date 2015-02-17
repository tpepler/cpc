BVD <- function(origdata, reps = 1000)
{
  ## Identifies the number of common eigenvectors using the bootstrap
  ## vector correlation distributions (BVD method)
  
  # Depends on functions: bootveccor, findcpc  
    
  k <- 2 # can handle only two groups at this stage
  p <- ncol(origdata[[1]])
  nvec <- rep(NA, times = k)
  covmats <- array(NA, dim = c(p, p, k))
  for(i in 1:k){
    nvec[i] <- nrow(origdata[[i]])
    covmats[, , i] <- cov(origdata[[i]])
  }
  
  B <- cpc::FG(covmats, nvec)$B
  findcpc.out <- findcpc(covmats, B = B, plotting = FALSE)
  commonvecnums <- findcpc.out$all.correlations[1:p, ]
  for(i in 2:(k + 1)){
    j <- 2
    while(j <= p){
      if(length(unique(commonvecnums[1:j, i])) == length(unique(commonvecnums[1:(j - 1), i]))){
        commonvecnums <- commonvecnums[-j, ]
        p <- p - 1
      }
      j <- j + 1
    }
  }
  
  commonvec.order <- commonvecnums[, "B"]
  
  bootreps <- bootveccor(origdata = origdata, veccormat = commonvecnums[, 1:k], nvec = nvec, reps = reps)  # depends on bootveccor function!
  
  commonvecs <- matrix(0, nrow = 7, ncol = p)
  
  for(j in 1:p){
    # BVD 1a: median > 0.71 AND median +- distance between median and 2.5th percentile >= 1
    llim <- quantile(bootreps[, j], probs = 0.025)
    temp.med <- median(bootreps[, j])
    temp.dist <- temp.med - llim
    temp.upper <- temp.med + temp.dist
    if((temp.med > 0.71) & (temp.upper >= 1)){commonvecs[1, j] <- 1}
    
    # BVD 1b: median > 0.75 AND median +- distance between median and 2.5th percentile >= 1
    if((temp.med > 0.75) & (temp.upper >= 1)){commonvecs[2, j] <- 1}
    
    # BVD 1c: median > 0.8 AND median +- distance between median and 2.5th percentile >= 1
    if((temp.med > 0.8) & (temp.upper >= 1)){commonvecs[3, j] <- 1}
    
    # BVD 2a: median on logarithmic sliding scale between 0.6-0.8 AND mediaan +- median +- distance between median and 2.5th percentile >= 1
    lnmax <- log(170)
    sampsize <- mean(nvec)
    mediancutoff <- min(log(max(sampsize - 30, 1)) / lnmax * 0.2 + 0.6, 0.8)
    if((temp.med > mediancutoff) & (temp.upper >= 1)){commonvecs[4, j] <- 1}
    
    # BVD 2b: median on logarithmic sliding scale between 0.65-0.8 AND mediaan +- median +- distance between median and 2.5th percentile >= 1
    mediancutoff <- min(log(max(sampsize - 30, 1)) / lnmax * 0.15 + 0.65, 0.8)
    if((temp.med > mediancutoff) & (temp.upper >= 1)){commonvecs[5, j] <- 1}
    
    # BVD 2c: median on logarithmic sliding scale between 0.7-0.8 AND mediaan +- median +- distance between median and 2.5th percentile >= 1
    mediancutoff <- min(log(max(sampsize - 30, 1)) / lnmax * 0.1 + 0.7, 0.8)
    if((temp.med > mediancutoff) & (temp.upper >= 1)){commonvecs[6, j] <- 1}
    
    # BVD 2d: median on logarithmic sliding scale between 0.75-0.8 AND mediaan +- median +- distance between median and 2.5th percentile >= 1
    mediancutoff <- min(log(max(sampsize - 30, 1)) / lnmax * 0.05 + 0.75, 0.8)
    if((temp.med > mediancutoff) & (temp.upper >= 1)){commonvecs[7, j] <- 1}
  }
  
  rownames(commonvecs) <- c("BVD 1a", "BVD 1b","BVD 1c",
                            "BVD 2a", "BVD 2b", "BVD 2c", "BVD 2d")
  commonvecs <- rbind("Common eigenvector" = commonvec.order, commonvecs)
  return(commonvecs)
}