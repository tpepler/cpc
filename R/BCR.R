BCR <- function(origdata, reps = 1000)
{
  ## Bootstrap confidence regions (BCR) method
  ## Calculates 95% bootstrap confidence regions for eigenvector pairs; if regions overlap, the eigenvectors are common
  
  k <- 2  # for two groups only!
  p <- ncol(origdata[[1]])  # for two groups only!
  
  nvec <- rep(NA, times = k)
  covmats <- array(NA, dim = c(p, p, k))
  E <- array(NA, dim = c(p, p, k))
  for(i in 1:k){
    nvec[i] <- nrow(origdata[[i]])
    covmats[, , i] <- cov(origdata[[i]])
    E[, , i] <- eigen(covmats[, , i])$vectors
  }
  
  B <- cpc::FG(covmats, nvec)$B
  findcpc.out <- findcpc(covmats, B = B, plotting = FALSE)
  commonvecnums <- findcpc.out$all.correlations[1:p, ]
  p.new <- p
  for(i in 2:(k + 1)){
    j <- 2
    while(j <= p.new){
      if(length(unique(commonvecnums[1:j, i])) == length(unique(commonvecnums[1:(j - 1), i]))){
        commonvecnums <- commonvecnums[-j, ]
        p.new <- p.new - 1
      }
      j <- j + 1
    }
  }
  
  commonvec.order <- commonvecnums[, "B"]
  eigenvec.order <- commonvecnums[, c("Group1", "Group2")]
  
  common.ind <- rep(1, times = p.new)
  
  for(j in 1:p.new){
    bootreps1 <- matrix(NA, ncol = p, nrow = reps)
    bootreps2 <- matrix(NA, ncol = p, nrow = reps)
    dotprod1 <- rep(NA, times = reps)
    dotprod2 <- rep(NA, times = reps)
    for(r in 1:reps){
      bootreps1[r, ] <- eigen(cov(origdata[[1]][sample(c(1:nvec[1]), size = nvec[1], replace = TRUE), ]))$vectors[, eigenvec.order[j, 1]]
      bootreps2[r, ] <- eigen(cov(origdata[[2]][sample(c(1:nvec[2]), size = nvec[2], replace = TRUE), ]))$vectors[, eigenvec.order[j, 2]]
      dotprod1[r] <- abs((bootreps1[r, , drop = FALSE])%*%(E[, , 1][, eigenvec.order[j, 1], drop = FALSE]))
      dotprod2[r] <- abs((bootreps2[r, , drop = FALSE])%*%(E[, , 2][, eigenvec.order[j, 2], drop = FALSE]))
    }
    bootreps1.cutoff <- quantile(dotprod1, probs = 0.05)
    bootreps1.trim <- bootreps1[dotprod1 > bootreps1.cutoff, ]
    temp <- abs(bootreps1.trim%*%E[, , 2][, eigenvec.order[j, 2]])    
    between.dotprod <- max(temp)    
    within.dotprod <- quantile(dotprod2, probs = 0.05)    
    if(within.dotprod > between.dotprod){common.ind[j] <- 0}
  }
  return(data.frame("Common eigenvector" = commonvec.order, common.ind = common.ind))
}