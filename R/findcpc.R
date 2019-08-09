findcpc <- function(covmats, B = NULL, cutoff = 0.95, plotting = TRUE, main = "Vector correlations for the permutations")
{
  # Identifies possibly common eigenvectors in k data groups by investigating the vectors correlations of all pairwise combinations of eigenvectors
  
  # covmats: array of covariance matrices from the groups to compare (made by command such as covmats<-array(NA,dim=c(row,col,nummatrices)))
  # B: modal matrix (orthogonal p x p matrix diagonalising the k covariance matrices simultaneously)
  # plotting: should a plot of the dot products be given?
  # cutoff: cutoff point for indicating significance of the geometric mean of the dot products for all pairwise comparisons of a combination of eigenvectors
  
  if(is.null(B)){k <- dim(covmats)[3]}
  else{k <- dim(covmats)[3] + 1}
  p <- dim(covmats)[1]
  
  PCA.array <- array(NA, dim = c(p, p, k))
  
  if(is.null(B)){
    for(i in 1:k){
      PCA.array[, , i] <- eigen(covmats[, , i])$vectors
    }
  }
  else{
    PCA.array[, , 1] <- B
    for(i in 2:k){
      PCA.array[, , i] <- eigen(covmats[, , (i - 1)])$vectors
    }
  }
  
  # Calculating the sum of the vector correlations of all pairwise vector comparisons per vector combination
  
  permsmat <- gtools::permutations(p, k, repeats.allowed = TRUE)  # Sets up matrix with all possible combinations of the columns of the k sets of eigenvectors
  numperms <- nrow(permsmat)  # Total number of combinations to test for commonnality of the vectors
  numcomparisons <- ncol(permsmat) - 1	# Total number of pairwise comparisons to be made per vector combination
  dotvec <- rep(0, times = numperms)
  
  for(i in 1:numperms){
    temp <- rep(NA, times = numcomparisons)
    for(j in 1:numcomparisons){
      testvecs <- cbind(PCA.array[, (permsmat[i, j]), j], PCA.array[, (permsmat[i, j + 1]), j + 1])
      temp[j] <- abs(t(testvecs[, 1]) %*% testvecs[, 2])
    }
    dotvec[i] <- exp(mean(log(temp)))		# Calculate the geometric mean of the vector correlations of all pairwise comparisons for the ith combination of eigenvectors
  }
  
  resultmat <- cbind(permsmat, dot.products = dotvec)
  resultmat <- resultmat[order(resultmat[, "dot.products"], decreasing = TRUE), ]
  
  dot.sig <- resultmat[resultmat[, "dot.products"] >= cutoff, , drop = FALSE]
  
  plotnum <- max(c((p + 2), 20))
  plotnum <- min(c(plotnum, numperms))
  resultmat <- resultmat[c(1:plotnum), ]
  if(plotting){
    resultmat <- resultmat[order(resultmat[, "dot.products"], decreasing = FALSE), ]
    
    par(pch = 20, xaxt = "n")
    plot(x = c(1:plotnum), y = resultmat[, "dot.products"], type = "b", main = main, ylab = "Absolute vector correlation", xlab = "", ylim = c(0, 1))
    abline(h = 1, lty = 3)
    par(adj = 0, srt = 90)
    for(i in 1:plotnum){
      if(resultmat[i, "dot.products"] >= cutoff){
        combtext <- toString(resultmat[i, 1:k])
        points(x = i, y = resultmat[i, "dot.products"], col = "red")
        text(x = i, y = resultmat[i, "dot.products"], labels = combtext, pos = 1, cex = 0.7, offset = 1)
      }
    }
    resultmat <- resultmat[order(resultmat[, "dot.products"], decreasing = TRUE), ]
  }
  
  if(is.null(B)){
    colnames(resultmat) <- c(paste("Group", 1:k, sep=""), "correlations")
  }
  else{
    colnames(resultmat) <- c("B", paste("Group", 1:(k - 1), sep = ""), "correlations")
  }
  
  if(!is.null(B)){
    commonvec.order <- resultmat[1:p, 1]
    commonvec.order <- append(unique(commonvec.order), sort(c(1:p)[-unique(commonvec.order)]))
  }
  else{commonvec.order <- NULL}
  
  par(adj = 0.5)
  #return(list(all.correlations=resultmat,significant=dot.sig,commonvec.order=commonvec.order))
  return(list(all.correlations = resultmat, commonvec.order = commonvec.order))
}
