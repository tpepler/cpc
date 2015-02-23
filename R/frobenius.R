frobenius <- function(datamat, targetmat)
{
  # Modified version of the Frobenius measure for a square symmetric matrix with p(p+1)/2 estimable parameters
  
  # datamat: symmetric square matrix for which to calculate the Frobenius measure
  # targetmat: target matrix of same size as datamat, to compare datamat against
  
  p <- ncol(datamat)
  frobtotal <- 0
  frobmat <- (datamat - targetmat)^2
  for(j in 1:p){
    for(h in j:p){
      frobtotal <- frobtotal + frobmat[j, h]
    }
  }
  return(sqrt(frobtotal))
}
