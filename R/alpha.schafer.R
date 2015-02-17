alpha.schafer<-function(datamat,B,reps=1000)
{
  # Estimates alpha weighting parameter by the method proposed in Schafer & Strimmer (2005), for "Target D",
  # for improved estimation of population covariance matrix
  
  # datamat: matrix containing sample data for the ith group
  # B: matrix of estimated common (and possibly non-common) eigenvectors
  # reps: number of bootstrap replications to use for estimation of the variances of the off-diagonal elements of the L_i
  
  p<-ncol(datamat)
  numobs<-nrow(datamat)
  S<-cov(datamat)
  L<-t(B)%*%S%*%B
  L.offdiagvals<-NULL
  
  for(i in 1:reps){
    bootdata<-datamat[sample(1:numobs,size=numobs,replace=T),]
    bootdata.cov<-cov(bootdata)
    L.boot<-t(B)%*%bootdata.cov%*%B
    L.offdiagvals<-rbind(L.offdiagvals,offdiag.vec(L.boot)) # Requires function offdiag.vec()!!
  }
  numer<-sum(apply(X=L.offdiagvals,MARGIN=2,FUN=var))
  denom<-sum((offdiag.vec(L))^2)
  return(1-min(numer/denom,1))
  #return(list(L.offdiagvals,L,numer,denom,numer/denom))
}
