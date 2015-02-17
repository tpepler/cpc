offdiag.vec<-function(datamat)
{
  # Stacks the rows of a square matrix (excluding diagonal elements) in a vector
  
  p<-ncol(datamat)
  datavec<-NULL
  for(j in 1:p){
    datavec<-c(datavec,datamat[j,(1:p)[-j]])
  }
  return(datavec)
}
