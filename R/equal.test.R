equal.test<-function(covmats,nvec)
{
  # covmats: array of covariance matrices to be tested for equality vs. complete heterogeneity
  # nvec: vector of sample sizes of the k groups
  
  k<-dim(covmats)[3]
  p<-dim(covmats)[1]
  
  covmats.pooltotal<-0
  
  for(i in 1:k){
    covmats.pooltotal<-covmats.pooltotal + covmats[,,i]*(nvec[i]-1)
  }
  
  covmats.pool<-covmats.pooltotal/(sum(nvec)-k)
  
  chi2total<-0
  covmats.equal<-array(NA,dim=c(p,p,k))
  
  for(i in 1:k){
    covmats.equal[,,i]<-covmats.pool
    chi2total<-chi2total + (nvec[i]-1)*log(det(covmats.equal[,,i])/det(covmats[,,i]))
  }
  
  df<-k*(0.5*p*(p-1)+p) - (0.5*p*(p-1)+p)
  
  return(list(chi.square=chi2total,df=df,covmats.equal=covmats.equal))
}
