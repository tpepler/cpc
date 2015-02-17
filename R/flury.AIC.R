flury.AIC<-function(covmats.high,covmats.low,nvec,df)
{
  # covmats.high: estimates of the covariance matrices under the "higher" model
  # covmats.low: estimates of the covariance matrices under the "lower" model (usually unrelated/individual covariance matrices)
  # nvec: vector of sample sizes
  # df: degrees of freedom of the higher model, versus unrelated covariance matrices
  
  k<-length(nvec)
  p<-dim(covmats.high)[1]
  aic.total<-0
  for(i in 1:k){
    aic.total<-aic.total + (nvec[i]-1)*(sum(diag(solve(covmats.high[,,i])%*%covmats.low[,,i])) + log(det(covmats.high[,,i])) - sum(diag(solve(covmats.low[,,i])%*%covmats.low[,,i])) - log(det(covmats.low[,,i])))
  }
  
  npar<-k*(0.5*p*(p-1)+p) - (0.5*p*(p-1)+p) - df
  
  aic.criterion<-aic.total+2*npar
  return(aic.criterion)
}
