prop.test<-function(covmats,nvec)
{
  # covmats: array of covariance matrices to be tested for proportionality vs. complete heterogeneity
  # nvec: vector of sample sizes of the k groups
  
  k<-dim(covmats)[3]
  p<-dim(covmats)[1]
  
  rvec<-(nvec-1)/(sum(nvec)-k)
  
  # Step PCM0
  
  rho<-rep(1,times=k)
  maxrho<-1
  prevmaxrho<-100
  
  while(abs(maxrho-prevmaxrho)>1e-09){
    
    # Step PCM1
    
    covmats.total<-0
    for(i in 1:k){
      covmats.total<-covmats.total + rvec[i]*covmats[,,i]/rho[i]
    }
    B<-eigen(covmats.total)$vectors
    avec<-matrix(NA,ncol=p,nrow=k)
    for(i in 1:k){
      avec[i,]<-diag(t(B)%*%covmats[,,i]%*%B)
    }
    
    # Step PCM2
    
    lambda<-rep(NA,times=p)
    for(j in 1:p){
      lambda[j]<-sum(rvec*avec[,j]/rho)
    }
    
    # Step PCM3
    
    for(i in 2:k){
      rho[i]<-1/p*sum(avec[i,]/lambda)
    }
    
    # Step PCM4
    
    prevmaxrho<-maxrho
    maxrho<-max(rho[2:k])
  }
  
  #  return(rho)
  
  covmats.prop<-array(NA,dim=c(p,p,k))
  for(i in 1:k){
    covmats.prop[,,i]<-B%*%(diag(lambda)*rho[i])%*%t(B)
  }
  
  chi2total<-0
  
  for(i in 1:k){
    chi2total<-chi2total + (nvec[i]-1)*log(det(covmats.prop[,,i])/det(covmats[,,i]))
  }
  
  df<-k*(0.5*p*(p-1)+p) - (0.5*p*(p-1)+p+k-1)
  
  return(list(chi.square=chi2total,df=df,covmats.prop=covmats.prop))
  
}
