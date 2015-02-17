cpcq.test<-function(covmats,nvec,B,q)
{
  # covmats: array of covariance matrices to be tested for CPC vs. complete heterogeneity
  # nvec: vector of sample sizes of the k groups
  # B: modal matrix with the q common eigenvectors in the first q columns
  # q: number of common eigenvectors (in B)
  
  k<-dim(covmats)[3]
  p<-dim(covmats)[1]
  
  B1<-B[,c(1:q)]
  B2<-B[,c((q+1):p)]
  
  covmats.cpcq<-array(NA,dim=c(p,p,k))
  for(i in 1:k){
    Q<-eigen(t(B2)%*%covmats[,,i]%*%B2)$vectors
    Bi<-cbind(B1,B2%*%Q)
    Li<-diag(diag(t(Bi)%*%covmats[,,i]%*%Bi))
    covmats.cpcq[,,i]<-Bi%*%Li%*%t(Bi)
  }
  
  chi2total<-0
  
  for(i in 1:k){
    chi2total<-chi2total + (nvec[i]-1)*log(det(covmats.cpcq[,,i])/det(covmats[,,i]))
  }
  
  df<-k*(0.5*p*(p-1)+p) - (0.5*p*(p-1) + k*p + 0.5*(k-1)*(p-q)*(p-q-1))
  
  return(list(chi.square=chi2total,df=df,covmats.cpcq=covmats.cpcq))
  
}
