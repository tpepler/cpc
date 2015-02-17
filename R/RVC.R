RVC<-function(covmats,reps=100000)
{
  # Randomisation test for common eigenvectors as proposed by Klingenberg (1996) -- Random vector correlations (RVC) method
  # H_0: eigenvector pair are not common
  # H_1: eigenvector pair are common
  
  # covmats: array of covariance matrices of the k groups
  # reps: number of randomisations to use
  
  # ONLY ABLE TO HANDLE 2 GROUPS AT THIS STAGE!!!
  
  k<-2  # works only for k = 2 groups at this stage!!
  p<-dim(covmats)[1]
  
  E<-array(NA,dim=c(p,p,k))
  for(i in 1:k){
    E[,,i]<-eigen(covmats[,,i])$vectors
  }
  
  rand.dotproducts<-rep(NA,times=reps)
  for(r in 1:reps){
    randvec1<-runif(n=p,min=-1,max=1)
    randvec1<-randvec1/sqrt(randvec1%*%randvec1)
    randvec2<-runif(n=p,min=-1,max=1)
    randvec2<-randvec2/sqrt(randvec2%*%randvec2)
    rand.dotproducts[r]<-abs(randvec1%*%randvec2)
  }
  
  commonvecnums<-findcpc(covmats,plotting=FALSE)$all.correlations[1:p,]
  for(i in 1:k){
    j<-2
    while(j <= p){
      if(length(unique(commonvecnums[1:j,i]))==length(unique(commonvecnums[1:(j-1),i]))){
        commonvecnums<-commonvecnums[-j,]
        p<-p-1
      }
      j<-j+1
    }
  }
  
  commonvec.order<-commonvecnums[,1:2]
  
  orig.dotproducts<-abs(diag(t(E[,,1][,commonvec.order[,1]])%*%(E[,,2][,commonvec.order[,2]])))
  
  pvals<-rep(NA,times=p)
  for(j in 1:p){
    pvals[j]<-length(rand.dotproducts[rand.dotproducts>=orig.dotproducts[j]])/reps
  }
  
  return(data.frame(commonvec.order,vec.correlations=orig.dotproducts,p.values=pvals))
}
