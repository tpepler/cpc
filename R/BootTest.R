BootTest<-function(origdata,q=ncol(origdata[[1]]),reps=1000)
{
  # Bootstrap test for common eigenvectors as proposed by Klingenberg (1996) -- Bootstrap hypothesis test (BootTest) method
  # H_0: eigenvector pair are common
  # H_1: eigenvector pair are not common
  
  # origdata: list of the original data for the k groups
  # q: number of common components to test for
  # reps: number of bootstrap replications to use
  
  # ONLY ABLE TO HANDLE 2 GROUPS AT THIS STAGE!!!
  
  k<-2  # works only for k = 2 groups at this stage!!
  p<-ncol(origdata[[1]])
  
  nvec<-rep(NA,times=k)
  covmats<-array(NA,dim=c(p,p,k))
  for(i in 1:k){
    nvec[i]<-nrow(origdata[[i]])
    covmats[,,i]<-cov(origdata[[i]])
  }
  
  B<-cpc::FG(covmats,nvec)$B
  findcpc.out<-findcpc(covmats,B=B,plotting=FALSE)
  commonvecnums<-findcpc.out$all.correlations[1:p,]
  
  temp.p<-p
  for(i in 1:(k+1)){
    j<-2
    while(j <= temp.p){
      if(length(unique(commonvecnums[1:j,i]))==length(unique(commonvecnums[1:(j-1),i]))){
        commonvecnums<-commonvecnums[-j,]
        temp.p<-temp.p-1
        j<-j-1
      }
      j<-j+1
    }
  }
  
  commonvec.order<-commonvecnums[,"B"]
  eigenvec.order<-commonvecnums[,c("Group1","Group2")]
  q<-min(q,nrow(eigenvec.order))
  eigenvec.order.full<-matrix(NA,nrow=p,ncol=k)
  if(q<p){
    for(i in 1:k){
      tempvec<-c(1:p)
      tempvec<-tempvec[-eigenvec.order[,i]]
      eigenvec.order.full[,i]<-c(eigenvec.order[,i],tempvec)
    }
  }
  else{eigenvec.order.full<-eigenvec.order}
  
  E<-array(NA,dim=c(p,p,k))
  data.rotated<-origdata
  for(i in 1:k){
    E[,,i]<-eigen(covmats[,,i])$vectors[,eigenvec.order.full[,i]]
    data.rotated[[i]]<-as.matrix(origdata[[i]])%*%E[,,i]%*%t(B)
  }
  orig.dotproducts<-abs(diag(t(E[,,1])%*%E[,,2])[1:q])
  
  
  bootreps<-matrix(NA,ncol=q,nrow=reps)
  for(r in 1:reps){
    rep.eigenvecs<-array(NA,dim=c(p,q,k))
    for(i in 1:k){
      bootdata<-data.rotated[[i]][sample(c(1:nvec[i]),size=nvec[i],replace=T),]
      rep.eigenvecs[,,i]<-eigen(cov(bootdata))$vectors[,eigenvec.order[1:q,i]]
    }
    for(j in 1:q){
      bootreps[r,j]<-abs(rep.eigenvecs[,j,1]%*%rep.eigenvecs[,j,2])
    }
  }
  
  pvals<-rep(NA,times=q)
  for(j in 1:q){
    pvals[j]<-nrow(bootreps[bootreps[,j]<=orig.dotproducts[j],,drop=F])/reps
  }
  
  return(data.frame(eigenvec.order[1:q,c("Group1","Group2"),drop=F],vec.correlations=orig.dotproducts[1:q],p.values=pvals[1:q]))
}
