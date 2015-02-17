betasim<-function(k, p, q)
{
  # Finds BETA matrices for the required number of groups (k), variables (p) and common principal components (q) for data simulation
  
  # k: number of groups
  # p: number of variables
  # q: number of common principal components, as in CPC(q)
  
  dotprod.cutoff<-cos(((90-2*k)*pi)/(180*k))  # maximum value for dot product of eigenvectors that are NOT common (0.766 corresponds to 40 degree angle)
  dotprod<-1
  BETA.matrices<-array(NA,dim=c(p,p,k))
  library(gtools)
  
  while(dotprod>dotprod.cutoff){
    mat<-matrix(c(runif(n=p*q,min=1,max=100)),nrow=p,ncol=q)
    for(i in 1:k){
      tempmat<-cbind(mat,matrix(c(runif(n=p*(p-q),min=1,max=100)),nrow=p,ncol=(p-q)))
      BETA.matrices[,,i]<-qr.Q(qr(tempmat))
    }
    if(q>=(p-1)){break}
    BETA.test<-BETA.matrices[,(q+1):p,]
    permsmat<-permutations(n=p-q,r=2,repeats.allowed=T)
    numperms<-nrow(permsmat)
    groupcombsmat<-combinations(n=k,r=2,repeats.allowed=F)
    numgroupcombs<-nrow(groupcombsmat)
    
    dotvec<-rep(NA,times=numperms*numgroupcombs)
    for(g in 1:numgroupcombs){
      testmats<-BETA.test[,,groupcombsmat[g,]]
      for(i in 1:numperms){
        dotvec[(g-1)*numperms+i]<-abs(testmats[,permsmat[i,1],1]%*%testmats[,permsmat[i,2],2])
      }
    }
    dotprod<-max(dotvec)
  }
  return(BETA.matrices)
}
