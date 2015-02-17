alpha.crossvalid<-function(datamat,B,reps=100)
{
  # Estimates alpha weighting parameter by cross-validation, for improved estimation of population covariance matrix
  
  # datamat: matrix containing sample data for the ith group
  # B: matrix of estimated common (and possibly non-common) eigenvectors
  # reps: number of replications to use in cross-validation
  
  #sample1<-datamat[[1]]
  #sample2<-datamat[[2]]
  p<-ncol(datamat)
  numobs<-nrow(datamat)
  train.n<-round(numobs*0.7,0)
  alpha.vals<-seq(from=0,to=1,by=0.01)
  alpha.n<-length(alpha.vals)
  min.error.alpha<-rep(NA,times=reps)
  
  for(i in 1:reps){
    sampledata<-datamat[sample(1:numobs,size=numobs,replace=F),]
    traindata<-sampledata[1:train.n,]
    testdata<-sampledata[(train.n+1):(numobs),]
    traindata.covmat<-cov(traindata)
    testdata.covmat<-cov(testdata)
    #S<-array(NA,dim=c(p,p,2))
    #S[,,1]<-cov(sample1)
    #S[,,2]<-cov(traindata)  
    #B<-FGALG(A=S,n.vek=c(nvek[1],train.n))$B
    #if(q<p){
    #  commonvec.order<-findcpc(S=S,B=B,plotting=F)$commonvec.order
    #  B<-partialcpc.estimate(S=S,nvek=c(nvek[1],train.n),B=B,commonvec.order=commonvec.order,q=q)[,,alpha.group]
    #}
    L.diag<-diag(diag(t(B)%*%traindata.covmat%*%B))
    S.cpc<-B%*%L.diag%*%t(B)
    frobvals<-rep(NA,times=alpha.n)
    for(j in 1:alpha.n){
      S.new<-alpha.vals[j]*traindata.covmat+(1-alpha.vals[j])*S.cpc
      frobvals[j]<-frobenius(S.new,testdata.covmat)
    }
    #return(frobvals)
    min.error.alpha[i]<-alpha.vals[which.min(frobvals)]
  }
  return(mean(min.error.alpha))
}
