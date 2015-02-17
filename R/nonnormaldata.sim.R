nonnormaldata.sim<-function(Sigma,n=100,df=rep(2,times=ncol(Sigma))){
  # Generates multivariate nonnormal data with the specified covariance structures (in Sigma)
  
  # Sigma: covariance structure for which to simulate the data
  # n: sample size
  # df: vector of chi square degrees of freedom - to control skewness of the variables (skew = sqrt(8/df))
  
  p<-ncol(Sigma)[1]
  
  xdata<-matrix(NA,nrow=n,ncol=p)
  zdata<-matrix(NA,nrow=n,ncol=p)
  for(j in 1:p){
    xdata[,j]<-rchisq(n=n,df=df[j])
    zdata[,j]<-xdata[,j]/sqrt(2*df[j])
  }
  
  wdata<-matrix(NA,nrow=n,ncol=p)
  Sigma.eigen<-eigen(Sigma)
  sqrt.Sigma<-Sigma.eigen$vectors%*%diag(Sigma.eigen$values^0.5)%*%solve(Sigma.eigen$vectors)
  wdata<-zdata%*%sqrt.Sigma
  
  return(wdata)
}
