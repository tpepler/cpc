flury.test<-function(covmats,nvec,B=cpc::FG(covmats,nvec)$B,p=dim(covmats)[1],qmax=p-2,commonvec.order=findcpc(covmats=covmats,B=B,plotting=F)$commonvec.order)
{
  # Calculates the partial chi square statistics and AIC values for all the models in Flury's (1988) hierarchy
  
  # Depends on: equal.test, prop.test, cpc.test, cpcq.test, AIC, findcpc functions!
  # Depends on: FGALG, GALG, EIGVEC functions!
  
  # covmats: array of covariance matrices to be tested, created by a command such as covmats<-array(NA,dim=c(p,p,k))
  # nvec: vector of sample sizes of the k groups
  # B: modal matrix (orthogonal p x p matrix diagonalising the k covariance matrices simultaneously)
  # commonvec.order: vector indicating the order of the most likely candidates of common eigenvectors - from 1 (most likely) to p (least likely)
  # p: number of variables
  # qmax: maximum for q when estimating the CPC(q) models
  
  if((qmax+2)>p){
    qmax<-p-2
    model.names<-c("Equality","Proportionality","CPC",paste("CPC(",seq(from=qmax,to=1),")",sep=""),"Heterogeneity")
    No.of.CPCs<-c(p,p,p,(p-2):1,0)
  }
  else if(qmax<1){
    qmax<-0
    model.names<-c("Equality","Proportionality","CPC","Heterogeneity")
    No.of.CPCs<-c(p,p,p,0)
  }
  else{
    model.names<-c("Equality","Proportionality","CPC",paste("CPC(",seq(from=qmax,to=1),")",sep=""),"Heterogeneity")
    No.of.CPCs<-c(p,p,p,(p-2):1,0)
  }
  
  nmodels<-length(model.names)
  chi.square<-rep(NA,times=nmodels)
  df<-rep(NA,times=nmodels)
  model.AIC<-rep(NA,times=nmodels)
  
  # Equality
  
  equal.test.output<-equal.test(covmats,nvec)
  chi.square[1]<-equal.test.output$chi.square
  df[1]<-equal.test.output$df
  model.AIC[1]<-flury.AIC(equal.test.output$covmats.equal,covmats,nvec,df=equal.test.output$df)
  
  # Proportionality
  
  prop.test.output<-prop.test(covmats,nvec)
  chi.square[2]<-prop.test.output$chi.square
  df[2]<-prop.test.output$df
  model.AIC[2]<-flury.AIC(prop.test.output$covmats.prop,covmats,nvec,df=prop.test.output$df)
  
  # CPC
  
  cpc.test.output<-cpc.test(covmats=covmats,nvec=nvec,B=B)
  chi.square[3]<-cpc.test.output$chi.square
  df[3]<-cpc.test.output$df
  model.AIC[3]<-flury.AIC(cpc.test.output$covmats.cpc,covmats,nvec,df=cpc.test.output$df)
  
  # CPC(q)
  
  if(qmax>0){
    B<-B[,commonvec.order]
    q<-qmax
    for(i in 1:qmax){
      cpcq.test.output<-cpcq.test(covmats,nvec,B,q=q)
      chi.square[3+i]<-cpcq.test.output$chi.square
      df[3+i]<-cpcq.test.output$df
      model.AIC[3+i]<-flury.AIC(cpcq.test.output$covmats.cpcq,covmats,nvec,df=cpcq.test.output$df)
      q<-q-1
    }
  }
  
  # Heterogeneity
  
  model.AIC[3+qmax+1]<-flury.AIC(covmats,covmats,nvec,df=0)
  
  chi.square[1:(nmodels-2)]<-chi.square[1:(nmodels-2)]-chi.square[2:(nmodels-1)]
  df[1:(nmodels-2)]<-df[1:(nmodels-2)]-df[2:(nmodels-1)]
  chi.div.df<-chi.square/df
  chi.square<-round(chi.square,2)
  chi.div.df<-round(chi.div.df,2)
  model.AIC<-round(model.AIC,2)
  
  resultmat<-data.frame(Model=model.names,Chi.square=chi.square,DF=df,Chi2.div.df=chi.div.df,AIC=model.AIC,No.of.CPCs=No.of.CPCs)
  return(resultmat)
}
