ensemble.test<-function(origdata,standardize=F)
{
  ## Function to identify common eigenvectors using the ensemble method
  
  # Ensemble method to identify common eigenvectors in k groups: majority vote on number 
  # of common eigenvectors from Flury's AIC, Bootstrap Vector Correlation Distribution (BVD), 
  # Bootstrap Confidence Regions (BCR), Random Vector Correlations (RVC) and Bootstrap hypothesis test (BootTest) methods
  
  # Note: ONLY ABLE TO HANDLE 2 GROUPS AT THIS STAGE!!!
  
  # origdata: list of the sample groups data
  # standardize: should the data be standardized (mean=0, stdev=1)?
  
  myFunc<-function(datavec){
    return((datavec-mean(datavec))/sd(datavec))
  }
  
  standcol<-function(datamat){
    return(apply(datamat,2,myFunc))
  }
  
  nvec<-c(nrow(origdata[[1]]),nrow(origdata[[2]]))
  p<-ncol(origdata[[1]])
  covmats<-array(NA,dim=c(p,p,2))
  if(standardize){
    origdata[[1]]<-standcol(origdata[[1]])
    origdata[[2]]<-standcol(origdata[[2]])
  }
  covmats[,,1]<-cov(origdata[[1]])
  covmats[,,2]<-cov(origdata[[2]])
  
  B<-cpc::FG(covmats,nvec)$B
  commonvec.order<-findcpc(covmats,B=B,plotting=F)$commonvec.order
  
  flury.out<-flury.test(covmats,nvec,B=B,p=dim(covmats)[1],qmax=p-2,commonvec.order=commonvec.order)
  fluryAIC.vote<-flury.out[which.min(flury.out[,"AIC"]),"No.of.CPCs"]
  fluryAIC.ind<-c(rep(1,times=fluryAIC.vote),rep(0,times=(p-fluryAIC.vote)))
  # Flury AIC common eigenvector order: same as findcpc
  
  BVD.out<-BVD(origdata,reps=1000)
  #BVD2c.ind<-BVD.out["BVD 2c",]
  #if(length(BVD2c.ind)<p){BVD2c.ind<-c(BVD2c.ind,rep(0,times=(p-length(BVD2c.ind))))}
  BVD.ind<-BVD.out["BVD 1a",]
  if(length(BVD.ind)<p){BVD.ind<-c(BVD.ind,rep(0,times=(p-length(BVD.ind))))}
  # BVD common eigenvector order: same as findcpc (check which ones are common!)
  
  BCR.ind<-BCR(origdata,reps=1000)[,"common.ind"]
  if(length(BCR.ind)<p){BCR.ind<-c(BCR.ind,rep(0,times=(p-length(BCR.ind))))}
  # BCR common eigenvector order: same as findcpc (check which ones are common!)
  
  bonferroni.sig<-1-0.95^(1/p)
  
  RVC.out<-RVC(covmats,reps=100000)[,"p.values"]
  RVC.ind<-rep(0,times=length(RVC.out))
  RVC.ind[RVC.out<=bonferroni.sig]<-1
  if(length(RVC.ind)<p){RVC.ind<-c(RVC.ind,rep(0,times=(p-length(RVC.ind))))}
  # RVC common eigenvector order: same as findcpc
  
  BootTest.out<-BootTest(origdata)[,"p.values"]
  BootTest.ind<-rep(0,times=length(BootTest.out))
  BootTest.ind[BootTest.out>bonferroni.sig]<-1
  if(length(BootTest.ind)<p){BootTest.ind<-c(BootTest.ind,rep(0,times=(p-length(BootTest.ind))))}
  # BootTest common eigenvector order: same as findcpc
  
  # Return majority vote on number of common eigenvectors (ties broken by choosing maximum mode)
  
  #resultmat<-rbind("Common eigenvector"=commonvec.order,"Flury AIC"=fluryAIC.ind,"BVD 2c"=BVD2c.ind,BCR=BCR.ind,RVC=RVC.ind,BootTest=BootTest.ind)
  resultmat<-rbind("Common eigenvector"=commonvec.order,"Flury AIC"=fluryAIC.ind,"BVD"=BVD.ind,BCR=BCR.ind,RVC=RVC.ind,BootTest=BootTest.ind)
  common.votes<-apply(resultmat[2:6,],2,sum)
  commonvecs<-rep(0,times=p)
  commonvecs[common.votes>2]<-1
  resultmat<-rbind(resultmat,"Common vectors"=commonvecs)
  commonvecnums<-commonvec.order[commonvecs>0]
  return(list(Results=resultmat,commonvecs=commonvecnums,commonvecmat=B[,commonvecnums,drop=F])) # Row 1: order of common eigenvectors in B; Row 2-5: results from AIC, BVD, BCR and RVC tests (1 = eigenvector common); Row 6: ensemble test common eigenvector indicator (1 = eigenvector common)
}
