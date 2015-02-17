rotationmat.maxcor<-function(p)
{
  # Finds rotation matrix with maximum correlation between the variables
  
  # p: number of variables
  
  theta<-pi/4 # 45 degree angle in radians
  rotatemat<-NULL
  
  for(i in 1:(p-1)){
    for(j in (i+1):p){
      tempmat<-diag(rep(1,times=p))#matrix(0,nrow=p,ncol=p)
      tempmat[i,i]<- cos(theta)
      tempmat[j,j]<- cos(theta)
      if((i%%2)==0){
        tempmat[i,j]<- sin(theta)
        tempmat[j,i]<- -sin(theta)
      }
      else{
        tempmat[i,j]<- -sin(theta)
        tempmat[j,i]<- sin(theta)
      }      
      if(is.null(rotatemat)){rotatemat<-tempmat}
      else{rotatemat<-rotatemat%*%tempmat}
    }
  }
  
  return(rotatemat)
}
