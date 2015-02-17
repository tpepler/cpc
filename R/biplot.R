biplot<-function(datalist,B,D3=FALSE,varex=1,plotvar=TRUE,main="CPC biplot",col=c("blue","red","green","orange","brown","purple"),radius=0.1,lwd=3)
{
  # Draws a 2- or 3-dimensional biplot of the data in datamat (with different colours indicating the different groups), rotated with the orthogonal matrix B.
  
  # datalist: list of the data from the k groups
  # B: orthogonal projection matrix
  # D3: should a 3-dimensional biplot be drawn?
  # varex: expansion factor for drawing the variables on the biplot
  # plotvar: should the variables be drawn on the biplot?
  # main: title of the biplot
  # col: colors for the data points of the k groups
  
  k<-length(datalist)
  datamat<-NULL
  nvec<-rep(NA,times=k)
  for(i in 1:k){
    datamat<-rbind(datamat,datalist[[i]])
    nvec[i]<-nrow(datalist[[i]])
  }
  p<-ncol(datamat)
  varnames<-colnames(datamat)
  
  # Standardize the columns of datamat by subtracting the column means
  datamat<-t(t(datamat)-apply(datamat,2,mean))
  
  plotpoints<-as.matrix(datamat)%*%B
  
  #t(t(plotpoints)-apply(plotpoints,2,mean))
  #return(plotpoints)
  
  # 3-dimensional biplot
  if(D3){
    library(rgl)
    rgl.open()
    rgl.bg(color="white",sphere=T)
    plot3d(x=plotpoints[,1],y=plotpoints[,2],z=plotpoints[,3],xlab="PC 1",ylab="PC 2",zlab="PC 3",type="n")
    decorate3d(main=main,xlab=NULL,ylab=NULL,zlab=NULL)
    begin<-1
    for(i in 1:k){
      end<-begin+nvec[i]-1
      #points3d(x=plotpoints[begin:end,1],y=plotpoints[begin:end,2],z=plotpoints[begin:end,3],col=col[i])
      spheres3d(x=plotpoints[begin:end,1],y=plotpoints[begin:end,2],z=plotpoints[begin:end,3],col=col[i],radius=radius)
      begin<-end+1    
    }
    
    if(plotvar){
      for(j in 1:p){
        lines3d(x=c(0,B[j,1]*varex),y=c(0,B[j,2]*varex),z=c(0,B[j,3]*varex),lwd=lwd)
      }
      text3d(x=B[,1]*varex,y=B[,2]*varex,z=B[,3]*varex,texts=varnames)
    }
    #cat(paste("Total amount of variation explained in the 3D biplot:",round(biplotvar.total(X=datamat,B=B,r=3)$explained*100,2),"%\n",sep=""))	# Depends on biplotvar.total function!
    #cat(paste("Amount of variation explained within the groups in the 3D biplot:",round(biplotvar.within(X=datamat,nvek=nvek,B=B,r=3)$explained*100,2),"%\n",sep=""))	# Depends on biplotvar.within function!
    #vars.quality<-biplotvar.vars(B,varnames=varnames,r=3)	# Depends on biplotvar.vars function!
    #cat("\nQuality of representation of the variables in the 3D biplot:\n\n")
    #print(round(vars.quality,3)*100)
  }
  
  # 2-dimensional biplot
  else{
    library(MASS)
    par(pch=20)
    eqscplot(x=plotpoints[,1],y=plotpoints[,2],type="n",xlab="PC 1",ylab="PC 2",main=main)
    
    begin<-1
    for(i in 1:k){
      end<-begin+nvec[i]-1
      points(x=plotpoints[begin:end,1],y=plotpoints[begin:end,2],type="p",col=col[i])
      begin<-end+1    
    }
    
    if(plotvar){
      arrows(x0=rep(0,times=p),x1=B[,1]*varex,y0=rep(0,times=p),y1=B[,2]*varex,lwd=lwd)
      text(x=B[,1]*varex,y=B[,2]*varex,pos=3,labels=varnames)
    }
    #cat(paste("Total amount of variation explained in the 2D biplot:",round(biplotvar.total(X=datamat,B=B,r=2)$explained*100,2),"%\n",sep=""))	# Depends on biplotvar.total function!
    #cat(paste("Amount of variation explained within the groups in the 2D biplot:",round(biplotvar.within(X=datamat,nvek=nvek,B=B,r=2)$explained*100,2),"%\n",sep=""))	# Depends on biplotvar.within function!
    #vars.quality<-biplotvar.vars(B,varnames=varnames,r=2)	# Depends on biplotvar.vars function!
    #cat("\nQuality of representation of the variables in the 2D biplot:\n\n")
    #print(round(vars.quality,3)*100)
  }
}
