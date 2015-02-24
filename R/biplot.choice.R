biplot.choice <- function(datalist, rdim, add.projectmats = NULL)
{
  # Gives biplot goodness of fit measures for different types of principal components biplots for grouped data
  
  # datalist: list containing all data
  # rdim: number of dimensions of biplot
  # add.projectmats: additional orthogonal projections matrices to compute biplot fit measures for
  
  k <- length(datalist)
  p <- ncol(datalist[[1]])
  nvec <- rep(NA, times = k)
  
  # Eigenvectors of pooled covariance matrix
  
  dfpooled <- 0
  SSpooled <- 0
  for(i in 1:k){
    nvec[i] <- nrow(datalist[[i]])
    dfpooled <- dfpooled + nvec[i] - 1
    SSpooled <- SSpooled + cov(as.matrix(datalist[[i]])) * (nvec[i] - 1)
  }
  Sp <- SSpooled / dfpooled
  Ep <- eigen(Sp)$vectors
  pooledcov.output <- biplot.measures(datalist = datalist, projectmat = Ep, rdim = rdim)
  
  # Eigenvectors of covariance matrix of pooled data
  
  datamat <- NULL
  for(i in 1:k){
    datamat <- rbind(datamat, as.matrix(datalist[[i]]))
  }
  E <- eigen(cov(datamat))$vectors
  pooleddata.output <- biplot.measures(datalist = datalist, projectmat = E, rdim = rdim)
  
  # CPC: FG algorithm
  
  S <- array(NA, dim = c(p, p, k))
  for(i in 1:k){
    S[, , i] <- cov(as.matrix(datalist[[i]]))
  }
  B.flury <- cpc::FG(covmats = S, nvec = nvec)$B
  flury.output <- biplot.measures(datalist = datalist, projectmat = B.flury, rdim = rdim)
  
  # CPC: Stepwise CPC
  
  B.stepwise <- stepwisecpc(covmats = S, nvec = nvec)$B#[[1]]
  stepwise.output <- biplot.measures(datalist = datalist, projectmat = B.stepwise, rdim = rdim)
  
  # CPC: JADE algorithm
  
  library(JADE)
  B.jade <- rjd(X = S)$V
  lvals <- rep(0, times = p)
  for(i in 1:k){
    lvals <- lvals + diag(t(B.jade) %*% S[, , i] %*% B.jade)
  }
  jade.order <- order(lvals, decreasing = TRUE)
  B.jade <- B.jade[, jade.order]
  jade.output <- biplot.measures(datalist = datalist, projectmat = B.jade, rdim = rdim)
  
  # Produce output table
  
  pooledcov <- c(pooledcov.output$overall.quality,
               pooledcov.output$within.quality.mean,
               pooledcov.output$between.quality,
               pooledcov.output$adequacies.median,
               pooledcov.output$mspe.mean,
               #pooledcov.output$axis.predictivities.mean,
               pooledcov.output$sample.predictivities.mean)
  
  pooleddata <- c(pooleddata.output$overall.quality,
                pooleddata.output$within.quality.mean,
                pooleddata.output$between.quality,
                pooleddata.output$adequacies.median,
                pooleddata.output$mspe.mean,
                #pooleddata.output$axis.predictivities.mean,
                pooleddata.output$sample.predictivities.mean)
  
  flury <- c(flury.output$overall.quality,
           flury.output$within.quality.mean,
           flury.output$between.quality,
           flury.output$adequacies.median,
           flury.output$mspe.mean,
           #flury.output$axis.predictivities.mean,
           flury.output$sample.predictivities.mean)
  
  stepwise <- c(stepwise.output$overall.quality,
              stepwise.output$within.quality.mean,
              stepwise.output$between.quality,
              stepwise.output$adequacies.median,
              stepwise.output$mspe.mean,
              #stepwise.output$axis.predictivities.mean,
              stepwise.output$sample.predictivities.mean)
  
  jade <- c(jade.output$overall.quality,
          jade.output$within.quality.mean,
          jade.output$between.quality,
          jade.output$adequacies.median,
          jade.output$mspe.mean,
          #jade.output$axis.predictivities.mean,
          jade.output$sample.predictivities.mean)
  
  resultsmat <- rbind(pooledcov, pooleddata, flury, stepwise, jade)
  rownames(resultsmat) <- c("Pooled S", "Pooled data", "Flury", "Stepwise CPC", "JADE")
  #colnames(resultsmat)<-c("Overall","Within","Between","Adequacy","Axis predictivities","Sample predictivities")
  colnames(resultsmat) <- c("Overall", "Within", "Between", "Adequacy", "MSPE", "Sample predictivities")
  
  return(resultsmat)
}
