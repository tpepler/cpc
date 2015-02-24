biplot.measures <- function(datalist, projectmat, rdim)
{
  # Calculates goodness of fit measures for r-dimensional principal component biplot of grouped data
  
  # datalist: list containing all data
  # projectmat: orthogonal projection matrix for biplot
  # rdim: number of dimensions of biplot
  
  k <- length(datalist)
  p <- ncol(datalist[[1]])
  if(rdim > p){
    cat("Number of biplot dimensions cannot be larger than number of variables in data!\n")
    return()
  }
  varnames <- colnames(datalist[[1]])
  
  # Overall quality of biplot display
  
  X <- NULL
  for(i in 1:k){
    X <- rbind(X, as.matrix(datalist[[i]]))
  }
  n <- nrow(X)
  Xmean <- apply(X, 2, mean)
  X <- t(t(X) - Xmean)
  Y <- X %*% projectmat[, 1:rdim]
  totalvar <- sum(diag(t(X) %*% X))
  fittedvar <- sum(diag(t(Y) %*% Y))
  overall.quality <- fittedvar / totalvar # Overall quality of display of points in biplot (within and between group variation)
  
  # Mean quality of biplot display of variation within groups
  
  withingroup.total <- rep(NA, times = k)
  withingroup.fitted <- rep(NA, times = k)
  within.total <- 0
  within.fitted <- 0
  between.total <- 0
  between.fitted <- 0
  Xmean.fit <- t(projectmat[, 1:rdim]) %*% Xmean
  for(i in 1:k){
    groupdata <- as.matrix(datalist[[i]])
    groupmean <- apply(groupdata, 2, mean)
    groupdata <- t(t(groupdata) - groupmean)
    fitteddata <- groupdata %*% projectmat[, 1:rdim]
    
    withingroup.total[i] <- sum(diag(t(groupdata) %*% groupdata))
    within.total <- within.total + withingroup.total[i]
    withingroup.fitted[i] <- sum(diag(t(fitteddata) %*% fitteddata))
    within.fitted <- within.fitted + withingroup.fitted[i]
    
    between.total <- between.total + t(groupmean - Xmean) %*% (groupmean - Xmean)
    groupfitted.mean <- t(projectmat[, 1:rdim]) %*% groupmean
    between.fitted <- between.fitted + t(groupfitted.mean - Xmean.fit) %*% (groupfitted.mean - Xmean.fit)
  }
  within.quality <- withingroup.fitted / withingroup.total # Within group variation quality of display
  within.quality <- matrix(within.quality, ncol = k, byrow = TRUE)
  colnames(within.quality) <- paste("Group ", c(1:k), sep = "")
  within.quality.mean <- within.fitted / within.total  # Mean quality of display of points (within group variation)
  
  # Quality of between group variation displayed in biplot
  
  between.quality <- as.numeric(between.fitted / between.total) # Overall quality of between group variation as represented in r-dimensional biplot
  
  # Adequacies of the variables
  
  adequacies <- diag(projectmat[, 1:rdim] %*% t(projectmat[, 1:rdim])) # adequacy of the variables as represented in a r-dimensional biplot
  adequacies.median <- median(adequacies)
  adequacies <- matrix(adequacies, ncol = p, byrow = TRUE)
  colnames(adequacies) <- varnames
  
  # Axis predictivities (predivities of the variables)
  
  J <- diag(c(rep(1, times = rdim), rep(0, times = (p - rdim))))
  X.fitted <- X %*% projectmat %*% J %*% t(projectmat)
  axis.predictivities <- diag(diag(diag(t(X.fitted) %*% X.fitted)) %*% solve(diag(diag(t(X) %*% X)))) # Axis predictivities
  axis.predictivities.mean <- sum(diag(diag(diag(t(X.fitted) %*% X.fitted)) %*% solve(diag(diag(t(X) %*% X))))) / p # Mean predictivity of axes (variables)
  
  axis.predictivities <- matrix(axis.predictivities, ncol = p, byrow = TRUE)
  colnames(axis.predictivities) <- varnames
  
  # Sample predictivities (predictivities of the observations)
  
  sample.predictivities <- diag(diag(diag(X.fitted %*% t(X.fitted))) %*% solve(diag(diag(X %*% t(X))))) # Sample predictivities
  sample.predictivities.mean <- sum(diag(diag(diag(X.fitted %*% t(X.fitted))) %*% solve(diag(diag(X %*% t(X)))))) / nrow(X) # Mean predictivity of samples (observations)
  
  # Mean standard predictive errors (Alves 2012)
  
  Xsd <- apply(X, 2, sd)
  onevec <- matrix(1, nrow = n, ncol = 1)
  mspe <- rep(NA, times = p)
  for(j in 1:p){
    mspe[j] <- (t(onevec) %*% abs(X[, j] - X %*% projectmat[, 1:rdim] %*% t(projectmat[j, 1:rdim, drop = FALSE]))) / (n * Xsd[j])
  }
  mspe.mean <- mean(mspe)
  
  return(list(overall.quality = overall.quality,
              within.quality = within.quality,
              within.quality.mean = within.quality.mean,
              between.quality = between.quality,
              adequacies = adequacies,
              adequacies.median = adequacies.median,
              axis.predictivities = axis.predictivities,
              axis.predictivities.mean = axis.predictivities.mean,
              sample.predictivities = sample.predictivities,
              sample.predictivities.mean = sample.predictivities.mean,
              mspe = mspe,
              mspe.mean = mspe.mean))
}
