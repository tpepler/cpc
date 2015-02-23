discriminant.qda <- function(origdata, group, method = c("unbiased", "pooled", "cpc", "fullcpccrossvalid"), B = NULL, standardize = FALSE)
{
  # origdata: matrix containing the sample data for two groups
  # group: vector (with values 1 and 2) indicating group membership for the rows in origdata
  # method: unbiased, cpc, fullcpccrossvalid, pooled
  # standardize: should the standardised covariance matrices (i.e. the correlation matrices) be used?
  
  n.all1 <- nrow(origdata[group == 1, ])
  n.all2 <- nrow(origdata[group == 2, ])
  n.train1 <- round(n.all1 * 0.7, 0)
  n.train2 <- round(n.all2 * 0.7, 0)
  n.test1 <- n.all1 - n.train1
  n.test2 <- n.all2 - n.train2
  p <- ncol(origdata)
  group.predict <- rep(NA, times = (n.test1 + n.test2))
  
  myFunc <- function(datavec){
    return(datavec / sd(datavec))
  }
  
  standcol <- function(datamat){
    return(apply(datamat, 2, myFunc))
  }
  
  if(standardize){
    group1data <- origdata[group == 1, ]
    group2data <- origdata[group == 2, ]
    group1data <- standcol(group1data)
    group2data <- standcol(group2data)
    temp.origdata <- rbind(group1data, group2data)
  }
  else{temp.origdata <- origdata}
  
  group1data <- temp.origdata[group == 1, ][1:n.train1, ]
  group2data <- temp.origdata[group == 2, ][1:n.train2, ]
  S1 <- cov(group1data)
  S2 <- cov(group2data)
  testdata <- rbind(temp.origdata[group == 1, ][(n.train1 + 1):n.all1, ], temp.origdata[group == 2, ][(n.train2 + 1):n.all2, ])
  
  if(method[1] == "pooled"){
    Sp <- (S1 * (nrow(group1data) - 1) + S2 * (nrow(group2data) - 1)) / (nrow(group1data) + nrow(group2data) - 2)
    S1 <- Sp
    S2 <- Sp
  }
  
  if(method[1] == "cpc"){
    S1 <- B %*% diag(diag(t(B) %*% S1 %*% B)) %*% t(B)
    S2 <- B %*% diag(diag(t(B) %*% S2 %*% B)) %*% t(B)
  }
  
  if(method[1] == "fullcpccrossvalid"){
    S1.cpc <- B %*% diag(diag(t(B) %*% S1 %*% B)) %*% t(B)
    S2.cpc <- B %*% diag(diag(t(B) %*% S2 %*% B)) %*% t(B)
    alpha1 <- alpha.crossvalid(group1data, B = B, reps = 100)
    alpha2 <- alpha.crossvalid(group2data, B = B, reps = 100)
    S1 <- alpha1 * S1 + (1 - alpha1) * S1.cpc
    S2 <- alpha2 * S2 + (1 - alpha2) * S2.cpc
  }
  
  S1.inv <- solve(S1)
  S2.inv <- solve(S2)
  xbar1 <- apply(group1data, 2, mean)
  xbar2 <- apply(group2data, 2, mean)
  
  c <- 0.5 * (log(det(S1) / det(S2)) + (xbar1 %*% S1.inv %*% xbar1 - xbar2 %*% S2.inv %*% xbar2))
  
  for(i in 1:(n.test1 + n.test2)){
    newobs <- as.matrix(testdata[i, , drop = FALSE])
    classvalue <- -0.5 * newobs %*% (S1.inv - S2.inv) %*% t(newobs) + (xbar1 %*% S1.inv - xbar2 %*% S2.inv) %*% t(newobs)
    if(classvalue >= c){group.predict[i] <- 1}
    else{group.predict[i] <- 2}
  }
  
  group.test <- c(group[group == 1][(n.train1 + 1):n.all1], group[group == 2][(n.train2 + 1):n.all2])
  discrep <- group.test - group.predict
  misclassrate <- length(discrep[discrep != 0])
  misclassrate <- misclassrate / (n.test1 + n.test2)
  return(list(inputdata = cbind(testdata, group = group.test, group.predict = group.predict), misclassrate))
}
