\name{ensemble.test}
\alias{ensemble.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ensemble.test(origdata, standardize = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{origdata}{
%%     ~~Describe \code{origdata} here~~
}
  \item{standardize}{
%%     ~~Describe \code{standardize} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (origdata, standardize = F) 
{
    myFunc <- function(datavec) {
        return((datavec - mean(datavec))/sd(datavec))
    }
    standcol <- function(datamat) {
        return(apply(datamat, 2, myFunc))
    }
    nvec <- c(nrow(origdata[[1]]), nrow(origdata[[2]]))
    p <- ncol(origdata[[1]])
    covmats <- array(NA, dim = c(p, p, 2))
    if (standardize) {
        origdata[[1]] <- standcol(origdata[[1]])
        origdata[[2]] <- standcol(origdata[[2]])
    }
    covmats[, , 1] <- cov(origdata[[1]])
    covmats[, , 2] <- cov(origdata[[2]])
    B <- FG(covmats, nvec)$B
    commonvec.order <- findcpc(covmats, B = B, plotting = F)$commonvec.order
    flury.out <- flury.test(covmats, nvec, B = B, p = dim(covmats)[1], 
        qmax = p - 2, commonvec.order = commonvec.order)
    fluryAIC.vote <- flury.out[which.min(flury.out[, "AIC"]), 
        "No.of.CPCs"]
    fluryAIC.ind <- c(rep(1, times = fluryAIC.vote), rep(0, times = (p - 
        fluryAIC.vote)))
    BVD.out <- BVD(origdata, reps = 1000)
    BVD2c.vote <- sum(BVD.out["BVD 2c", ])
    BVD2c.ind <- BVD.out["BVD 2c", ]
    if (length(BVD2c.ind) < p) {
        BVD2c.ind <- c(BVD2c.ind, rep(0, times = (p - length(BVD2c.ind))))
    }
    BCR.ind <- BCR(origdata, reps = 1000)[, "common.ind"]
    BCR.vote <- sum(BCR.ind)
    if (length(BCR.ind) < p) {
        BCR.ind <- c(BCR.ind, rep(0, times = (p - length(BCR.ind))))
    }
    RVC.out <- RVC(covmats, reps = 1e+05)[, "p.values"]
    RVC.ind <- rep(0, times = p)
    RVC.ind[RVC.out <= 0.05] <- 1
    BootTest.out <- BootTest(origdata)[, "p.values"]
    BootTest.ind <- rep(0, times = p)
    BootTest.ind[BootTest.out > 0.05] <- 1
    resultmat <- rbind(`Common eigenvector` = commonvec.order, 
        `Flury AIC` = fluryAIC.ind, `BVD 2c` = BVD2c.ind, BCR = BCR.ind, 
        RVC = RVC.ind, BootTest = BootTest.ind)
    common.votes <- apply(resultmat[2:6, ], 2, sum)
    commonvecs <- rep(0, times = p)
    commonvecs[common.votes > 2] <- 1
    resultmat <- rbind(resultmat, `Common vectors` = commonvecs)
    commonvecnums <- commonvec.order[commonvecs > 0]
    return(list(Results = resultmat, commonvecs = commonvecnums, 
        commonvecmat = B[, commonvecnums]))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
