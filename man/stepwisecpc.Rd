\name{stepwisecpc}
\alias{stepwisecpc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Stepwise common principal components (CPC)
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Estimates stepwise common principal components to simultaneously diagonalise several covariance matrices.
}
\usage{
stepwisecpc(covmats, nvec)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{covmats}{
%%     ~~Describe \code{covmats} here~~
array containing the covariance matrices for the groups, created with a command such as covmats<-array(NA,dim=c(p,p,k)), where p refers to the number of rows/columns of each covariance matrix, and k is the number of groups (or covariance matrices)
}
  \item{nvec}{
%%     ~~Describe \code{nvec} here~~
vector containing the sample sizes of the k groups
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
Returns a list with the values:
\item{B}{orthogonal matrix containing the estimated common eigenvectors in the columns}
\item{eigenvals}{matrix of which the columns contain the estimated eigenvalues of the k covariance matrices under the CPC model}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Sarel Steel, Theo Pepler
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
function (covmats, nvec) 
{
    p <- dim(covmats)[1]
    k <- length(nvec)
    ntot <- sum(nvec)
    Spooled <- matrix(0, nrow = p, ncol = p)
    for (j in 1:k) {
        Spooled <- Spooled + (nvec[j] - 1) * S[, , j]
    }
    Spooled <- Spooled/(ntot - k)
    Qmat <- matrix(0, nrow = p, ncol = p)
    Qtilde <- eigen(Spooled, symmetric = TRUE)$vectors
    pimat <- diag(p)
    muvec <- rep(0, k)
    elmax <- 10
    for (j in 1:p) {
        xvec <- as.vector(Qtilde[, j])
        xvec <- pimat \%*\% xvec
        for (i in 1:k) {
            muvec[i] <- t(xvec) \%*\% covmats[, , i] \%*\% xvec
        }
        for (el in 1:elmax) {
            Smat <- matrix(0, nrow = p, ncol = p)
            for (i in 1:k) {
                Smat <- Smat + (nvec[i] - 1) * covmats[, , i]/muvec[i]
            }
            yvec <- pimat \%*\% Smat \%*\% xvec
            xvec <- yvec/as.numeric(sqrt(t(yvec) \%*\% yvec))
            for (i in 1:k) {
                muvec[i] <- t(xvec) \%*\% covmats[, , i] \%*\% xvec
            }
        }
        Qmat[, j] <- xvec
        pimat <- pimat - xvec \%*\% t(xvec)
    }
    eigenvals <- matrix(0, p, k)
    for (i in 1:k) {
        eigenvals[, i] <- t((diag(t(Qmat) \%*\% covmats[, , i] \%*\% 
            Qmat)))
    }
    results <- list(B = Qmat, eigenvals = eigenvals)
    return(results)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
