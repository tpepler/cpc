\name{FG}
\alias{FG}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Flury-Gautschi (FG) algorithm}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
The Flury-Gautschi (FG) algorithm for the simultaneous diagonalisation of several covariance matrices.
}
\usage{
FG(covmats, nvec)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{covmats}{
%%     ~~Describe \code{covmats} here~~
array of covariance matrices to be simultaneously diagonalized, created by a command such as covmats<-array(NA,dim=c(p,p,k))
}
  \item{nvec}{
%%     ~~Describe \code{nvec} here~~
vector of sample sizes for the covariance matrices in covmats
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
A modal matrix (B) is calculated for the covariance matrices under the assumption of common principal components in the populations. The common eigenvectors in B are ordered according to the sums of the associated eigenvalues in the k groups, from the largest sum to the smallest.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
Returns a list with the values:
\item{B}{orthogonal matrix containing the estimated common eigenvectors in the columns}
\item{diagvals}{matrix of which the columns contain the estimated eigenvalues of the k covariance matrices under the CPC model}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Danie Uys, Theo Pepler
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
    p <- dim(covmats)[2]
    B <- diag(p)
    DIFF <- 100
    k <- dim(covmats)[3]
    while (DIFF > 1e-09) {
        B.old <- B
        T.mat <- array(NA, dim = c(2, 2, k))
        for (m in 1:(p - 1)) {
            for (j in (m + 1):p) {
                vek <- c(m, j)
                for (i in 1:k) {
                  T.mat[, , i] <- t(B[, vek]) \%*\% covmats[, , 
                    i] \%*\% B[, vek]
                }
                J <- G.algorithm(T.mat, nvec)
                B[, vek] <- B[, vek] \%*\% J
            }
        }
        for (i in 1:p) {
            for (j in 1:p) {
                DIFF <- abs(B[i, j] - B.old[i, j])
            }
        }
    }
    diagvals <- 0
    for (i in 1:k) {
        diagvals <- diagvals + diag(t(B) \%*\% covmats[, , i] \%*\% 
            B)
    }
    B <- B[, order(diagvals, decreasing = TRUE)]
    diagvalsmat <- matrix(NA, nrow = p, ncol = k)
    for (i in 1:k) {
        diagvalsmat[, i] <- diag(t(B) \%*\% covmats[, , i] \%*\% 
            B)
    }
    return(list(B = B, diagvals = diagvalsmat))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
