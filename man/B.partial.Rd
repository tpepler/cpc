\name{B.partial}
\alias{B.partial}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Modal matrices under partial CPC model
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Estimates the modal matrices under the partial CPC model, to diagonalise several covariance matrices.
}
\usage{
B.partial(covmats, nvec, B = FG(covmats = covmats, nvec = nvec)$B, commonvec.order, q)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{covmats}{
%%     ~~Describe \code{covmats} here~~
Array of sample covariance matrices for the k groups.
}
  \item{nvec}{
%%     ~~Describe \code{nvec} here~~
Vector of sample sizes of the k groups.
}
  \item{B}{
%%     ~~Describe \code{B} here~~
Matrix of common eigenvectors estimated under the assumption of full CPC. Defaults to the modal matrix obtained with the FG algorithm.
}
  \item{commonvec.order}{
%%     ~~Describe \code{commonvec.order} here~~
Vector containing the order of the common eigenvectors in B (with the q truly common eigenvectors in the first q positions).
}
  \item{q}{
%%     ~~Describe \code{q} here~~
Number of eigenvectors common to all k groups.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
Estimates the matrices of common (and non-common) eigenvectors for each of the groups, according to the method described in Flury (1988).
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
Returns an array containing the modal matrices for the k groups, estimated under the CPC(q) model.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Theo Pepler
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
function (covmats, nvec, B = FG(covmats = covmats, nvec = nvec)$B, 
    commonvec.order, q) 
{
    k <- dim(covmats)[3]
    p <- dim(covmats)[1]
    B <- B[, commonvec.order]
    Bmats <- array(NA, dim = c(p, p, k))
    for (i in 1:k) {
        B1 <- B[, 1:q, drop = F]
        B2 <- B[, (q + 1):p]
        Q1 <- eigen(t(B2) \%*\% covmats[, , i] \%*\% B2)$vectors
        B21 <- B2 \%*\% Q1
        Bmats[, , i] <- cbind(B1, B21)
    }
    return(Bmats)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
