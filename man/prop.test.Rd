\name{prop.test}
\alias{prop.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Likelihood ratio test of proportionality of several covariance matrices
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Calculates the likelihood ratio statistic and its degrees of freedom for the hypothesis of proportional covariance matrices against the alternative of unrelated covariance matrices.
}
\usage{
prop.test(covmats, nvec)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{covmats}{
%%     ~~Describe \code{covmats} here~~
Array of covariance matrices.
}
  \item{nvec}{
%%     ~~Describe \code{nvec} here~~
Vector of sample sizes of the k groups.
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
Returns a list with the following:
\item{chi.square}{The likelihood ratio test statistic.}
\item{df}{Degrees of freedom of the test statistic under the null hypothesis.}
\item{covmats.prop}{Estimated covariance matrices under the null hypothesis model.}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
Theo Pepler
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{flury.test}}, \code{\link{equal.test}}, \code{\link{cpc.test}} and \code{\link{cpcq.test}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (covmats, nvec) 
{
    k <- dim(covmats)[3]
    p <- dim(covmats)[1]
    rvec <- (nvec - 1)/(sum(nvec) - k)
    rho <- rep(1, times = k)
    maxrho <- 1
    prevmaxrho <- 100
    while (abs(maxrho - prevmaxrho) > 1e-09) {
        covmats.total <- 0
        for (i in 1:k) {
            covmats.total <- covmats.total + rvec[i] * covmats[, 
                , i]/rho[i]
        }
        B <- eigen(covmats.total)$vectors
        avec <- matrix(NA, ncol = p, nrow = k)
        for (i in 1:k) {
            avec[i, ] <- diag(t(B) \%*\% covmats[, , i] \%*\% B)
        }
        lambda <- rep(NA, times = p)
        for (j in 1:p) {
            lambda[j] <- sum(rvec * avec[, j]/rho)
        }
        for (i in 2:k) {
            rho[i] <- 1/p * sum(avec[i, ]/lambda)
        }
        prevmaxrho <- maxrho
        maxrho <- max(rho[2:k])
    }
    covmats.prop <- array(NA, dim = c(p, p, k))
    for (i in 1:k) {
        covmats.prop[, , i] <- B \%*\% (diag(lambda) * rho[i]) \%*\% 
            t(B)
    }
    chi2total <- 0
    for (i in 1:k) {
        chi2total <- chi2total + (nvec[i] - 1) * log(det(covmats.prop[, 
            , i])/det(covmats[, , i]))
    }
    df <- k * (0.5 * p * (p - 1) + p) - (0.5 * p * (p - 1) + 
        p + k - 1)
    return(list(chi.square = chi2total, df = df, covmats.prop = covmats.prop))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
