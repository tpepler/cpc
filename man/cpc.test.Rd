\name{cpc.test}
\alias{cpc.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Likelihood ratio test of common principal components in of several groups
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Calculates the likelihood ratio statistic and its degrees of freedom for the hypothesis of common eigenvectors in the k groups against the alternative of unrelated covariance matrices.
}
\usage{
cpc.test(covmats, nvec, B)
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
  \item{B}{
%%     ~~Describe \code{B} here~~
Modal matrix simultaneously diagonalising the covariance matrices, estimated under the assumption of common eigenvectors in the k groups.
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
\item{covmats.cpc}{Estimated covariance matrices under the null hypothesis model.}
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
\code{\link{flury.test}}, \code{\link{equal.test}}, \code{\link{prop.test}} and \code{\link{cpcq.test}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (covmats, nvec, B) 
{
    k <- dim(covmats)[3]
    p <- dim(covmats)[1]
    covmats.cpc <- array(NA, dim = c(p, p, k))
    chi2total <- 0
    for (i in 1:k) {
        lambda <- diag(t(B) \%*\% covmats[, , i] \%*\% B)
        covmats.cpc[, , i] <- B \%*\% diag(lambda) \%*\% t(B)
        chi2total <- chi2total + (nvec[i] - 1) * log(det(covmats.cpc[, 
            , i])/det(covmats[, , i]))
    }
    df <- k * (0.5 * p * (p - 1) + p) - (0.5 * p * (p - 1) + 
        k * p)
    return(list(chi.square = chi2total, df = df, covmats.cpc = covmats.cpc))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
