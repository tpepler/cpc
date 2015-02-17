\name{flury.AIC}
\alias{flury.AIC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Flury's AIC statistics
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Calculates the Akaike Information Criterion (AIC) for the higher model vs. the lower model in Flury's hierarchy of covariance matrices.
}
\usage{
flury.AIC(covmats.high, covmats.low, nvec, df)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{covmats.high}{
%%     ~~Describe \code{covmats.high} here~~
Array of estimated covariance matrices under the higher model in Flury's hierarchy.
}
  \item{covmats.low}{
%%     ~~Describe \code{covmats.low} here~~
Array of estimated covariance matrices under the lower model in Flury's hierarchy.
}
  \item{nvec}{
%%     ~~Describe \code{nvec} here~~
Vector of sample sizes of the k groups.
}
  \item{df}{
%%     ~~Describe \code{df} here~~
Degrees of freedom of the higher model, versus the model of unrelated covariance matrices.
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
Returns the AIC value (scalar).
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
\code{\link{flury.test}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (covmats.high, covmats.low, nvec, k, p, df) 
{
    aic.total <- 0
    for (i in 1:k) {
        aic.total <- aic.total + (nvec[i] - 1) * (sum(diag(solve(covmats.high[, 
            , i]) \%*\% covmats.low[, , i])) + log(det(covmats.high[, 
            , i])) - sum(diag(solve(covmats.low[, , i]) \%*\% covmats.low[, 
            , i])) - log(det(covmats.low[, , i])))
    }
    npar <- k * (0.5 * p * (p - 1) + p) - (0.5 * p * (p - 1) + 
        p) - df
    aic.criterion <- aic.total + 2 * npar
    return(aic.criterion)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
