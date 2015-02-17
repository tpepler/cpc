\name{box.mtest}
\alias{box.mtest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Box's M test
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Chi-square approximation to Box's M test for equality of k covariance matrices, as described in Rencher (2002, p255).
}
\usage{
box.mtest(covmats, nvec)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{covmats}{
%%     ~~Describe \code{covmats} here~~
Array of k covariance matrices to be tested for equality.
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
Returns the p value for the test.
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
function (covmats, nvec) 
{
    k <- dim(covmats)[3]
    p <- dim(covmats)[2]
    ntot <- sum(nvec)
    Sp <- matrix(0, nrow = p, ncol = p)
    for (i in 1:k) {
        Sp <- Sp + (nvec[i] - 1) * covmats[, , i]
    }
    Sp <- Sp/(ntot - k)
    temp <- 0
    for (i in 1:k) {
        temp <- temp + 1/(nvec[i] - 1)
    }
    c1 <- (temp - (1/(ntot - k))) * (2 * (p^2) + 3 * p - 1)/(6 * 
        (p + 1) * (k - 1))
    temp2 <- 0
    for (i in 1:k) {
        temp2 <- temp2 + (nvec[i] - 1) * log(det(covmats[, , 
            i]))
    }
    lnM <- 0.5 * (temp2) - 0.5 * (ntot - k) * log(det(Sp))
    chi2 <- -2 * (1 - c1) * lnM
    v <- 0.5 * (k - 1) * p * (p + 1)
    pval <- pchisq(q = chi2, df = v, lower.tail = F)
    return(pval)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
