\name{bootveccor}
\alias{bootveccor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
bootveccor(origdata, veccormat, nvec, reps = 1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{origdata}{
%%     ~~Describe \code{origdata} here~~
}
  \item{veccormat}{
%%     ~~Describe \code{veccormat} here~~
}
  \item{nvec}{
%%     ~~Describe \code{nvec} here~~
}
  \item{reps}{
%%     ~~Describe \code{reps} here~~
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
function (origdata, veccormat, nvec, reps = 1000) 
{
    numcomb <- nrow(veccormat)
    p <- ncol(origdata[[1]])
    k <- length(nvec)
    bootreps <- matrix(NA, ncol = numcomb, nrow = reps)
    for (r in 1:reps) {
        group.PCA <- array(NA, dim = c(p, p, k))
        for (i in 1:k) {
            bootdata <- origdata[[i]][sample(c(1:nvec[i]), size = nvec[i], 
                replace = T), ]
            group.PCA[, , i] <- eigen(cov(bootdata))$vectors
        }
        for (c in 1:numcomb) {
            bootreps[r, c] <- abs(t(group.PCA[, veccormat[c, 
                1], 1]) \%*\% group.PCA[, veccormat[c, 2], 2])
        }
    }
    return(bootreps)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
