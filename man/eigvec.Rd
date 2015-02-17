\name{eigvec}
\alias{eigvec}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Internal function of the FG algorithm.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
%eigvec(U)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
%  \item{U}{
%%     ~~Describe \code{U} here~~
%}
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
function (U) 
{
    if (U[1, 2] != 0) {
        ratio <- (U[2, 2] - U[1, 1])/U[1, 2]
    }
    else {
        ratio = 0
    }
    discr <- sqrt(ratio^2 + 4)
    T1 <- (ratio + discr)/2
    T2 <- (ratio - discr)/2
    if (abs(T1) > abs(T2)) {
        T.waarde <- T2
    }
    else {
        T.waarde <- T1
    }
    C.waarde <- 1/sqrt(1 + T.waarde^2)
    S <- T.waarde * C.waarde
    J <- matrix(NA, nrow = 2, ncol = 2)
    J[1, 1] <- C.waarde
    J[2, 1] <- S
    J[1, 2] <- -S
    J[2, 2] <- C.waarde
    return(J)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
