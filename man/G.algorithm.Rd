\name{G.algorithm}
\alias{G.algorithm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{G algorithm}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Implementation of the G algorithm as internal function of the FG algorithm.
}
\usage{
%G.algorithm(T.mat, nvec)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
%  \item{T.mat}{
%%     ~~Describe \code{T.mat} here~~
%}
%  \item{nvec}{
%%     ~~Describe \code{nvec} here~~
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
function (T.mat, nvec) 
{
    k <- dim(T.mat)[3]
    Q <- matrix(rep(200, 4), nrow = 2, ncol = 2)
    Q.nuut <- diag(2)
    delta <- matrix(NA, nrow = k, ncol = 2)
    while (abs(Q[2, 1] - Q.nuut[2, 1]) > 1e-07) {
        Q <- Q.nuut
        for (i in 1:k) {
            for (j in 1:2) {
                delta[i, j] <- t(Q[, j]) \%*\% T.mat[, , i] \%*\% 
                  Q[, j]
            }
        }
        U <- matrix(rep(0, 4), nrow = 2, ncol = 2)
        for (i in 1:k) {
            U <- U + nvec[i] * ((delta[i, 1] - delta[i, 2])/(delta[i, 
                1] * delta[i, 2])) * T.mat[, , i]
        }
        Q.nuut <- eigvec(U)
    }
    return(Q.nuut)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
