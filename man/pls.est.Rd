\name{pls.est}
\alias{pls.est}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Partial least squares (PLS) algorithm
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Implementation of the partial least squares (PLS) algorithm as described in Hastie et al. (2009).
}
\usage{
pls.est(X, Y)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X }{Matrix containing the independent variables.}
  \item{Y }{Vector containing the dependent/response values.}
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
Returns a list with the following items:
\item{pls.scores }{The partial least squares scores.}
\item{pls.loadings }{The partial least squares loadings.}
}
\references{
%% ~put references to the literature/web site here ~
Hastie, T., Tibshirani, R., and Friedman, J. (2009). The Elements of Statistical Learning. Springer.
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
}
\examples{
# The Iris data
data(iris)
indep <- iris[, 1:3]
dep <- iris[, 4]
pls.est(X = indep, Y = dep)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
