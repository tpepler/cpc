\name{frobenius}
\alias{frobenius}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
(Modified) Frobenius matrix norm
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Calculates a modified version of the Frobenius matrix norm for a square symmetric matrix with p(p + 1) / 2 estimable parameters.
}
\usage{
frobenius(datamat, targetmat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{datamat }{Symmetric square matrix for which to calculate the Frobenius norm.}
  \item{targetmat }{Target matrix of same size as \code{datamat}, to compare \code{datamat} against.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
This is a utility function used by \code{\link{alpha.crossvalid}}, but can also be called directly if required. This version of the Frobenius matrix norm gives equal weight to all parameters (variances and covariances) when used with covariance matrices.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
Returns the modified Frobenius matrix norm (scalar).
}
\references{
%% ~put references to the literature/web site here ~
Pepler, P.T. (2014). The identification and application of common principal components. PhD dissertation in the Department of Statistics and Actuarial Science, Stellenbosch University.
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
\code{\link{alpha.crossvalid}}
}
\examples{
# Compare the covariance matrix of the versicolor group to
# an identity matrix
data(iris)
versicolor <- iris[51:100, 1:4]
frobenius(datamat = cov(versicolor), targetmat = diag(4))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
