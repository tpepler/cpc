\name{betasim}
\alias{betasim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Simulate common eigenvector matrices
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Simulates common eigenvector matrices of a specified size and CPC or CPC(q) structure.
}
\usage{
betasim(k, p, q)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{k }{The number of common eigenvector matrices required.}
  \item{p }{Number of rows/columns (variables), i.e. the size of the eigenvector matrices.}
  \item{q }{Number of eigenvectors that should be common to all of the matrices.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
The purpose of this function is for the simulation of common eigenvector matrices for use in simulation studies. CPC or CPC(q) matrices can be generated for several groups.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
Returns an array containing the common eigenvector matrices.
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
}
\examples{
# Generate three (partially) common eigenvector matrices of
# dimension p = 5, with a CPC(2) structure
betasim(k = 3, p = 5, q = 2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
