\name{offdiag.vec}
\alias{offdiag.vec}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Off-diagonal element stacking
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Stacks the rows of a square matrix (excluding diagonal elements) in a vector.
}
\usage{
offdiag.vec(datamat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{datamat }{Any square matrix.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
This is a utility function used by \code{\link{alpha.schafer}}, but can also be called directly if required.
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
\code{\link{alpha.schafer}}
}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
