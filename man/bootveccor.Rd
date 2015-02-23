\name{bootveccor}
\alias{bootveccor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Calculate bootstrap vector correlations
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Utility function to calculate the bootstrap vector correlations for the BVD method to identify common eigenvectors.
}
\usage{
bootveccor(origdata, veccormat, nvec, reps = 1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{origdata }{List of the grouped sample data.}
  \item{veccormat }{Matrix of the p eigenvector combinations with the largest inner products.}
  \item{nvec }{Vector containing the group sizes.}
  \item{reps }{Number of bootstrap replications to use.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
This is a utility function used by \code{\link{BVD}}, but can also be called directly if required.
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
Note that this function can currently handle only two groups of data.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{BVD}}
}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
