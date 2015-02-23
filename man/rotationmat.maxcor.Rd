\name{rotationmat.maxcor}
\alias{rotationmat.maxcor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Maximum correlation rotation matrix
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Finds a rotation matrix with maximum correlation between the variables.
}
\usage{
rotationmat.maxcor(p)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p }{Number of variables (rows/columns) required for the rotation matrix.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
Finding rotation matrices (approximately) maximising the correlations between the variables can be useful in some simulation study settings, for example when the purpose is to study the effect of high correlations on the performance of some statistical method.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
Returns the p x p square rotation matrix.
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
rotationmat.maxcor(p = 5)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
