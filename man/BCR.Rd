\name{BCR}
\alias{BCR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Bootstrap confidence regions (BCR) method
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Identifies the number of common eigenvectors in several groups using the bootstrap confidence regions (BCR) method.
}
\usage{
BCR(origdata, reps = 1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{origdata}{
%%     ~~Describe \code{origdata} here~~
List of the sample data sets.
}
  \item{reps}{
%%     ~~Describe \code{reps} here~~
Number of bootstrap replications to use.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
Calculates 95\% bootstrap confidence regions for eigenvector pairs; if the regions overlap, the eigenvectors are considered to be common.

%Tests the hypothesis:
%H_0: the pair of eigenvectors are common
%H_1: the pair of eigenvectors are NOT common
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
Returns a data frame with the columns:
\item{Common eigenvector}{Order of the eigenvectors in the model matrix.}
\item{common.ind}{Indicator of eigenvector commonness (1 = common eigenvector).}
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
Note that this implementation of the BCR method can currently handle only two groups of data.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{ensemble.test}}
}
\examples{
# Determine number of common eigenvectors in the covariance matrices of the 
# versicolor and virginica groups

data(iris)
versicolor <- iris[51:100, 1:4]
virginica <- iris[101:150, 1:4]
BCR(origdata = list(versicolor, virginica))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
