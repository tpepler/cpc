\name{BVD}
\alias{BVD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Bootstrap vector correlation distribution (BVD) method
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Identifies the number of common eigenvectors in several groups using the bootstrap vector correlation distributions (BVD) method.
}
\usage{
BVD(origdata, reps = 1000)
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
Identifies the number of common eigenvectors using the BVD 1a, BVD 1b, BVD 1c, BVD 2a, BVD 2b, BVD 2c and BVD 2d methods.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
Returns a matrix with the following rows:
\item{Common eigenvector}{Number of the eigenvector in the modal matrix.}
\item{BVD 1a}{Vector indicating the common eigenvectors (=1) according to this method.}
\item{BVD 1b}{Vector indicating the common eigenvectors (=1) according to this method.}
\item{BVD 1c}{Vector indicating the common eigenvectors (=1) according to this method.}
\item{BVD 2a}{Vector indicating the common eigenvectors (=1) according to this method.}
\item{BVD 2b}{Vector indicating the common eigenvectors (=1) according to this method.}
\item{BVD 2c}{Vector indicating the common eigenvectors (=1) according to this method.}
\item{BVD 2d}{Vector indicating the common eigenvectors (=1) according to this method.}
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
Note that this implementation of the BVD method can currently handle only two groups of data.
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
BVD(origdata = list(versicolor, virginica))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
