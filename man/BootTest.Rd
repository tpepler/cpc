\name{BootTest}
\alias{BootTest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Bootstrap hypothesis test (BootTest) method
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Identifies the number of common eigenvectors in several groups using the bootstrap hypothesis test (BootTest) method, adapted from Klingenberg and McIntyre (1998).
}
\usage{
BootTest(origdata, q = ncol(origdata[[1]]), reps = 1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{origdata }{List of the sample data sets.}
  \item{q }{Number of common eigenvectors to test for.}
  \item{reps}{Number of bootstrap replications to use.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
Tests the hypothesis, H_0: eigenvector pair are common, against the alternative, H_1: eigenvector pair are NOT common.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
Returns a data frame with the columns:
\item{Group1 }{Order of the eigenvector for Group 1.}
\item{Group2 }{Order of the eigenvector for Group 2.}
\item{vec.correlations }{Vector correlations of the eigenvector pairs.}
\item{p.values }{P-values for the null hypothesis of commonness of the eigenvector pairs.}
}
\references{
%% ~put references to the literature/web site here ~
Klingenberg, C. P. and McIntyre, G. S. (1998). Geometric morphometrics of developmental instability: Analyzing patterns of fluctuating asymmetry with Procrustes methods. Evolution, 52(5): 1363-1375.

Pepler, P.T. (2014). The identification and application of common principal components. PhD dissertation in the Department of Statistics and Actuarial Science, Stellenbosch University.
}
\author{
%%  ~~who you are~~
Theo Pepler
}
\note{
%%  ~~further notes~~
Note that this implementation of the BootTest method can currently handle only two groups of data.
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
BootTest(origdata = list(versicolor, virginica))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
