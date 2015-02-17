\name{RVC}
\alias{RVC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Random vector correlations (RVC) method
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Identifies the number of common eigenvectors in several groups using the random vector correlations (RVC) method, adapted from Klingenberg and McIntyre (1998).
}
\usage{
RVC(covmats, reps = 100000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{covmats}{
%%     ~~Describe \code{covmats} here~~
Array of covariance matrices for the k groups.
}
  \item{reps}{
%%     ~~Describe \code{reps} here~~
Number of randomisations to use.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
Vector correlations between sample eigenvectors are compared to a distribution of vector correlations between of pairs of vectors randomly generated on the unit sphere, in order to determine whether the associated eigenvectors of the population covariance matrices are common.

%Tests the hypothesis:
%H_0: eigenvector pair are NOT common
%H_1: eigenvector pair are common
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
Returns a data frame with the columns:
\item{commonvec.order}{Order of the eigenvectors in the two groups.}
\item{vec.correlations}{Vector correlations of the eigenvector pairs.}
\item{p.values}{P values for the null hypothesis of commonness of the eigenvector pairs.}
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
Note that this implementation of the RVC method can currently handle only two groups of data.
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

# Create array containing the covariance matrices
S <- array(NA, dim = c(4, 4, 2))
S[, , 1] <- cov(versicolor)
S[, , 2] <- cov(virginica)

RVC(covmats = S)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
