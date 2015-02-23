\name{alpha.crossvalid}
\alias{alpha.crossvalid}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Crossvalidation method to estimate shrinkage intensity
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Estimates alpha shrinkage intensity parameter by the method proposed in Pepler (2014), for improved estimation of population covariance matrices.
}
\usage{
alpha.crossvalid(datamat, B, reps = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{datamat }{Matrix containing sample data for the ith group.}
  \item{B }{Matrix of estimated common (and possibly non-common) eigenvectors. Can be estimated using simultaneous diagonalisation algorithms such as the Flury-Gautschi (implemented in \code{\link{FG}} or the stepwise CPC (implemented in \code{\link{stepwisecpc}}) algorithms.}
  \item{reps }{Number of replications to use in cross-validation.}
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
Returns the estimated shrinkage intensity (scalar), a value between 0 and 1.
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
# Versicolor and virginica groups of the Iris data
data(iris)
versicolor <- iris[51:100, 1:4]
virginica <- iris[101:150, 1:4]

# Create array containing the two covariance matrices
S <- array(NA, c(4, 4, 2))
S[, , 1] <- cov(versicolor)
S[, , 2] <- cov(virginica)
nvec <- c(nrow(versicolor), nrow(virginica))

# Estimate the modal matrix using the FG algorithm
B <- FG(covmats = S, nvec = nvec)$B

# Estimate optimal shrinkage intensity for the versicolor covariance matrix
alpha.crossvalid(datamat = versicolor, B = B, reps = 1000)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
