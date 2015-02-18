\name{stepwisecpc}
\alias{stepwisecpc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Stepwise common principal components (CPC)
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Estimates stepwise common principal components to simultaneously diagonalise several covariance matrices.
}
\usage{
stepwisecpc(covmats, nvec)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{covmats }{Array containing the covariance matrices for the groups, created with a command such as covmats <- array(NA, dim = c(p, p, k)), where p refers to the number of rows/columns of each covariance matrix, and k is the number of groups (or covariance matrices).}
  \item{nvec }{Vector containing the sample sizes of the k groups.}
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
Returns a list with the values:
\item{B }{Orthogonal matrix containing the estimated common eigenvectors in the columns.}
\item{eigenvals }{Matrix of which the columns contain the estimated eigenvalues of the k covariance matrices under the CPC model.}
}
\references{
%% ~put references to the literature/web site here ~
Trendafilov, N. T. (2010). Stepwise estimation of common principal components. Computational Statistics and Data Analysis, 54(12): 3446-3457.
}
\author{
Sarel Steel, Theo Pepler
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{FG}}
}
\examples{
# Estimate the modal matrix for the versicolor and virginica groups
# under the common eigenvector assumption.

data(iris)
versicolor <- iris[51:100, 1:4]
virginica <- iris[101:150, 1:4]

# Create array containing the covariance matrices
S <- array(NA, dim = c(4, 4, 2))
S[, , 1] <- cov(versicolor)
S[, , 2] <- cov(virginica)

stepwisecpc(covmats = S, nvec = c(nrow(versicolor), nrow(virginica)))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
