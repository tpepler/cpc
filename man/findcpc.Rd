\name{findcpc}
\alias{findcpc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Find common principal components
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Provides descriptive measures to decide whether there are common eigenvectors in several covariance matrices.
}
\usage{
findcpc(covmats, B = NULL, cutoff = 0.95, plotting = TRUE, main = "Vector correlations for the permutations")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{covmats }{Array of covariance matrices of the k groups.}
  \item{B }{Modal matrix p x p matrix diagonalising the k covariance matrices simultaneously, estimated under the assumption of common eigenvectors in the groups. Can be estimated using simultaneous diagonalisation algorithms such as the Flury-Gautschi (implemented in \code{\link{FG}} or the stepwise CPC (implemented in \code{\link{stepwisecpc}}) algorithms.}
  \item{cutoff }{Cut-off value to use in the vector correlation scree plot.}
  \item{plotting }{Logical, indicating whether a scree plot of the vector correlations should be constructed (default = TRUE).}
  \item{main }{Title of the scree plot, if \code{plotting = TRUE}.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
Identifies possibly common eigenvectors in k data groups by investigating the vectors correlations of all combinations of eigenvectors from the groups. These sets may be tested further for commonness.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
Produces a scree plot of the vector correlations (if \code{plotting = TRUE}) and returns a list with the values:
\item{all.correlations }{Summary of all eigenvector combinations from the k groups, and the geometric means of the vector correlations.}
\item{commonvec.order }{Order of the (possibly) common eigenvectors in the modal matrix (if an estimate was supplied).}
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
\code{\link{ensemble.test}}, \code{\link{flury.test}}
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

findcpc(covmats = S)

# Estimate the modal matrix with the FG algorithm
nvec <- c(nrow(versicolor), nrow(virginica))
B <- cpc::FG(covmats = S, nvec = nvec)$B
findcpc(covmats = S, B = B)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
