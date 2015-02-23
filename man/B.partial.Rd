\name{B.partial}
\alias{B.partial}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Eigenvector matrices under partial CPC model
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Estimates the eigenvector matrices under the partial CPC model, to diagonalise several covariance matrices.
}
\usage{
B.partial(covmats, nvec, B = cpc::FG(covmats = covmats, nvec = nvec)$B, commonvec.order, q)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{covmats }{Array of sample covariance matrices for the k groups.}
  \item{nvec }{Vector of sample sizes of the k groups.}
  \item{B }{Matrix of common eigenvectors estimated under the assumption of full CPC. Defaults to the modal matrix obtained with the FG algorithm.}
  \item{commonvec.order }{Vector containing the order of the common eigenvectors in B (with the q truly common eigenvectors in the first q positions).}
  \item{q }{Number of eigenvectors common to all k groups.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
Estimates the matrices of common (and non-common) eigenvectors for each of the groups, according to the method described in Flury (1988).
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
Returns an array containing the eigenvector matrices for the k groups, estimated under the CPC(q) model.
}
\references{
%% ~put references to the literature/web site here ~
Flury, B. (1988). Common Principal Components and Related Multivariate Models. Wiley.
}
\author{
Theo Pepler
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{FG}}, \code{\link{stepwisecpc}}
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

# Estimate the eigenvector matrices under the CPC(1) model
B.partial(covmats = S, nvec = nvec, q = 1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
