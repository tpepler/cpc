\name{cpc.test}
\alias{cpc.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Likelihood ratio test of common principal components in the covariance matrices of several groups
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Calculates the likelihood ratio statistic and its degrees of freedom for the hypothesis of common eigenvectors in the k groups against the alternative of unrelated covariance matrices.
}
\usage{
cpc.test(covmats, nvec, B = cpc::FG(covmats = covmats, nvec = nvec)$B)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{covmats }{Array of covariance matrices.}
  \item{nvec }{Vector of sample sizes of the k groups.}
  \item{B }{Modal matrix simultaneously diagonalising the covariance matrices, estimated under the assumption of common eigenvectors in the k groups. Can be estimated using simultaneous diagonalisation algorithms such as the Flury-Gautschi (implemented in \code{\link{FG}} or the stepwise CPC (implemented in \code{\link{stepwisecpc}}) algorithms.}
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
Returns a list with the following:
\item{chi.square }{The likelihood ratio test statistic.}
\item{df }{Degrees of freedom of the test statistic under the null hypothesis.}
\item{covmats.cpc }{Estimated covariance matrices under the null hypothesis model.}
}
\references{
%% ~put references to the literature/web site here ~
Flury, B. (1988). Common Principal Components and Related Multivariate Models. Wiley.
}
\author{
%%  ~~who you are~~
Theo Pepler
}
\note{
%%  ~~further notes~~
This test is based on the assumption that the populations from which the data originated are distributed multivariate normal.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{FG}}, \code{\link{flury.test}}, \code{\link{equal.test}}, \code{\link{prop.test}} and \code{\link{cpcq.test}}
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

cpc.test(covmats = S, nvec = nvec)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
