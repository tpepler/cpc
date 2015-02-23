\name{FG}
\alias{FG}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Flury-Gautschi (FG) algorithm}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
The Flury-Gautschi (FG) algorithm for the simultaneous diagonalisation of several covariance matrices.
}
\usage{
FG(covmats, nvec)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{covmats }{Array of covariance matrices to be simultaneously diagonalized, created by a command such as \code{covmats <- array(NA, dim = c(p, p, k))}.}
  \item{nvec }{Vector of sample sizes for the covariance matrices in \code{covmats}.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
A modal matrix (B) is calculated for the covariance matrices under the assumption of common principal components in the populations. The common eigenvectors in B are ordered according to the sums of the associated eigenvalues in the k groups, from the largest sum to the smallest.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
Returns a list with the values:
\item{B }{Orthogonal matrix containing the estimated common eigenvectors in the columns.}
\item{diagvals }{Matrix of which the columns contain the estimated eigenvalues of the k covariance matrices under the CPC model.}
}
\references{
%% ~put references to the literature/web site here ~
Flury, B. N. and Gautschi, W. (1986). An algorithm for simultaneous orthogonal transformation of several positive definite symmetric matrices to nearly diagonal form. SIAM Journal on Scientific and Statistical Computing, 7(1): 169-184.
}
\author{
Danie Uys, Theo Pepler
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{G.algorithm}}, \code{\link{stepwisecpc}}
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

cpc::FG(covmats = S, nvec = nvec)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
