\name{flury.phi}
\alias{flury.phi}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Flury's phi measure
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Calculates phi measure of goodness of diagonalisation as proposed by Flury (1988).
}
\usage{
flury.phi(datamat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{datamat}{A square symmetric matrix.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
A smaller phi value points to better diagonalisation, and a larger value to worse diagonalisation. Phi can be a minimum of one, if \code{datamat} is perfectly diagonal.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
Returns the value (scalar) of the phi measure for the matrix.
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
Bmats <- B.partial(covmats = S, nvec = nvec, q = 1)

# Check the goodness of diagonalisation for the versicolor
# sample covariance matrix
S1.diag <- t(Bmats[, , 1]) \%*\% S[, , 1] \%*\% Bmats[, , 1]
flury.phi(S1.diag)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
