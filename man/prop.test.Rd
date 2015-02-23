\name{prop.test}
\alias{prop.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Likelihood ratio test of proportionality of several covariance matrices
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Calculates the likelihood ratio statistic and its degrees of freedom for the hypothesis of proportional covariance matrices against the alternative of unrelated covariance matrices.
}
\usage{
prop.test(covmats, nvec)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{covmats }{Array of covariance matrices.}
  \item{nvec }{Vector of sample sizes of the k groups.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
This is an implementation of the algorithm described in Flury (1988).
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
\item{covmats.prop }{Estimated covariance matrices under the null hypothesis model.}
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
\code{\link{flury.test}}, \code{\link{equal.test}}, \code{\link{cpc.test}} and \code{\link{cpcq.test}}
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

cpc::prop.test(covmats = S, nvec = nvec)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line