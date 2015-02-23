\name{nonnormaldata.sim}
\alias{nonnormaldata.sim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Simulate multivariate non-normal data
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Generates multivariate non-normal data for a specified covariance matrix structure.
}
\usage{
nonnormaldata.sim(Sigma, n = 100, df = rep(2, times = ncol(Sigma)))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Sigma }{Covariance matrix of the population from which to simulate the data.}
  \item{n }{Sample size.}
  \item{df }{Vector of chi square degrees of freedom, to control skewness of the variables (skew = sqrt(8/df)).}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
Simulates data from a multivariate distribution of which the variables are marginally distributed chi-squared with two degrees of freedom (default). This function can be useful in simulation studies when the purpose is to determine the effect of non-normality in the population on some statistical method, given a specific population covariance structure.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
Returns a matrix containing the simulated data.
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
}
\examples{
# Simulate 30 observations from a multivariate non-normally distributed
# population with the same covariance structure as the versicolor group
# in the Iris data set
data(iris)
versicolor <- iris[51:100, 1:4]
nonnormaldata.sim(Sigma = cov(versicolor), n = 30)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
