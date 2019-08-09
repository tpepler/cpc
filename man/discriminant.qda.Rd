\name{discriminant.qda}
\alias{discriminant.qda}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Discriminant analysis under various assumed covariance matrix structures
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Performs discriminant analysis for two groups, under the assumption of equal, CPC or unrelated population covariance matrices.
}
\usage{
discriminant.qda(origdata, group, method = c("unbiased", "pooled", "cpc", "fullcpccrossvalid"), B = NULL, standardize = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{origdata }{Matrix containing the sample data for two groups.}
  \item{group }{Vector (with values 1 and 2) indicating group membership for the rows in \code{origdata}.}
  \item{method }{Options: \code{unbiased} (default), \code{cpc}, \code{fullcpccrossvalid}, \code{pooled}. Each method estimates the covariance matrices of the two groups in a different way (under different assumptions). See Pepler (2014) for details.}
  \item{B }{Modal matrix, if available.}
  \item{standardize }{Logical, indicating whether the covariance matrices of the two groups should be calculated on standardised data before performing the discriminant analysis (default = FALSE). Discriminant analysis on the standardised covariance matrices is equivalent to discriminant analysis on the correlation matrices of the two groups.}
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
Returns a list with the following components:
\item{inputdata }{The observations in \code{origdata}, together with an indication of the group membership and the predicted group membership.}
\item{misclassrate }{The misclasssification error rate.}
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
# Discriminant analysis to distinguish between the versicolor 
# and virginica groups, under the CPC assumption.
data(iris)
versicolor <- iris[51:100, 1:4]
virginica <- iris[101:150, 1:4]
S <- array(NA, dim = c(4, 4, 2))
S[, , 1] <- cov(versicolor)
S[, , 2] <- cov(virginica)
nvec <- c(nrow(versicolor), nrow(virginica))
B <- cpc::FG(covmats = S, nvec = nvec)$B
discriminant.qda(origdata = rbind(versicolor, virginica), group = c(rep(1, times = nrow(versicolor)), rep(2, times = nrow(virginica))), method = "cpc", B = B)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
