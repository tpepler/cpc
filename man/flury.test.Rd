\name{flury.test}
\alias{flury.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Flury's tests for common covariance structures
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Calculates the AIC and chi square statistics for all the applicable models in Flury's hierarchy of covariance structures for several covariance matrices.
}
\usage{
flury.test(covmats, nvec, B = FG(covmats = covmats, nvec = nvec)$B, qmax = NULL, commonvec.order = findcpc(covmats = covmats, B = B, plotting = F)$commonvec.order)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{covmats }{Array of covariance matrices to be tested, created by a command such as \code{covmats <- array(NA, dim = c(p, p, k))}.}
  \item{nvec }{Vector of sample sizes of the k groups.}
  \item{B }{Modal matrix, i.e. the orthogonal p x p matrix diagonalising the k covariance matrices simultaneously under the assumption of common principal components. Can be estimated using simultaneous diagonalisation algorithms such as the Flury-Gautschi (implemented in \code{\link{FG}} or the stepwise CPC (implemented in \code{\link{stepwisecpc}}) algorithms.}
  \item{qmax }{Maximum number of common eigenvectors for which to estimate the partial CPC model statistics.}
  \item{commonvec.order }{Vector indicating the order of the most likely candidates of common eigenvectors -- from 1 (most likely) to p (least likely).}
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
Returns a data frame with the following columns, containing a value for each of the models in Flury's hierarchy:
\item{Model }{Name of the model in Flury's hierarchy.}
\item{Chi.square }{Chi square statistic for the model.}
\item{DF }{Degrees of freedom associated with the chi square statistic.}
\item{Chi2.div.df }{The chi square statistic divided by its degrees of freedom.}
\item{AIC }{Akaike Information Criterion for the model.}
\item{No.of.CPCs }{Number of common eigenvectors for the model.}
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
Flury (1988) suggested picking the model with the smallest AIC value. Alternatively one may select the model for which the chi square statistic divided by its degrees of freedom is closest to one, although Flury (and others) cautions against this practice as the chi square statistics for the models are probably not independent of one another.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{flury.AIC}}, \code{\link{equal.test}}, \code{\link{prop.test}}, \code{\link{cpc.test}} and \code{\link{cpcq.test}}
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

flury.test(covmats = S, nvec = nvec)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
