\name{ensemble.test}
\alias{ensemble.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Ensemble test for the identification of common eigenvectors
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Identifies the number of common eigenvectors in several groups using the ensemble test.
}
\usage{
ensemble.test(origdata, standardize = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{origdata }{List of the sample data sets.}
  \item{standardize }{Logical, indicating whether the data columns should be standardized (mean=0, stdev=1) before performing the ensemble test (default = FALSE).}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
Ensemble method to identify common eigenvectors in k groups: majority vote on number of common eigenvectors from Flury's AIC, Bootstrap Vector Correlation Distribution (BVD), Bootstrap Confidence Regions (BCR), Random Vector Correlations (RVC) and Bootstrap hypothesis test (BootTest) methods.  
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
Returns a list with the following components:
\item{Results }{Row 1: order of common eigenvectors in B; Row 2-5: results from AIC, BVD, BCR and RVC tests (1 = eigenvector common); Row 6: ensemble test common eigenvector indicator (1 = eigenvector common).}
\item{commonvecs }{Positions of the common eigenvectors in the modal matrix.}
\item{commonvecmat }{The estimated common eigenvectors, extracted from the modal matrix.}
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
Note that this implementation of the Ensemble test can currently handle only two groups of data.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{flury.AIC}}, \code{\link{BVD}}, \code{\link{BCR}}, \code{\link{RVC}}, and \code{\link{BootTest}}
}
\examples{
# Determine number of common eigenvectors in the covariance matrices of the 
# versicolor and virginica groups

data(iris)
versicolor <- iris[51:100, 1:4]
virginica <- iris[101:150, 1:4]
ensemble.test(origdata = list(versicolor, virginica))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
