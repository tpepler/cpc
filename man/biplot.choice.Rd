\name{biplot.choice}
\alias{biplot.choice}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Biplot choice tool
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Provides quality measures for different types of principal component biplots of grouped data, to select the best method for constructing a principal component biplot for any specific data set.
}
\usage{
biplot.choice(datalist, rdim, add.projectmats = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{datalist }{List of the data for which the biplot is to be constructed, created with a command such as \code{list(group1, group2, ...)}.}
  \item{rdim }{Number of dimensions of the biplot representation.}
  \item{add.projectmats }{(Optional) Additional orthogonal projection matrices for which to compute the biplot quality measures [NOT YET IMPLEMENTED].}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
This function calculates a number of quality measures relevant in judging the appropriateness of a principal component biplot for grouped data. The user would then choose the method leading to the best biplot for the purpose of interest.

Quality measures are automatically calculated for biplots constructed from the following orthogonal projection matrices:

'Pooled S' -- eigenvectors of the pooled data (groups uncentred)

'Pooled data' -- eigenvectors obtained from the pooled covariance matrix

'Flury' -- common eigenvector matrix estimated with the FG algorithm

'Stepwise CPC' -- common eigenvector matrix estimated with the stepwise CPC algorithm

'JADE' -- common eigenvector matrix estimated with the \code{rjd} function in the \code{JADE} package
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
Returns a matrix of results with the columns:
\item{Overall }{Overall quality of the biplot display.}
\item{Within }{Mean quality of the within-group variation as displayed in the biplot.}
\item{Between }{Quality of the between-group variation as displayed in the biplot.}
\item{Adequacy }{Median adequacy of the variables.}
\item{MSPE }{Mean of the mean standard predictive errors (MSPE) of the variables.}
\item{Sample predictivities }{Mean predictivity of the observations as displayed in the biplot.}
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
See \code{\link{biplot.measures}} for more details on the different quality measures for principal component biplots of grouped data.
}
\examples{
# The Iris data
data(iris)
setosa <- iris[1:50, 1:4]
versicolor <- iris[51:100, 1:4]
virginica <- iris[101:150, 1:4]

# For a 2-dimensional biplot
biplot.choice(datalist = list(setosa, versicolor, virginica), rdim = 2)

# For a 3-dimensional biplot
biplot.choice(datalist = list(setosa, versicolor, virginica), rdim = 3)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
