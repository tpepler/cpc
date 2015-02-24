\name{biplot.measures}
\alias{biplot.measures}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Biplot quality measures
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Calculates a number of quality measures for principal component biplots of grouped data.
}
\usage{
biplot.measures(datalist, projectmat, rdim)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{datalist }{List of the data for which the biplot is to be constructed, created with a command such as \code{list(group1, group2, ...)}.}
  \item{projectmat }{Orthogonal projection matrix used in constructing the biplot.}
  \item{rdim }{Number of dimensions of the biplot representation.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
For the plain PCA biplot, use the eigenvectors of the pooled data (not centred per group before pooling) as the projection matrix.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
Returns a list with the values:
\item{overall.quality }{Overall quality of the biplot display, a scalar value in the[0-1] range.}
\item{within.quality }{Quality of the within-group variation per group, a vector of values in the [0-1] range.}
\item{within.quality.mean }{Mean quality of the within-group variation displayed in the biplot, a scalar value in the [0-1] range.}
\item{between.quality }{Quality of the between-group variation as displayed in the biplot, a scalar value in the [0-1] range.}
\item{adequacies }{Vector containing the adequacies of the variables (values in the [0-1] range).}
\item{adequacies.median }{Median of the variable adequacies.}
\item{axis.predictivities }{Vector containing the axis predictivities (values in the [0-1] range, but see the notes below).}
\item{axis.predictivities.mean }{Mean of the axis predicitivities.}
\item{sample.predictivities }{Vector containing the sample predictivities (values in the [0-1] range).}
\item{sample.predictivities.mean }{Mean of the sample predictivities.}
\item{mspe }{Vector containing the mean standard predictive errors (MSPE) of the variables.}
\item{mspe.mean }{Mean of the MSPE values.}
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
The median (instead of the mean) adequacy of the variables is calculated, as the mean adequacy of the variables will always be equal to r/p (r = number of dimensions; p = number of variables), and is therefore uninformative as a quality measure.

The axis predictivities of the variables are only valid if the Type B orthogonality condition holds. It is thus a valid measure for the plain PCA biplot, but not for other types of principal component biplots.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{biplot.choice}}
}
\examples{
# The Iris data
data(iris)
setosa <- iris[1:50, 1:4]
versicolor <- iris[51:100, 1:4]
virginica <- iris[101:150, 1:4]

project.matrix <- eigen(cov(rbind(setosa, versicolor, virginica)))$vectors

# For a 2-dimensional biplot
biplot.measures(datalist = list(setosa, versicolor, virginica), projectmat = project.matrix, rdim = 2)

# For a 3-dimensional biplot
biplot.measures(datalist = list(setosa, versicolor, virginica), projectmat = project.matrix, rdim = 3)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
