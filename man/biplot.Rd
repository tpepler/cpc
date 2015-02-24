\name{biplot}
\alias{biplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Biplot
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Produces a 2- or 3-dimensional biplot of the data according to the orthogonal projection matrix supplied.
}
\usage{
biplot(datalist, B, D3 = FALSE, varex = 1, plotvar = TRUE, main = "CPC biplot", col = c("blue", "red", "green", "orange", "brown", "purple"), radius = 0.1, lwd = 3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{datalist }{List of the data from the k groups.}
  \item{B }{Orthogonal projection matrix.}
  \item{D3 }{Logical, indicating whteher a 3-dimensional biplot should be produced (default = FALSE).}
  \item{varex }{Expansion factor for drawing the variables on the biplot.}
  \item{plotvar }{Logical, indicating whether the variables should be drawn on the biplot (default = TRUE).}
  \item{main }{Title of the biplot, to be passed to \code{\link{plot}}.}
  \item{col }{Colors for the data points of the k groups.}
  \item{radius }{Radius of the spheres in the case of a 3D biplot being drawn, to be passed to \code{\link{spheres3d}}.}
  \item{lwd }{Line width of the variable vectors in the case of a 3D biplot being drawn, to be passed to \code{\link{lines3d}}.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
This function may be used to produce a biplot in Euclidean space using any orthogonal projection matrix \code{B}. The obvious choice for \code{B} is the eigenvectors of the pooled data, producing the usual PCA biplot. Other choices for \code{B} includes the common eigenvector matrices estimated with the FG algorithm or stepwise CPC.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
Produces a biplot in a graphics device window.
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
\code{\link{biplot.measures}}, \code{\link{biplot.choice}}
}
\examples{
# The Iris data
data(iris)
setosa <- iris[1:50, 1:4]
versicolor <- iris[51:100, 1:4]
virginica <- iris[101:150, 1:4]
project.matrix <- eigen(cov(rbind(setosa, versicolor, virginica)))$vectors

# For a 2-dimensional biplot
cpc::biplot(datalist = list(setosa, versicolor, virginica), B = project.matrix, varex = 3, main = "CPC biplot: Iris data")

# For a 3-dimensional biplot
cpc::biplot(datalist = list(setosa, versicolor, virginica), B = project.matrix, D3 = TRUE, main = "CPC biplot: Iris data")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
