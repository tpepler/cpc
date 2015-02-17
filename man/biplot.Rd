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
biplot(datamat, nvec, B, D3 = FALSE, varex = 1, plotvar = TRUE, main = "CPC biplot", col = c("blue", "red", "green", "orange", "brown", "purple"), radius = 0.1, lwd = 3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{datamat}{
%%     ~~Describe \code{datamat} here~~
Matrix of the data from the k groups combined.
}
  \item{nvec}{
%%     ~~Describe \code{nvec} here~~
Vector of the group sizes in the datamat matrix (in the correct order!).
}
  \item{B}{
%%     ~~Describe \code{B} here~~
Orthogonal projection matrix.
}
  \item{D3}{
%%     ~~Describe \code{D3} here~~
Should a 3-dimensional biplot be drawn? (\code{TRUE}/\code{FALSE})
}
  \item{varex}{
%%     ~~Describe \code{varex} here~~
Expansion factor for drawing the variables on the biplot.
}
  \item{plotvar}{
%%     ~~Describe \code{plotvar} here~~
Should the variables be drawn on the biplot? (\code{TRUE}/\code{FALSE})
}
  \item{main}{
%%     ~~Describe \code{main} here~~
Title of the biplot, to be passed to 'plot'.
}
  \item{col}{
%%     ~~Describe \code{col} here~~
Colors for the data points of the k groups.
}
  \item{radius}{
%%     ~~Describe \code{radius} here~~
Radius of the spheres in the case of a 3D biplot being drawn, to be passed to the function 'spheres3d'.
}
  \item{lwd}{
%%     ~~Describe \code{lwd} here~~
Line width of the variable vectors in the case of a 3D biplot being drawn, to be passed to the function 'lines3d'.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
This function may be used to produce a biplot in Euclidean space using any orthogonal projection matrix \code{B}. The obvious choice for \code{B} is the eigenvectors of \code{datamat}, producing the usual PCA biplot. Other choices for \code{B} includes the common eigenvector matrices estimated with the FG algorithm or stepwise CPC.
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (datamat, nvec, B, D3 = FALSE, varex = 1, plotvar = TRUE, 
    main = "CPC biplot", col = c("blue", "red", "green", "orange", 
        "brown", "purple"), radius = 0.1, lwd = 3) 
{
    for (i in 1:ncol(datamat)) {
        datamat[, i] <- datamat[, i] - mean(datamat[, i])
    }
    varnames <- colnames(datamat)
    p <- ncol(datamat)
    k <- length(nvec)
    plotpoints <- as.matrix(datamat) \%*\% B
    if (D3) {
        library(rgl)
        rgl.open()
        rgl.bg(color = "white", sphere = T)
        plot3d(x = plotpoints[, 1], y = plotpoints[, 2], z = plotpoints[, 
            3], xlab = "PC 1", ylab = "PC 2", zlab = "PC 3", 
            type = "n")
        decorate3d(main = main, xlab = NULL, ylab = NULL, zlab = NULL)
        begin <- 1
        for (i in 1:k) {
            end <- begin + nvec[i] - 1
            spheres3d(x = plotpoints[begin:end, 1], y = plotpoints[begin:end, 
                2], z = plotpoints[begin:end, 3], col = col[i], 
                radius = radius)
            begin <- end + 1
        }
        if (plotvar) {
            for (j in 1:p) {
                lines3d(x = c(0, B[j, 1] * varex), y = c(0, B[j, 
                  2] * varex), z = c(0, B[j, 3] * varex), lwd = lwd)
            }
            text3d(x = B[, 1] * varex, y = B[, 2] * varex, z = B[, 
                3] * varex, texts = varnames)
        }
    }
    else {
        library(MASS)
        par(pch = 20)
        eqscplot(x = plotpoints[, 1], y = plotpoints[, 2], type = "n", 
            xlab = "PC 1", ylab = "PC 2", main = main)
        begin <- 1
        for (i in 1:k) {
            end <- begin + nvec[i] - 1
            points(x = plotpoints[begin:end, 1], y = plotpoints[begin:end, 
                2], type = "p", col = col[i])
            begin <- end + 1
        }
        if (plotvar) {
            arrows(x0 = rep(0, times = p), x1 = B[, 1] * varex, 
                y0 = rep(0, times = p), y1 = B[, 2] * varex)
            text(x = B[, 1] * varex, y = B[, 2] * varex, pos = 3, 
                labels = varnames)
        }
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
