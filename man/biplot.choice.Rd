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
  \item{datalist}{
%%     ~~Describe \code{datalist} here~~
List of the data for which the biplot is to be constructed, created with a command such as \code{list(group1, group2, ...)}.
}
  \item{rdim}{
%%     ~~Describe \code{rdim} here~~
Number of dimensions of the biplot representation.
}
  \item{add.projectmats}{
%%     ~~Describe \code{add.projectmats} here~~
(Optional) Additional orthogonal projection matrices for which to compute the biplot quality measures [NOT YET IMPLEMENTED].
}
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
\item{Overall}{Overall quality of the biplot display.}
\item{Within}{Mean quality of the within-group variation as displayed in the biplot.}
\item{Between}{Quality of the between-group variation as displayed in the biplot.}
\item{Adequacy}{Median adequacy of the variables.}
\item{MSPE}{Mean of the mean standard predictive errors (MSPE) of the variables.}
\item{Sample predictivities}{Mean predictivity of the observations as displayed in the biplot.}
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
See \code{\link{biplot.measures}} for more details on the different quality measures for principal component biplots of grouped data.
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (datalist, rdim, add.projectmats = NULL) 
{
    k <- length(datalist)
    p <- ncol(datalist[[1]])
    nvec <- rep(NA, times = k)
    dfpooled <- 0
    SSpooled <- 0
    for (i in 1:k) {
        nvec[i] <- nrow(datalist[[i]])
        dfpooled <- dfpooled + nvec[i] - 1
        SSpooled <- SSpooled + cov(as.matrix(datalist[[i]])) * 
            (nvec[i] - 1)
    }
    Sp <- SSpooled/dfpooled
    Ep <- eigen(Sp)$vectors
    pooledcov.output <- biplot.measures(datalist = datalist, 
        projectmat = Ep, rdim = rdim)
    datamat <- NULL
    for (i in 1:k) {
        datamat <- rbind(datamat, as.matrix(datalist[[i]]))
    }
    E <- eigen(cov(datamat))$vectors
    pooleddata.output <- biplot.measures(datalist = datalist, 
        projectmat = E, rdim = rdim)
    S <- array(NA, dim = c(p, p, k))
    for (i in 1:k) {
        S[, , i] <- cov(as.matrix(datalist[[i]]))
    }
    B.flury <- FG(covmats = S, nvec = nvec)$B
    flury.output <- biplot.measures(datalist = datalist, projectmat = B.flury, 
        rdim = rdim)
    B.stepwise <- stepwisecpc(covmats = S, nvec = nvec)$B
    stepwise.output <- biplot.measures(datalist = datalist, projectmat = B.stepwise, 
        rdim = rdim)
    library(JADE)
    B.jade <- rjd(X = S)$V
    jade.output <- biplot.measures(datalist = datalist, projectmat = B.jade, 
        rdim = rdim)
    pooledcov <- c(pooledcov.output$overall.quality, pooledcov.output$within.quality.mean, 
        pooledcov.output$between.quality, pooledcov.output$adequacies.median, 
        pooledcov.output$mspe.mean, pooledcov.output$sample.predictivities.mean)
    pooleddata <- c(pooleddata.output$overall.quality, pooleddata.output$within.quality.mean, 
        pooleddata.output$between.quality, pooleddata.output$adequacies.median, 
        pooleddata.output$mspe.mean, pooleddata.output$sample.predictivities.mean)
    flury <- c(flury.output$overall.quality, flury.output$within.quality.mean, 
        flury.output$between.quality, flury.output$adequacies.median, 
        flury.output$mspe.mean, flury.output$sample.predictivities.mean)
    stepwise <- c(stepwise.output$overall.quality, stepwise.output$within.quality.mean, 
        stepwise.output$between.quality, stepwise.output$adequacies.median, 
        stepwise.output$mspe.mean, stepwise.output$sample.predictivities.mean)
    jade <- c(jade.output$overall.quality, jade.output$within.quality.mean, 
        jade.output$between.quality, jade.output$adequacies.median, 
        jade.output$mspe.mean, jade.output$sample.predictivities.mean)
    resultsmat <- rbind(pooledcov, pooleddata, flury, stepwise, 
        jade)
    rownames(resultsmat) <- c("Pooled S", "Pooled data", "Flury", 
        "Stepwise CPC", "JADE")
    colnames(resultsmat) <- c("Overall", "Within", "Between", 
        "Adequacy", "MSPE", "Sample predictivities")
    return(resultsmat)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
