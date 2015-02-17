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
  \item{datalist}{
%%     ~~Describe \code{datalist} here~~
List of the data for which the biplot is to be constructed, created with a command such as \code{list(group1, group2, ...)}.
}
  \item{projectmat}{
%%     ~~Describe \code{projectmat} here~~
Orthogonal projection matrix used in constructing the biplot.
}
  \item{rdim}{
%%     ~~Describe \code{rdim} here~~
Number of dimensions of the biplot representation.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
For the plain PCA biplot, use the eigenvectors of the pooled data (not centred) as the projection matrix.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
Returns a list with the values:
\item{overall.quality}{Overall quality of the biplot display, a scalar value in the[0-1] range.}
\item{within.quality}{Quality of the within-group variation per group, a vector of values in the [0-1] range.}
\item{within.quality.mean}{Mean quality of the within-group variation displayed in the biplot, a scalar value in the [0-1] range.}
\item{between.quality}{Quality of the between-group variation as displayed in the biplot, a scalar value in the [0-1] range.}
\item{adequacies}{Vector containing the adequacies of the variables (values in the [0-1] range).}
\item{adequacies.median}{Median of the variable adequacies.}
\item{axis.predictivities}{Vector containing the axis predictivities (values in the [0-1] range, but see the notes below).}
\item{axis.predictivities.mean}{Mean of the axis predicitivities.}
\item{sample.predictivities}{Vector containing the sample predictivities (values in the [0-1] range).}
\item{sample.predictivities.mean}{Mean of the sample predictivities.}
\item{mspe}{Vector containing the mean standard predictive errors (MSPE) of the variables.}
\item{mspe.mean}{Mean of the MSPE values.}
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
The median (instead of the mean) adequacy of the variables is calculated, as the mean adequacy of the variables will always be equal to r/p (r = number of dimensions; p = number of variables), and is therefore uninformative as a quality measure.

The axis predictivities of the variables are only valid if the Type B orthogonality condition holds. It is thus a valid measure for the plain PCA biplot, but not for other types of principal component biplots.
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
function (datalist, projectmat, rdim) 
{
    k <- length(datalist)
    p <- ncol(datalist[[1]])
    if (rdim > p) {
        cat("Number of biplot dimensions cannot be larger than number of variables in data!\n")
        return()
    }
    varnames <- colnames(datalist[[1]])
    X <- NULL
    for (i in 1:k) {
        X <- rbind(X, as.matrix(datalist[[i]]))
    }
    n <- nrow(X)
    Xmean <- apply(X, 2, mean)
    X <- t(t(X) - Xmean)
    Y <- X \%*\% projectmat[, 1:rdim]
    totalvar <- sum(diag(t(X) \%*\% X))
    fittedvar <- sum(diag(t(Y) \%*\% Y))
    overall.quality <- fittedvar/totalvar
    withingroup.total <- rep(NA, times = k)
    withingroup.fitted <- rep(NA, times = k)
    within.total <- 0
    within.fitted <- 0
    between.total <- 0
    between.fitted <- 0
    Xmean.fit <- t(projectmat[, 1:rdim]) \%*\% Xmean
    for (i in 1:k) {
        groupdata <- as.matrix(datalist[[i]])
        groupmean <- apply(groupdata, 2, mean)
        groupdata <- t(t(groupdata) - groupmean)
        fitteddata <- groupdata \%*\% projectmat[, 1:rdim]
        withingroup.total[i] <- sum(diag(t(groupdata) \%*\% groupdata))
        within.total <- within.total + withingroup.total[i]
        withingroup.fitted[i] <- sum(diag(t(fitteddata) \%*\% fitteddata))
        within.fitted <- within.fitted + withingroup.fitted[i]
        between.total <- between.total + t(groupmean - Xmean) \%*\% 
            (groupmean - Xmean)
        groupfitted.mean <- t(projectmat[, 1:rdim]) \%*\% groupmean
        between.fitted <- between.fitted + t(groupfitted.mean - 
            Xmean.fit) \%*\% (groupfitted.mean - Xmean.fit)
    }
    within.quality <- withingroup.fitted/withingroup.total
    within.quality <- matrix(within.quality, ncol = k, byrow = T)
    colnames(within.quality) <- paste("Group ", c(1:k), sep = "")
    within.quality.mean <- within.fitted/within.total
    between.quality <- as.numeric(between.fitted/between.total)
    adequacies <- diag(projectmat[, 1:rdim] \%*\% t(projectmat[, 
        1:rdim]))
    adequacies.median <- median(adequacies)
    adequacies <- matrix(adequacies, ncol = p, byrow = T)
    colnames(adequacies) <- varnames
    J <- diag(c(rep(1, times = rdim), rep(0, times = (p - rdim))))
    X.fitted <- X \%*\% projectmat \%*\% J \%*\% t(projectmat)
    axis.predictivities <- diag(diag(diag(t(X.fitted) \%*\% X.fitted)) \%*\% 
        solve(diag(diag(t(X) \%*\% X))))
    axis.predictivities.mean <- sum(diag(diag(diag(t(X.fitted) \%*\% 
        X.fitted)) \%*\% solve(diag(diag(t(X) \%*\% X)))))/p
    axis.predictivities <- matrix(axis.predictivities, ncol = p, 
        byrow = T)
    colnames(axis.predictivities) <- varnames
    sample.predictivities <- diag(diag(diag(X.fitted \%*\% t(X.fitted))) \%*\% 
        solve(diag(diag(X \%*\% t(X)))))
    sample.predictivities.mean <- sum(diag(diag(diag(X.fitted \%*\% 
        t(X.fitted))) \%*\% solve(diag(diag(X \%*\% t(X))))))/nrow(X)
    Xsd <- apply(X, 2, sd)
    onevec <- matrix(1, nrow = n, ncol = 1)
    mspe <- rep(NA, times = p)
    for (j in 1:p) {
        mspe[j] <- (t(onevec) \%*\% abs(X[, j] - X \%*\% projectmat[, 
            1:rdim] \%*\% t(projectmat[j, 1:rdim, drop = F])))/(n * 
            Xsd[j])
    }
    mspe.mean <- mean(mspe)
    return(list(overall.quality = overall.quality, within.quality = within.quality, 
        within.quality.mean = within.quality.mean, between.quality = between.quality, 
        adequacies = adequacies, adequacies.median = adequacies.median, 
        axis.predictivities = axis.predictivities, axis.predictivities.mean = axis.predictivities.mean, 
        sample.predictivities = sample.predictivities, sample.predictivities.mean = sample.predictivities.mean, 
        mspe = mspe, mspe.mean = mspe.mean))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
