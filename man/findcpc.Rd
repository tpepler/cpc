\name{findcpc}
\alias{findcpc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Find common principal components
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Provides descriptive measures to decide whether there are common eignevectors in several groups of data.
}
\usage{
findcpc(covmats, B = NULL, cutoff = 0.95, plotting = TRUE, main = "Vector correlations for the permutations")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{covmats}{
%%     ~~Describe \code{covmats} here~~
Array of covariance matrices of the k groups.
}
  \item{B}{
%%     ~~Describe \code{B} here~~
Modal matrix p x p matrix diagonalising the k covariance matrices simultaneously, estimated under the assumption of common eigenvectors in the groups.
}
  \item{cutoff}{
%%     ~~Describe \code{cutoff} here~~
Cut-off value to use in the vector correlation scree plot.
}
  \item{plotting}{
%%     ~~Describe \code{plotting} here~~
Should a scree plot of the vector correlations be constructed? (\code{TRUE}/\code{FALSE}).
}
  \item{main}{
%%     ~~Describe \code{main} here~~
Title of the scree plot, if \code{plotting = TRUE}.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
Identifies possibly common eigenvectors in k data groups by investigating the vectors correlations of all combinations of eigenvectors from the groups. These sets may be tested further for commonness.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
Produces a scree plot of the vector correlations (if \code{plotting = TRUE}) and returns a list with the values:
\item{all.correlations}{Summary of all eigenvector combinations from the k groups, and the geometric means of the vector correlations.}
\item{commonvec.order}{Order of the (possibly) common eigenvectors in the modal matrix (if an estimate was supplied).}
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
function (covmats, B = NULL, cutoff = 0.95, plotting = TRUE, 
    main = "Vector correlations for the permutations") 
{
    library(gtools)
    if (is.null(B)) {
        k <- dim(covmats)[3]
    }
    else {
        k <- dim(covmats)[3] + 1
    }
    p <- dim(covmats)[1]
    PCA.array <- array(NA, dim = c(p, p, k))
    if (is.null(B)) {
        for (i in 1:k) {
            PCA.array[, , i] <- eigen(covmats[, , i])$vectors
        }
    }
    else {
        PCA.array[, , 1] <- B
        for (i in 2:k) {
            PCA.array[, , i] <- eigen(covmats[, , (i - 1)])$vectors
        }
    }
    permsmat <- permutations(p, k, repeats.allowed = T)
    numperms <- nrow(permsmat)
    numcomparisons <- ncol(permsmat) - 1
    dotvec <- rep(0, times = numperms)
    for (i in 1:numperms) {
        temp <- rep(NA, times = numcomparisons)
        for (j in 1:numcomparisons) {
            testvecs <- cbind(PCA.array[, (permsmat[i, j]), j], 
                PCA.array[, (permsmat[i, j + 1]), j + 1])
            temp[j] <- abs(t(testvecs[, 1]) \%*\% testvecs[, 2])
        }
        dotvec[i] <- exp(mean(log(temp)))
    }
    resultmat <- cbind(permsmat, dot.products = dotvec)
    resultmat <- resultmat[order(resultmat[, "dot.products"], 
        decreasing = T), ]
    dot.sig <- resultmat[resultmat[, "dot.products"] >= cutoff, 
        , drop = F]
    plotnum <- max(c((p + 2), 20))
    plotnum <- min(c(plotnum, numperms))
    resultmat <- resultmat[c(1:plotnum), ]
    if (plotting) {
        resultmat <- resultmat[order(resultmat[, "dot.products"], 
            decreasing = F), ]
        par(pch = 20, xaxt = "n")
        plot(x = c(1:plotnum), y = resultmat[, "dot.products"], 
            type = "b", main = main, ylab = "Absolute vector correlation", 
            xlab = "", ylim = c(0, 1))
        abline(h = 1, lty = 3)
        par(adj = 0, srt = 90)
        for (i in 1:plotnum) {
            if (resultmat[i, "dot.products"] >= cutoff) {
                combtext <- toString(resultmat[i, 1:k])
                points(x = i, y = resultmat[i, "dot.products"], 
                  col = "red")
                text(x = i, y = resultmat[i, "dot.products"], 
                  labels = combtext, pos = 1, cex = 0.7, offset = 1)
            }
        }
        resultmat <- resultmat[order(resultmat[, "dot.products"], 
            decreasing = T), ]
    }
    if (is.null(B)) {
        colnames(resultmat) <- c("Group1", "Group2", "correlations")
    }
    else {
        colnames(resultmat) <- c("B", "Group1", "Group2", "correlations")
    }
    if (!is.null(B)) {
        commonvec.order <- resultmat[1:p, 1]
        commonvec.order <- append(unique(commonvec.order), sort(c(1:p)[-unique(commonvec.order)]))
    }
    else {
        commonvec.order <- NULL
    }
    par(adj = 0.5)
    return(list(all.correlations = resultmat, commonvec.order = commonvec.order))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
