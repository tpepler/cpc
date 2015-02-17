\name{BootTest}
\alias{BootTest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Bootstrap hypothesis test (BootTest) method
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Identifies the number of common eigenvectors in several groups using the bootstrap hypothesis test (BootTest) method, as proposed by Klingenberg (1996).
}
\usage{
BootTest(origdata, q = ncol(origdata[[1]]), reps = 1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{origdata}{
%%     ~~Describe \code{origdata} here~~
List of the sample data sets.
}
  \item{q}{
%%     ~~Describe \code{q} here~~
Number of common eigenvectors to test for.
}
  \item{reps}{
%%     ~~Describe \code{reps} here~~
Number of bootstrap replications to use.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
Tests the hypothesis:
H_0: eigenvector pair are common
H_1: eigenvector pair are NOT common
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
Returns a data frame with the columns:
\item{Group1}{Order of the eigenvector for Group 1.}
\item{Group2}{Order of the eigenvector for Group 2.}
\item{vec.correlations}{Vector correlations of the eigenvector pairs.}
\item{p.values}{P values for the null hypothesis of commonness of the eigenvector pairs.}
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
Note that this implementation of the BootTest method can handle only two groups of data currently.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{ensemble.test}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (origdata, q = ncol(origdata[[1]]), reps = 1000) 
{
    k <- 2
    p <- ncol(origdata[[1]])
    nvec <- rep(NA, times = k)
    covmats <- array(NA, dim = c(p, p, k))
    for (i in 1:k) {
        nvec[i] <- nrow(origdata[[i]])
        covmats[, , i] <- cov(origdata[[i]])
    }
    B <- FG(covmats, nvec)$B
    findcpc.out <- findcpc(covmats, B = B, plotting = FALSE)
    commonvecnums <- findcpc.out$all.correlations[1:p, ]
    temp.p <- p
    for (i in 1:k) {
        j <- 2
        while (j <= temp.p) {
            if (length(unique(commonvecnums[1:j, i])) == length(unique(commonvecnums[1:(j - 
                1), i]))) {
                commonvecnums <- commonvecnums[-j, ]
                temp.p <- temp.p - 1
            }
            j <- j + 1
        }
    }
    commonvec.order <- commonvecnums[, "B"]
    eigenvec.order <- commonvecnums[, c("Group1", "Group2")]
    E <- array(NA, dim = c(p, p, k))
    data.rotated <- origdata
    for (i in 1:k) {
        E[, , i] <- eigen(covmats[, , i])$vectors[, eigenvec.order[, 
            i]]
        data.rotated[[i]] <- as.matrix(origdata[[i]]) \%*\% E[, 
            , i] \%*\% t(B)
    }
    orig.dotproducts <- abs(diag(t(E[, , 1]) \%*\% E[, , 2])[1:q])
    bootreps <- matrix(NA, ncol = q, nrow = reps)
    for (r in 1:reps) {
        rep.eigenvecs <- array(NA, dim = c(p, p, k))
        for (i in 1:k) {
            bootdata <- data.rotated[[i]][sample(c(1:nvec[i]), 
                size = nvec[i], replace = T), ]
            rep.eigenvecs[, , i] <- eigen(cov(bootdata))$vectors[, 
                eigenvec.order[, i]]
        }
        for (j in 1:q) {
            bootreps[r, j] <- abs(rep.eigenvecs[, j, 1] \%*\% rep.eigenvecs[, 
                j, 2])
        }
    }
    pvals <- rep(NA, times = q)
    for (j in 1:q) {
        pvals[j] <- nrow(bootreps[bootreps[, j] <= orig.dotproducts[j], 
            , drop = F])/reps
    }
    return(data.frame(eigenvec.order[1:q, c("Group1", "Group2"), 
        drop = F], vec.correlations = orig.dotproducts[1:q], 
        p.values = pvals[1:q]))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
