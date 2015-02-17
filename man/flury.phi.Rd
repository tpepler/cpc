\name{flury.phi}
\alias{flury.phi}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Flury's phi measure
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Calculates phi measure of goodness of diagonalisation as proposed by Flury (1988).
}
\usage{
flury.phi(datamat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{datamat}{
%%     ~~Describe \code{datamat} here~~
A square symmetric matrix.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
A smaller phi value points to better diagonalisation, and a larger value to worse diagonalisation.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
Returns the value of the phi measure for the matrix.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
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
function (datamat) 
{
    return(det(diag(diag(datamat)))/det(datamat))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
