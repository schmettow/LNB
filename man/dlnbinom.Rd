\name{dlnbinom}
\alias{dlnbinom}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
dlnbinom(x, size, m, s)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{size}{
%%     ~~Describe \code{size} here~~
}
  \item{m}{
%%     ~~Describe \code{m} here~~
}
  \item{s}{
%%     ~~Describe \code{s} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
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
function (x, size, m, s) 
{
    x = as.integer(x)
    size = as.integer(size)
    f = function(P, x, size, m, s) (1 - P)^(size - x - 1) * P^(x - 
        1) * exp(-(logit(P) - m)^2/(2 * s^2))
    F = function(x, size, m, s) {
        if (s > 0) 
            return(integrate(f, 0, 1, x, size, m, s, rel.tol = .Machine$double.eps^0.25, 
                subdivisions = 100)$value * choose(size, x)/(sqrt(2 * 
                pi * s^2)))
        if (s == 0) 
            return(dbinom(x, size, exp(m)))
        if (s < 0) 
            return(0)
    }
    tryCatch(mapply(F, x, size, m, s), error = function(cond) {
        return(dlnbinom.mc(x, size, m, s, nruns = 10000))
    })
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
