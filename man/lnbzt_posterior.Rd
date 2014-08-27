\name{lnbzt_posterior}
\alias{lnbzt_posterior}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
lnbzt_posterior(param, K, MS)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{param}{
%%     ~~Describe \code{param} here~~
}
  \item{K}{
%%     ~~Describe \code{K} here~~
}
  \item{MS}{
%%     ~~Describe \code{MS} here~~
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
function (param, K, MS) 
{
    if (param[2] >= 0) {
        occ = unique(MS)
        frq = freq(MS, occ)
        ml = log(dlnbinom.zt(occ, K, param[1], param[2]))
        loglik = sum(ml * frq)
    }
    else {
        loglik = NaN
    }
    mprior = dnorm(param[1], mean = -2, sd = 5, log = T)
    sprior = dunif(param[2], 0, 10, log = T)
    loglik + mprior + sprior
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
