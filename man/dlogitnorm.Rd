\name{dlogitnorm}
\alias{dlogitnorm}
\alias{rlogitnorm}
\alias{mean.logitnorm}
\alias{var.logitnorm}

\title{
logit-normal distribution
}
\description{
distribution functions for the logit-normal distribution
}
\usage{
dlogitnorm(x, m, s)
rlogitnorm(r, x, m, s)
mean.logitnorm(m, s)
var.logitnorm(m, s)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  quantile
}
  \item{m}{
  mean
%%     ~~Describe \code{m} here~~
}
  \item{s}{
  standard deviation
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
function (x, m, s) 
dnorm(logit(x), m, s)/(x * (1 - x))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
