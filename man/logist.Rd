\name{logist}
\alias{logist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
logistic transformation
}
\description{
applies a logistic transformation
}
\usage{
logist(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
numeric
}
}
\details{
 The logistic transformation is applied, mapping x onto the interval [0;1] with a sigmoid curve. This is the reverse of the logit transformation and the same as the cumulative probability density of the logistic distribution.
}
\value{
numeric in the interval [0;1]
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Martin Schmettow
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{logit}}
\code{\link{plogis}}
}
\examples{
logist(rnorm(10,0,20))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
