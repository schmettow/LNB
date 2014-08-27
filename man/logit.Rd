\name{logit}
\alias{logit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{logit}
\description{
Taking the logarithm of an odds
}
\usage{
logit(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{numeric value between 0 and 1}
}
\details{The logit function log(p/(1-p)) is the inverse of the logistic function. It takes a value between 0 and 1 (e.g., a probability) as input and returns a real value. It is often use to map bounded values to real numbers.}
\value{numeric}

\author{
Martin Schmettow
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{logist}}
\code{\link{logitnorm}}
\code{\link{lnbinom}}
}
\examples{
logit(runif(10))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
