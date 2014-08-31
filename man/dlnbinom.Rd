\name{dlnbinom}
\alias{dlnbinom}
\alias{rlnbinom}
\alias{plnbinom}
\title{
logit-normal binomial distribution
}
\description{
distribution functions (probability, random number generation and cumulative) for the logit-normal binomial distribution
}
\usage{
dlnbinom(x, size, m, s)
rlnbinom(n, size, m, s)
dlnbinom(p, size, m, s)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
vector of quantiles
}
  \item{size}{
number of trials (zero or more)
}
  \item{m}{
distribution parameter (central tendency)
}
  \item{s}{
distribution parameter (dispersion)
}
}
\details{
The logit-normal distribution is created by imposing a logit-normal prior on the binomial parameter p.
}
\value{
dlnbinom gives the density, rlnbinom produces random numbers and plnbinom gives the cumulative probability. 

The length of the result is determined by n for rbinom, and is the maximum of the lengths of the numerical parameters for the other functions. The numerical parameters other than n are recycled to the length of the result.
}
\references{
Schmettow, M. (2009). Controlling the usability evaluation process under varying defect visibility. In BCS-HCI ’09: Proceedings of the 23rd British HCI Group Annual Conference on People and Computers: Celebrating People and Technology (pp. 188–197). Swinton, UK: British Computer Society.
}
\author{
Martin Schmettow
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
dlnbinom.zt
dlngeom
}
\examples{
  dlnbinom(3, 10, 0.5, 2)
  rlnbinom(10, 10, 0.5, 2)
  plnbinom(0.5, 10, 0.5, 2)
  
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
