\name{dlnbinom.zt}
\alias{dlnbinom.zot}
\title{
zero-truncated logit-normal binomial distribution
}
\description{
distribution functions (probability, random number generation and cumulative) for the tzero-ztruncated and zero-one trauncated logit-normal binomial distribution
}
\usage{
dlnbinom.zt(x, size, m, s)
dlnbinom.zot(x, size, m, s)
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
The logit-normal distribution is created by imposing a logit-normal prior on the binomial parameter p. The zero-truncated LNB is derived by restricting the random variable to a strictly positive range and re-scaling the probability mass accordingly. The zero-one truncated distribution is created respectively.
}
\value{
dlnbinom.zt gives zero-truncated probability, dlnbinom.zt returns the zero-one trauncated probability
}
\references{
Schmettow, M. (2009). Controlling the usability evaluation process under varying defect visibility. In BCS HCI 09: Proceedings of the 23rd British HCI Group Annual Conference on People and Computers: Celebrating People and Technology (pp. 188-197). Swinton, UK: British Computer Society.
}
\author{
Martin Schmettow
}
\note{
%%  ~~further notes~~
}


\seealso{
dlnbinom
dlngeom
}
\examples{
  dlnbinom.zt(3, 10, 0.5, 2)
  dlnbinom.zt(0, 10, 0.5, 2)
  }