\name{dlngeom}
\alias{dlngeom}
\alias{plngeom}
\alias{hlngeom}
\title{
logit-normal geometric distribution
}
\description{
distribution functions (probability, quantile and hazard) for the logit-normal geometric distribution
}
\usage{
dlngeom(x, size, m, s)
plngeom(x, size, m, s)
qlngeom(q, size, m, s)
hlngeom(x, size, m, s)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, q}{
vector of random variable realizations or quantiles
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
The logit-normal geometric is created by imposing a logit-normal prior on the parameter p of the geometric distribution.
}
\value{
probability (dlngeom), cumulative probabilty (qlngeom), quantile (qlngeom) or hazard rate (hlngeom)
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
dlnbinom
dlnbinom.zt
dgeom
}
\examples{
  dlngeom(3, 10, 0.5, 2)
  qlngeom(3, 10, 0.5, 2)
  }
  
  