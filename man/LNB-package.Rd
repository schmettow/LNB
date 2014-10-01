\name{LNB-package}
\alias{LNB-package}
\alias{LNB}
\docType{package}
\title{Logit-normal binomial distribution
}
\description{
provides functions to work with the logit-normal binomial distribution and derivates
}
\details{
\tabular{ll}{
Package: \tab LNB\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-08-26\cr
License: \tab GPL-3\cr
}

}
\author{
Martin Schmettow

Maintainer: Martin Schmettow <schmettow@web.de>
}
\references{
Schmettow, M. (2009). Controlling the usability evaluation process under varying defect visibility. In BCS HCI 09: Proceedings of the 23rd British HCI Group Annual Conference on People and Computers: Celebrating People and Technology (pp. 188-197). Swinton, UK: British Computer Society.
}

\keyword{ package }
\seealso{

\code{\link[logitnorm:logitnorm-package]{logitnorm}}
}
\examples{
  dlnbinom(3, 10, 0.5, 2)
  dlngeom(3, 0.5, 2)
}
