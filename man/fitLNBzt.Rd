\name{fitLNBzt}
\alias{fitLNBzt}
\alias{Dhat}
\alias{Dnull}
\title{
ML estimation of zero-truncated LNB
}
\description{
takes a margin sum data set and the size (number of sessions) and fits by ML
for extreme sparse or strongly dispersed data, the fit may fail.
Experiment with the start values (startval) in that case.
}
\usage{
fitLNBzt(ms,n,startval=c(-1,2))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ms}{
vector with margin sum
}
  \item{n}{
number of trials
}
  \item{K}{
starting values for mu and sd
}
}
\details{
takes a margin sum data set and the size (number of sessions) and fits by ML
for extreme sparse or strongly dispersed data, the fit may fail.
Experiment with the start values (startval) in that case.
}
\value{
fitLNB
\item{mu}{}
\item{sd}{}
\item{n}{}
\item{discovered}{}
\item{item.ms}{}
\item{nlogLik}{}
\item{AIC}{}

}
\references{
Schmettow, M. (2009). Controlling the usability evaluation process under varying defect visibility. In BCS HCI 09: Proceedings of the 23rd British HCI Group Annual Conference on People and Computers: Celebrating People and Technology (pp. 188-197). Swinton, UK: British Computer Society.
}
\author{
Martin Schmettow
}
\note{
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{dlnbinom.zt}}
}
\examples{
  ms <- rlnbinom(100, 15, 0.5, 1.2)
  ms <- ms[ms > 0]
  fit <- fitLNBzt(ms, 15)
  fit
  D.hat(fit)
  D.null
}
