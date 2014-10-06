\name{fitLNB}
\alias{fitLNB}
\alias{fitLNBzt}
\alias{dhat}
\alias{dnull}
\alias{dtotal}
\alias{print}
\alias{summary}
\alias{coef}
\alias{plot}
\title{
ML estimation of (zero-truncated) LNB
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
\item{mu}{mean of logitnormal prior}
\item{sd}{standard deviation of logitnormal prior}
\item{n}{number of sessions/trials}
\item{discovered}{number of responses > 0}
\item{ms}{margin sum}
\item{nlogLik}{negative log-likelihood}
\item{AIC}{Akaike Information Criterion}
\item{zt}{zero-truncated estimation (logical, TRUE with fitLNBzt)}
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
  summary(fit)
  dhat(fit)
  dnull(fit)
  plot(fit)
}
