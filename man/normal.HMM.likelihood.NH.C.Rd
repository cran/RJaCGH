\name{normal.HMM.likelihood.NH.C}
\alias{normal.HMM.likelihood.NH.C}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Likelihood for non-homogeneous hidden Markov model}
\description{
This function returns the log-likelihood for RJaCGH model, a hidden
Markov model with normal distributed emissions and a non-homogeneous
transition matrix as computed by Q.NH.
}
\usage{
normal.HMM.likelihood.NH.C(y, x, mu, sigma.2, beta, stat = NULL,
q=-beta)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{Log Ratios observed}
  \item{x}{Vector of distances between genes}
  \item{mu}{Vector of means for the hidden states}
  \item{sigma.2}{Vector of variances for the hidden states}
  \item{beta}{\code{beta} in transition matrix}
  \item{stat}{Vector of initial probabilities. If NULL, a uniforma
    distribution is assumed.}
  \item{q}{\code{q} parameter in transition matrix. In RJaCGH model,
    \code{q}=\code{-beta}}
}
\details{
This function is just an interface for the C routine to compute
log-likelihood in RJaCGH model. 

}
\value{
It returns  a list with the same components passed plus:
  \item{loglik }{Log-likelihood}
}
\references{  Oscar M. Rueda and Ramon Diaz Uriarte. A flexible, accurate and
  extensible statistical method for detecting genomic copy-number
  changes. \url{http://biostats.bepress.com/cobra/ps/art9/}.
  {http://biostats.bepress.com/cobra/ps/art9/}.
}
\author{Oscar M. Rueda and Ramon Diaz-Uriarte}

\examples{

## create data
y <- c(rnorm(100, 0, 1), rnorm(50, 3, 1), rnorm(20, -3, 1),
rnorm(60, 0, 1))
x <- sample(1:1000, 229, replace=FALSE)
x <- x/max(x)
Chrom <- rep(1:23, rep(10, 23))
## same model for all genome
loglik <- 0
for (i in 1:23) {
loglik <- loglik + normal.HMM.likelihood.NH.C(y=y, x =x, mu=c(-3, 0, 3),
sigma.2=c(1,1,1), beta=matrix(c(0, 1, 1, 1, 0, 1, 1, 1, 0), 3))$loglik
}
loglik
}
\keyword{models}% at least one, from doc/KEYWORDS


