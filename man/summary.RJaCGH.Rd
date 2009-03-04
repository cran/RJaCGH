\name{summary.RJaCGH}
\alias{summary.RJaCGH}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Summarizing RJaCGH models}
\description{
'summary' method for objects of class 'RJaCGH'.
}
\usage{
\method{summary}{RJaCGH}(object, array=NULL, Chrom=NULL, 
k = NULL, point.estimator = "median",
quantiles=NULL, \ldots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{RJaCGH objects}
  \item{array}{vector of names of arrays to summarize. If NULL, all of them.}
  \item{Chrom}{vector of chromosomes to summarize. If NULL, all of them.}
  \item{k}{Model to summarize (i.e., number of hidden states). If NULL,
    the most visited is taken.}
  \item{point.estimator}{Type of point estimator for \emph{mu},
    \emph{sigma.2} and \emph{beta}. It can be "mean", "median" or
    "mode".}
  \item{quantiles}{A vector of probabilities for the quantiles of the
    posterior distribution of means and variances.}
  \item{...}{Additional arguments passed to summary.}
}
\details{
Depending of the arguments passed, a list with contains sublists can be
returned, similarly to RJaCGH and similar objects of the family.
The point estimator "mode" is simply the max value obtained in a kernel
density estimation through the function \code{\link{density}}
}
\value{
  \item{k}{Frequencies of the hidden states visited by the sampler.}
  \item{mu}{Quantiles of the posterior distribution of \emph{mu}}
  \item{sigma.2}{Quantiles of the posterior distribution of \emph{sigma.2}}
  \item{beta}{Point estimator of \emph{beta}}
  \item{stat}{Initial distribution of the hidden states.}
}
\references{
  Rueda OM, Diaz-Uriarte R.
  Flexible and Accurate Detection of Genomic Copy-Number Changes from
  aCGH.
  PLoS Comput Biol. 2007;3(6):e122
}
\author{Oscar M. Rueda and Ramon Diaz Uriarte}

\seealso{\code{\link{RJaCGH}},
  \code{\link{states}}, \code{\link{modelAveraging}},
  \code{\link{plot.RJaCGH}}, \code{\link{trace.plot}},
%   \code{\link{gelman.rubin.plot}}, \code{\link{collapseChain}}
}

\examples{
y <- c(rnorm(100, 0, 1), rnorm(10, -3, 1), rnorm(20, 3, 1), rnorm(100,
0, 1))
Pos <- sample(x=1:500, size=230, replace=TRUE)
Pos <- cumsum(Pos)
Chrom <- rep(1:23, rep(10, 23))
jp <- list(sigma.tau.mu=rep(0.5, 5), sigma.tau.sigma.2=rep(0.3, 5),
sigma.tau.beta=rep(0.7, 5), tau.split.mu=0.5, tau.split.beta=0.5)
fit.chrom <- RJaCGH(y=y, Pos=Pos, Chrom=Chrom, model="Chrom",
burnin=10, TOT=100, jump.parameters=jp, k.max = 5)
summary(fit.chrom)
}
\keyword{models}% at least one, from doc/KEYWORDS

