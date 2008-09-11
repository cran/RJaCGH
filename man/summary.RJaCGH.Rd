\name{summary.RJaCGH}
\alias{summary.RJaCGH}
\alias{summary.RJaCGH.Chrom}
\alias{summary.RJaCGH.Genome}
\alias{summary.RJaCGH.array}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Summarizing RJaCGH models}
\description{
'summary' method for objects of class 'RJaCGH', 'RJaCGH.Chrom',
'RJaCGH.Genome' and 'RJaCGH.array'.
}
\usage{
\method{summary}{RJaCGH}(object, k = NULL, point.estimator = "median",
quantiles=NULL, \ldots)
\method{summary}{RJaCGH.Chrom}(object, point.estimator="median",
Chrom=NULL,quantiles=NULL,  \ldots)
\method{summary}{RJaCGH.Genome}(object, k=NULL,
point.estimator="median",
quantiles=NULL, \ldots)
\method{summary}{RJaCGH.array}(object, point.estimator="median",
quantiles=NULL, \ldots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{any of RJaCGH, RJaCGH.Chrom, RJaCGH.Genome, RJaCGH.array objects}
  \item{k}{Model to summarize (i.e., number of hidden states). If NULL,
    the most visited is taken.}
  \item{point.estimator}{Type of point estimator for \code{mu},
    \code{sigma.2} and \code{beta}. It can be "mean", "median" or
    "mode".}
  \item{Chrom}{If \code{object} has a different model for every
    chromosome, this argument is the chromosome to get summaries from.
    If NULL, all chromosomes will be summarized.}
  \item{quantiles}{A vector of probabilities for the quantiles of the
    posterior distribution of means and variances.}
  \item{...}{Additional arguments passed to summary.}
}
\details{
Depending of the type of object, a list with contains sublists can be
returned, similarly to RJaCGH and similar objects of the family.
The point estimator "mode" is simply the max value obtained in a kernel
density estimation through the function \code{\link{density}}
}
\value{
  \item{k}{Frequencies of the hidden states visited by the sampler.}
  \item{mu}{Quantiles of the posterior distribution of \code{mu}}
  \item{sigma.2}{Quantiles of the posterior distribution of \code{sigma.2}}
  \item{beta}{Point estimator of \code{beta}}
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
  \code{\link{states}}, \code{\link{model.averaging}},
  \code{\link{plot.RJaCGH}}, \code{\link{trace.plot}},
  \code{\link{gelman.rubin.plot}}, \code{\link{collapseChain}}}

\examples{
y <- c(rnorm(100, 0, 1), rnorm(10, -3, 1), rnorm(20, 3, 1), rnorm(100,
0, 1))
Pos <- sample(x=1:500, size=230, replace=TRUE)
Pos <- cumsum(Pos)
Chrom <- rep(1:23, rep(10, 23))
jp <- list(sigma.tau.mu=rep(0.5, 5), sigma.tau.sigma.2=rep(0.3, 5),
sigma.tau.beta=rep(0.7, 5), tau.split.mu=0.5)
fit.chrom <- RJaCGH(y=y, Pos=Pos, Chrom=Chrom, model="Chrom",
burnin=10, TOT=100, jump.parameters=jp, k.max = 5)
summary(fit.chrom)
}
\keyword{models}% at least one, from doc/KEYWORDS

