\name{print.summary.RJaCGH}
\alias{print.summary.RJaCGH}
\alias{print.summary.RJaCGH.Chrom}
\alias{print.summary.RJaCGH.Genome}
\alias{print.summary.RJaCGH.array}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ print summary of RJaCGH fit}
\description{
The list stored in \code{\link{summary.RJaCGH}} is printed.
}
\usage{
\method{print}{summary.RJaCGH}(x, ...)
\method{print}{summary.RJaCGH.Chrom}(x, ...)
\method{print}{summary.RJaCGH.Genome}(x, ...)
\method{print}{summary.RJaCGH.array}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a \code{\link{summary.RJaCGH}} object}
  \item{\ldots}{Additional arguments passed to \code{\link{print}}}
}
\value{
  The summary is printed
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

