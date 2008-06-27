\name{trace.plot}
\alias{trace.plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Trace plot for 'RJaCGH' object}
\description{
A trace plot with the trajectory of the Markov Chain.
}
\usage{
trace.plot(obj, k = NULL, array = NULL, Chrom = NULL, main.text = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{any of RJaCGH, RJaCGH.Chrom, RJaCGH.Genome, RJaCGH.array objects}
  \item{k}{Model to plot (i.e., number of hidden states). If NULL,
    the most visited is taken.}
  \item{array}{if \code{obj} is 'RJaCGH.array', the name of the array to
  plot must be given.}
  \item{Chrom}{if \code{obj} is 'RJaCGH.Chrom', the number of the chromosome to
  plot must be given.}
  \item{main.text}{Main text of the plot}
}
\details{
This is simply a call to \code{\link{matplot}} to show the values
sampled in the chain. \newline
The colors does not correspond to any particular level of gain/loss.
}
\value{
  A plot is drawn.
}
\references{
    Rueda OM, Diaz-Uriarte R.
  Flexible and Accurate Detection of Genomic Copy-Number Changes from
  aCGH.
  PLoS Comput Biol. 2007;3(6):e122
}
\author{Oscar M. Rueda and Ramon Diaz}

\seealso{\code{\link{RJaCGH}},
  \code{\link{summary.RJaCGH}}, \code{\link{model.averaging}},
  \code{\link{plot.RJaCGH}}, \code{\link{states}},
  \code{\link{gelman.rubin.plot}}, \code{\link{collapseChain}}}

\examples{
y <- c(rnorm(100, 0, 1), rnorm(10, -3, 1), rnorm(20, 3, 1), rnorm(100,
0, 1))
Pos <- sample(x=1:500, size=230, replace=TRUE)
Pos <- cumsum(Pos)
Chrom <- rep(1:23, rep(10, 23))
jp <- list(sigma.tau.mu=rep(0.5, 5), sigma.tau.sigma.2=rep(0.3, 5),
sigma.tau.beta=rep(0.7, 5), tau.split.mu=0.5, tau.split.beta=0.5)
fit.genome <- RJaCGH(y=y, Pos=Pos, Chrom=Chrom, model="Genome",
burnin=10, TOT=100, jump.parameters=jp, k.max = 5)
trace.plot(fit.genome)
}
\keyword{models}% at least one, from doc/KEYWORDS

