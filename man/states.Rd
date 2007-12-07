\name{states}
\alias{states}
\alias{states.RJaCGH}
\alias{states.RJaCGH.Chrom}
\alias{states.RJaCGH.genome}
\alias{states.RJaCGH.array}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{'states' method for RJaCGH objects}
\description{
Methods for estimating the hidden state sequence of a RJaCGH model.
}
\usage{
states(obj, k)
\method{states}{RJaCGH}(obj, k=NULL)
\method{states}{RJaCGH.Chrom}(obj, k=NULL)
\method{states}{RJaCGH.genome}(obj, k=NULL)
\method{states}{RJaCGH.array}(obj, k=NULL)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{any of RJaCGH, RJaCGH.Chrom, RJaCGH.genome, RJaCGH.array objects}
  \item{k}{Model to summarize (i.e., number of hidden states). If NULL,
    the most visited is taken.}
}
\details{
  The posterior probability of the hidden state sequence
  is computed via viterbi.



The state with more observatios is called 'Normal'. Those with bigger
means than it are called 'Gain', 'Gain1'... and those with lesser means
are called 'Loss', 'Loss1',...


Depending on the hierarchy of the object, it can return lists with
sublists, as in \code{\link{RJaCGH}}.
}
\value{
  \item{states}{Factor with the hidden state sequence}
  \item{prob.states}{Matrix with the probabilities associated to every
    states for every observation.}
}
\references{
    Rueda OM, Diaz-Uriarte R.
  Flexible and Accurate Detection of Genomic Copy-Number Changes from
  aCGH.
  PLoS Comput Biol. 2007;3(6):e122
}

  \author{Oscar M. Rueda and Ramon Diaz Uriarte}
  \seealso{\code{\link{RJaCGH}},
  \code{\link{summary.RJaCGH}}, \code{\link{model.averaging}},
  \code{\link{plot.RJaCGH}}, \code{\link{trace.plot}},
  \code{\link{gelman.brooks.plot}}, \code{\link{collapseChain}}}

\examples{
\dontrun{y <- c(rnorm(100, 0, 1), rnorm(10, -3, 1), rnorm(20, 3, 1), rnorm(100,
0, 1))
Pos <- sample(x=1:500, size=230, replace=TRUE)
Pos <- cumsum(Pos)
Chrom <- rep(1:23, rep(10, 23))
jp <- list(sigma.tau.mu=rep(0.05, 4), sigma.tau.sigma.2=rep(0.03, 4),
           sigma.tau.beta=rep(0.07, 4), tau.split.mu=0.1, tau.split.beta=0.1)
fit.genome <- RJaCGH(y=y, Pos=Pos, Chrom=Chrom, model="genome",
burnin=100, TOT=1000, jump.parameters=jp, k.max = 4)
states(fit.genome)}
}
\keyword{models}% at least one, from doc/KEYWORDS

