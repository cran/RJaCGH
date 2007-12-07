\name{chainsSelect}
\alias{chainsSelect}
\alias{chainsSelect.RJaCGH}
\alias{chainsSelect.RJaCGH.genome}
\alias{chainsSelect.RJaCGH.array}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Select between several parallel chains ('RJaCGH' objects)}
\description{
  method to remove outlier chains ('RJaCGH' objects') from a set of
  parallel RJMCMC chains.
}
\usage{
chainsSelect(obj, nutrim=4, trim=NULL)
\method{chainsSelect}{RJaCGH}(obj, nutrim=4, trim=NULL)
\method{chainsSelect}{RJaCGH.genome}(obj, nutrim = 4, trim = NULL)
\method{chainsSelect}{RJaCGH.array}(obj, nutrim = 4, trim = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{a list containing several parallel chains; that is objects
    of any of RJaCGH, RJaCGH.genome, RJaCGH.array classes
    (obviously, all of the same class).}
  \item{nutrim}{Number of chains to remove.}
  \item{trim}{Proportion of chains to remove. This or \code{nutrim}
  must be passed to the method.}
}
\details{
  In Reversible Jump MCMC, there may be occasions in that a chain
  can get trapped in a model. If this happens, Gelman-brooks - 
  see \code{\link{gelman.brooks.plot}} diagnostics
  for \code{k}, the number of states is not defined. This method removes
  the chains whose k trajectories are most unusual compared to the rest.
  In the current version of RJaCGH, this method is not defined for
  chromosome models (\code{RJaCGH.Chrom}).
  After removing otilers chains, one can run \code{gelman.brooks.plot}
  to check for convergnce of the remaining chains and later
  \code{collapseChain} to join them, or directly this last function.
}

\value{
A list of less elements of the same class as\code{obj}.
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
  \code{\link{plot.RJaCGH}}, \code{\link{states}},
  \code{\link{trace.plot}}, \code{\link{collapseChain}}}

\examples{
y <- c(rnorm(100, 0, 1), rnorm(10, -3, 1), rnorm(20, 3, 1),
       rnorm(100,0, 1)) 
Pos <- sample(x=1:500, size=230, replace=TRUE)
Pos <- cumsum(Pos)
Chrom <- rep(1:23, rep(10, 23))

jp <- list(sigma.tau.mu=rep(0.5, 4), sigma.tau.sigma.2=rep(0.3, 4),
           sigma.tau.beta=rep(0.7, 4), tau.split.mu=0.5, tau.split.beta=0.5)

fit.chrom <- RJaCGH(y=y, Pos=Pos, Chrom=Chrom, model="Chrom",
                    burnin=10, TOT=1000, k.max = 4,
                    jump.parameters=jp)
fit.genome <- list()
for (i in 1:4) {
fit.genome[[i]] <- RJaCGH(y=y, Pos=Pos, Chrom=Chrom, model="genome",
burnin=10, TOT=1000, jump.parameters=jp, k.max = 4)
fit.genome <- chainsSelect(fit.genome, nutrim=1)
}
}
\keyword{models}% at least one, from doc/KEYWORDS
