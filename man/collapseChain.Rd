\name{collapseChain}
\alias{collapseChain}
\alias{collapseChain.RJaCGH}
\alias{collapseChain.RJaCGH.Chrom}
\alias{collapseChain.RJaCGH.Genome}
\alias{collapseChain.RJaCGH.array}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Collapse several parallel chains ('RJaCGH' objects)}
\description{
method to join or collapse several 'RJaCGH' objects, for use in every
method of class 'RJaCGH'.
}
\usage{
collapseChain(obj)
\method{collapseChain}{RJaCGH}(obj)
\method{collapseChain}{RJaCGH.Genome}(obj)
\method{collapseChain}{RJaCGH.Chrom}(obj)
\method{collapseChain}{RJaCGH.array}(obj)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{a list containing several parallel chains; that is objects
    of any of RJaCGH, RJaCGH.Chrom, RJaCGH.Genome, RJaCGH.array classes
    (obviously, all of the same class).}
}
\details{
If several parallel chains are run and if they converge (see
\code{\link{gelman.rubin.plot}}) they should be joined in one. This is
what this method does.
Please note that this function returns only one object, so one can not
call later \code{gelman.rubin.plot}.
}
\value{
An object of the same class as any of the list \code{obj}.
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
fit.genome[[i]] <- RJaCGH(y=y, Pos=Pos, Chrom=Chrom, model="Genome",
burnin=10, TOT=1000, jump.parameters=jp, k.max = 4)
}

\dontrun{gelman.rubin.plot(fit.genome)}
##If all R seem to be round 1
\dontrun{fit.genome <- collapseChain(fit.genome)}
}
\keyword{models}% at least one, from doc/KEYWORDS
