\name{print.pREC_S}
\alias{print.pREC_S}
\alias{print.pREC_S.none}
\alias{print.pREC_S.Chromosomes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Method for printing probabilistic common regions}
\description{
  A print method for \code{\link{pREC_S}} objects.
}
\usage{
\method{print}{pREC_S}(x,...)
\method{print}{pREC_S.Chromosomes}(x,...)
\method{print}{pREC_S.none}(x,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{An object of class \code{\link{pREC_S}},
    \code{pREC_S.Chromosomes} or
      \code{pREC_S.none}.}
  \item{\dots}{Additional arguments passed to \code{\link{print}}.
    Currently ignored.}
}
\value{
  A data.frame is printed with as many rows as regions found and with columns
  containing chromosome where the region is, position of start and end
  of the region and the arrays that belong to the region.

}
\references{
  Rueda OM, Diaz-Uriarte R.
  Flexible and Accurate Detection of Genomic Copy-Number Changes from
  aCGH.
  PLoS Comput Biol. 2007;3(6):e122

}

\author{Oscar M. Rueda and Ramon Diaz Uriarte}

\seealso{\code{\link{pREC_S}}}
\examples{
y <- c(rnorm(100, 0, 1), rnorm(10, -3, 1), rnorm(20, 3, 1),
       rnorm(100,0, 1)) 
Pos <- sample(x=1:500, size=230, replace=TRUE)
Pos <- cumsum(Pos)
Chrom <- rep(1:23, rep(10, 23))

jp <- list(sigma.tau.mu=rep(0.05, 4), sigma.tau.sigma.2=rep(0.03, 4),
           sigma.tau.beta=rep(0.07, 4), tau.split.mu=0.1)

z <- c(rnorm(110, 0, 1), rnorm(20, 3, 1),
       rnorm(100,0, 1)) 
zz <- c(rnorm(90, 0, 1), rnorm(40, 3, 1),
       rnorm(100,0, 1)) 

fit.array.genome <- RJaCGH(y=cbind(y,z,zz),
Pos=Pos, Chrom=Chrom, model="genome",
burnin=1000, TOT=1000, jump.parameters=jp, k.max = 4)
pREC_S(fit.array.genome, p=0.4, freq.array=2,
alteration="Gain")
pREC_S(fit.array.genome, p=0.4, freq.array=2, alteration="Loss")
}
\keyword{models}% at least one, from doc/KEYWORDS

