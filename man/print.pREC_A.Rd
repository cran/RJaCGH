\name{print.pREC_A}
\alias{print.pREC_A}
\alias{print.pREC_A.none}
\alias{print.pREC_A.Chromosomes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Method for printing probabilistic common region. }
\description{
  A print method for \code{\link{pREC_A}} objects
}
\usage{
\method{print}{pREC_A}(x,...)
\method{print}{pREC_A.none}(x,...)
\method{print}{pREC_A.Chromosomes}(x,...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ An object of class \code{pREC_A},
    \code{pREC_A.Chromosomes}, \code{pREC_A.RJaCGH.none}.
    
  }
  \item{...}{ Additional arguments passed to \code{print}. Currently ignored.}
}
\value{
  A data.frame is printed with as many rows as regions found and with columns
  containing chromosome where the region is, position of start and end
  of the region, number of genes in it and joint probability.

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
  \code{\link{pREC_A}}
}
\examples{

## MCR for a single array:
y <- c(rnorm(100, 0, 1), rnorm(10, -3, 1), rnorm(20, 3, 1),
       rnorm(100,0, 1)) 
Pos <- sample(x=1:500, size=230, replace=TRUE)
Pos <- cumsum(Pos)
Chrom <- rep(1:23, rep(10, 23))

jp <- list(sigma.tau.mu=rep(0.05, 4), sigma.tau.sigma.2=rep(0.03, 4),
           sigma.tau.beta=rep(0.07, 4), tau.split.mu=0.1, tau.split.beta=0.1)

fit.genome <- RJaCGH(y=y, Pos=Pos, Chrom=Chrom, model="Genome",
burnin=100, TOT=1000, jump.parameters=jp, k.max = 4)
pREC_A(fit.genome, p=0.8, alteration="Gain")
pREC_A(fit.genome, p=0.8, alteration="Loss")

##MCR for two arrays:
z <- c(rnorm(110, 0, 1), rnorm(20, 3, 1),
       rnorm(100,0, 1)) 
fit.array.genome <- RJaCGH(y=cbind(y,z), Pos=Pos, Chrom=Chrom, model="Genome",
burnin=100, TOT=1000, jump.parameters=jp, k.max = 4)
pREC_A(fit.array.genome, p=0.4, alteration="Gain")
pREC_A(fit.array.genome, p=0.4, alteration="Loss")

}
\keyword{models}% at least one, from doc/KEYWORDS
