\name{print.pMCR.RJaCGH}
\alias{print.pMCR.RJaCGH}
\alias{print.pMCR.RJaCGH.Chrom}
\alias{print.pMCR.RJaCGH.genome}
\alias{print.pMCR.RJaCGH.array}
\alias{print.pMCR.RJaCGH.array.Chrom}
\alias{print.pMCR.RJaCGH.array.genome}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Method for printing probabilistic minimal common region. }
\description{
  A print method for \code{\link{pMCR}} objects
}
\usage{
\method{print}{pMCR.RJaCGH}(x,...)
\method{print}{pMCR.RJaCGH.Chrom}(x,...)
\method{print}{pMCR.RJaCGH.genome}(x,...)
\method{print}{pMCR.RJaCGH.array}(x,...)
\method{print}{pMCR.RJaCGH.array.Chrom}(x,...)
\method{print}{pMCR.RJaCGH.array.genome}(x,...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ An object of class \code{pMCR.RJaCGH},
    \code{pMCR.RJaCGH.Chrom}, \code{pMCR.RJaCGH.genome},
    \code{pMCR.RJaCGH.array}, \code{pMCR.RJaCGH.array.Chrom} or
    \code{pMCR.RJaCGH.array.genome}.
    
  }
  \item{...}{ Additional arguments passed to 'print'. Currently ignored.}
}
\details{

}
\value{
  A matrix is printed with as many rows as regions found and with columns
  containing chromosome where the region is, position of start and end
  of the region, number of genes in it and joint probability.

}
\references{
  Oscar M. Rueda and Ramon Diaz Uriarte. A flexible, accurate and
  extensible statistical method for detecting genomic copy-number
  changes. \url{http://biostats.bepress.com/cobra/ps/art9/}
  {http://biostats.bepress.com/cobra/ps/art9/}.
}

\author{Oscar M. Rueda and Ramon Diaz Uriarte}
  
\seealso{\code{\link{RJaCGH}},
  \code{\link{states}}, \code{\link{model.averaging}},
  \code{\link{pMCR}}
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

fit.genome <- RJaCGH(y=y, Pos=Pos, Chrom=Chrom, model="genome",
burnin=100, TOT=1000, jump.parameters=jp, k.max = 4)
pMCR(fit.genome, p=0.8, alteration="Gain")
pMCR(fit.genome, p=0.8, alteration="Loss")

##MCR for two arrays:
z <- c(rnorm(110, 0, 1), rnorm(20, 3, 1),
       rnorm(100,0, 1)) 
fit.array.genome <- RJaCGH(y=cbind(y,z), Pos=Pos, Chrom=Chrom, model="genome",
burnin=100, TOT=1000, jump.parameters=jp, k.max = 4)
pMCR(fit.array.genome, p=0.4, alteration="Gain")
pMCR(fit.array.genome, p=0.4, alteration="Loss")

}
\keyword{models}% at least one, from doc/KEYWORDS
