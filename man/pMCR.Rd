\name{pMCR}
\alias{pMCR}
\alias{pMCR.RJaCGH}
\alias{pMCR.RJaCGH.Chrom}
\alias{pMCR.RJaCGH.genome}
\alias{pMCR.RJaCGH.array}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Probabilistic Minimal Common Regions for copy number alteration. }
\description{
This method compute regions of gain/lost copy number with a joint
probability of alteration greater than a given threshold.

}
\usage{
pMCR(obj, p, alteration = "Gain", array.weights = NULL)
\method{pMCR}{RJaCGH}(obj, p, alteration = "Gain", array.weights = NULL)
\method{pMCR}{RJaCGH.Chrom}(obj, p, alteration = "Gain", array.weights = NULL)
\method{pMCR}{RJaCGH.genome}(obj, p, alteration = "Gain", array.weights
= NULL)
\method{pMCR}{RJaCGH.array}(obj, p, alteration = "Gain", array.weights = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{ An object of class 'RJaCGH', 'RJaCGH.Chrom',
    'RJaCGH.genome' or 'RJaCGH.array'.
  }
  \item{p}{ Threshold for the minimum joint probability of alteration of
    the region. }
  \item{alteration}{ Either 'Gain' or 'Lost' }
  \item{array.weights}{ When 'obj' contains several arrays, the user can
    give a weight to each of them according to their reliability or precision.}
}
\details{
  RJaCGH can compute minimal common regions taking into account the
  probability of every gene to have an altered copy number. The result
  is a set of genes whose joint probability (not the product of their
  marginal probabilities, as returned by 'states' or 'model.averaging'
  is at least as 'p' or greater. 

  Please note that if the method returns several sets or regions, the
  probability of alteration of all of them doesn't have to be over the
  probability threshold; in other words 'p' is computed for every
  region, not for all the sequence of regions.
  
}
\value{
  An object of class \texttt{pMCR.RJaCGH}, \texttt{pMCR.RJaCGH.Chrom}
  or \texttt{RJaCGH.genome}, as corresponding.
  They are lists with a sublist for every region encountered and
  elements:
  \item{start}{Start position of the region.}
  \item{indexStart}{index position of the start of the region.}
  \item{indexEnd}{index position of the end of the region.}
  \item{end}{End position of the region.}
  \item{genes}{Number of genes in the region.}
  \item{prob}{Joint probability of gain/loss of the region.}

  If there are chromosome information (that is, the object inputed is
  of class \texttt{RJaCGH.Chrom}, \texttt{RJaCGH.genome} or
  \texttt{RJaCGH.array} with each array of any of these classes), then
  those information will be enclosed in a list for each chromosome.
  
  
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
  \code{\link{print.pMCR.RJaCGH}}
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
burnin=1000, TOT=1000, jump.parameters=jp, k.max = 4)
pMCR(fit.genome, p=0.8, alteration="Gain")
pMCR(fit.genome, p=0.8, alteration="Loss")

##MCR for two arrays:
z <- c(rnorm(110, 0, 1), rnorm(20, 3, 1),
       rnorm(100,0, 1)) 
fit.array.genome <- RJaCGH(y=cbind(y,z), Pos=Pos, Chrom=Chrom, model="genome",
burnin=1000, TOT=1000, jump.parameters=jp, k.max = 4)
pMCR(fit.array.genome, p=0.4, alteration="Gain")
pMCR(fit.array.genome, p=0.4, alteration="Loss")

}
\keyword{models}% at least one, from doc/KEYWORDS

