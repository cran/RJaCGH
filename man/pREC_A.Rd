\name{pREC_A}
\alias{pREC_A}
\alias{pREC_A.RJaCGH}
\alias{pREC_A.RJaCGH.Chrom}
\alias{pREC_A.RJaCGH.genome}
\alias{pREC_A.RJaCGH.array}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Probabilistic Common Regions for copy number alteration. }
\description{
This method compute regions of gain/lost copy number with a joint
probability of alteration greater than a given threshold.

}
\usage{
pREC_A(obj, p, alteration = "Gain", array.weights = NULL)
\method{pREC_A}{RJaCGH}(obj, p, alteration = "Gain", array.weights = NULL)
\method{pREC_A}{RJaCGH.Chrom}(obj, p, alteration = "Gain", array.weights = NULL)
\method{pREC_A}{RJaCGH.genome}(obj, p, alteration = "Gain", array.weights
= NULL)
\method{pREC_A}{RJaCGH.array}(obj, p, alteration = "Gain", array.weights = NULL)
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
  RJaCGH can compute common regions taking into account the
  probability of every probe to have an altered copy number. The result
  is a set of probes whose joint probability (not the product of their
  marginal probabilities, as returned by \code{\link{states}}
  or \code{\link{model.averaging}})
  is at least as \code{p} or greater. 

  Please note that if the method returns several sets or regions, the
  probability of alteration of all of them doesn't have to be over the
  probability threshold; in other words \code{p} is computed for every
  region, not for all the sequence of regions.

  \code{pREC_A} calls the function
  \code{\link{getSequence}} that creates temporal files
  in the working directory containing the sequence of hidden states
  for every MCMC sample, so there should be writing permisssions in
  that directory.
  Then it calls repeatedly \code{\link{prob.seq}} to compute the
  joint probability of sets of probes over the MCMC samples.
}
\value{
  An object of class \code{pREC_A.RJaCGH}, \code{pREC_A.RJaCGH.Chrom},
  \code{pREC_A.RJaCGH.genome}, \code{pREC_A.RJaCGH.array},
  \code{pREC_A.RJaCGH.array.Chrom} or
  \code{pREC_A.RJaCGH.array.genome}, as corresponding.
  They are lists with a sublist for every region encountered and
  elements:
  \item{start}{Start position of the region.}
  \item{indexStart}{index position of the start of the region.}
  \item{indexEnd}{index position of the end of the region.}
  \item{end}{End position of the region.}
  \item{genes}{Number of genes in the region.}
  \item{prob}{Joint probability of gain/loss of the region.}

  If there are chromosome information (that is, the object inputed is
  of class \code{RJaCGH.Chrom}, \code{RJaCGH.genome} or
  \code{RJaCGH.array} with each array of any of these classes), then
  this information will be enclosed in a list for each chromosome.
  
}

\references{
  Rueda OM, Diaz-Uriarte R.
  Flexible and Accurate Detection of Genomic Copy-Number Changes from
  aCGH.
  PLoS Comput Biol. 2007;3(6):e122
}

\author{Oscar M. Rueda and Ramon Diaz Uriarte}
  
\note{
  This is a preliminary method, so it can be slow.
  This class supersedes the class \code{pMCR}. Objects created with
  function \code{pMCR} in older versions of RJaCGH can be converted into
  the new one in order to print them. Two steps must be performed:
  1.-Change the class apropriately; if for example \code{obj} has
  class 'pMCR.RJaCGH.Chrom'  make
  \code{class(obj) <- 'pREC_A.RJaCGH.Chrom'}
  2.-If \code{obj} has information about chromosomes, run
  \code{names(obj) <- 1:length(obj)}
}

\seealso{\code{\link{RJaCGH}},
  \code{\link{states}}, \code{\link{model.averaging}},
  \code{\link{print.pREC_A}}
  \code{\link{getSequence}}
  \code{\link{prob.seq}}
  \code{\link{pREC_S}}
}
\examples{
\dontrun{
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
pREC_A(fit.genome, p=0.8, alteration="Gain")
pREC_A(fit.genome, p=0.8, alteration="Loss")

##MCR for two arrays:
z <- c(rnorm(110, 0, 1), rnorm(20, 3, 1),
       rnorm(100,0, 1)) 
fit.array.genome <- RJaCGH(y=cbind(y,z), Pos=Pos, Chrom=Chrom, model="genome",
burnin=1000, TOT=1000, jump.parameters=jp, k.max = 4)
pREC_A(fit.array.genome, p=0.4, alteration="Gain")
pREC_A(fit.array.genome, p=0.4, alteration="Loss")
}
}
\keyword{models}% at least one, from doc/KEYWORDS

