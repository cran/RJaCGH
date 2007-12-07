\name{pREC_S}
\alias{pREC_S}
\alias{pREC_S.RJaCGH}
\alias{pREC_S.RJaCGH.array}
\alias{pREC_S.RJaCGH.array.Chrom}
\alias{pREC_S.RJaCGH.array.genome}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Subgroups of arrays that share common alterations }
\description{
  An algorithm to find regions of gain/lost copy number
  shared by a given proportion of arrays over a 
  probability threshold.
}
\usage{
pREC_S(obj, p, freq.array, alteration = "Gain")
\method{pREC_S}{RJaCGH.array}(obj, p, freq.array, alteration = "Gain")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{An object of class 'RJaCGH.array'. }
  \item{p}{ Threshold for the minimum joint probability
    of the region on every array.}
  \item{freq.array}{ Minimum number of arrays that share every
    region.}
  \item{alteration}{ Either 'Gain' or 'Loss'. }
}
\details{
  This algorithm, as \code{\link{pREC_A}}  computes
  probabilistic common regions
  but instead of finding regions that have
  a joint probability of alteration over all arrays, \code{pREC_S}
  searches for regions that have a probability of alteration
  higher than a threshold in at least a minimum number of
  arrays.
  So, \code{pREC_S} finds subsets of arrays that share subsets
  of alterations.
  Please note that if the method returns several sets or regions, the
  probability of alteration of all of them doesn't have to be over the
  probability threshold; in other words \code{p} is computed for every
  region, not for all the sequence of regions.
}
\value{
  An object of class \code{pREC_S.RJaCGH.array},
  \code{pREC_S.RJaCGH.array.Chrom} or
  \code{pREC_S.RJaCGH.array.genome}, as corresponding.
  They are lists with a sublist for every region encountered and
  elements:
  \item{start}{Start position of the region.}
  \item{indexStart}{index position of the start of the region.}
  \item{indexEnd}{index position of the end of the region.}
  \item{end}{End position of the region.}
  \item{members}{Arrays that share the region.}

  If there are chromosome information,
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
  This class supersedes the class \code{pMCR2}. Objects created with
  function \code{pMCR2} in older versions of RJaCGH can be converted into
  the new one in order to print them, or plot them.
  Two steps must be performed:
  1.-Change the class apropriately; if for example \code{obj} has
  class 'pMCR2.RJaCGH.array'  make
  \code{class(obj) <- 'pREC_S.RJaCGH.array'}
  2.-If \code{obj} has information about chromosomes, run
  \code{names(obj) <- 1:length(obj)}
}
\seealso{\code{\link{RJaCGH}},
  \code{\link{states}}, \code{\link{model.averaging}},
  \code{\link{print.pREC_S.RJaCGH.array}}
  \code{\link{plot.pREC_S.RJaCGH.array}}
  \code{\link{getSequence}}
  \code{\link{prob.seq}}
  \code{\link{pREC_A}}
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
}
\keyword{models}% at least one, from doc/KEYWORDS

