\name{pREC_A}
\alias{pREC_A}
% \alias{pREC_A.RJaCGH}
% \alias{pREC_A.RJaCGH.Chrom}
% \alias{pREC_A.RJaCGH.genome}
% \alias{pREC_A.RJaCGH.array}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Probabilistic Common Regions for copy number alteration. }
\description{
This method compute regions of gain/lost copy number with a joint
probability of alteration greater than a given threshold.

}
\usage{
pREC_A(obj, p, alteration = "Gain", array.weights = NULL,
       verbose = FALSE)
% \method{pREC_A}{RJaCGH}(obj, p, alteration = "Gain", array.weights = NULL)
% \method{pREC_A}{RJaCGH.Chrom}(obj, p, alteration = "Gain", array.weights = NULL)
% \method{pREC_A}{RJaCGH.genome}(obj, p, alteration = "Gain", array.weights
% = NULL)
% \method{pREC_A}{RJaCGH.array}(obj, p, alteration = "Gain", array.weights = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{ An object of class 'RJaCGH', 'RJaCGH.Chrom',
    'RJaCGH.Genome' or 'RJaCGH.array'.
  }
  \item{p}{ Threshold for the minimum joint probability of alteration of
    the region. }
  \item{alteration}{ Either 'Gain' or 'Lost' }
  \item{array.weights}{ When 'obj' contains several arrays, the user can
    give a weight to each of them according to their reliability or precision.}
%   \item{force.write.files}{If TRUE, the gzipped files with the Viterbi
%     sequence are written to disk even if they are, supposedly, already
%     there. }
  \item{verbose}{If TRUE provide more details on what is being done, including
    intermediate output from the C functions themselves. Only helpful
    for debugging or if you are bored; this will often write more output than
    you want.}
%   \item{delete.rewritten}{If TRUE, delete the files with the Viterbi
%     sequence from disk.}
}
\details{
  RJaCGH can compute common regions taking into account the
  probability of every probe to have an altered copy number. The result
  is a set of probes whose joint probability (not the product of their
  marginal probabilities, as returned by \code{\link{states}}
  or \code{\link{modelAveraging}})
  is at least as \emph{p} or greater. 

  Please note that if the method returns several sets or regions, the
  probability of alteration of all of them doesn't have to be over the
  probability threshold; in other words \code{p} is computed for every
  region, not for all the sequence of regions.
  
  Writing the files with the Vitterbi sequence to disk will be done by
  default if RJaCGH was run with the default "delete\_gzipped = TRUE"
  (which will happen if the global flag ".\_\_DELETE\_GZIPPED" is TRUE). To
  preserve disk space, as soon as the gzipped files are read into R (and
  incorporated into the RJaCGH object), by default they are deleted from
  disk. Sometimes, however, you might want not to delete them from disk;
  for instance, if you will continue working from this directory, and
  you want to save some CPU time. If the files exist in the directory
  when you call pREC there is no need to write them from R to disk,
  which allows you to save the time in the "writeBin" calls inside
  pREC. In this case, you would run RJaCGH with "delete\_gzipped =
  FALSE". Now, if for some reason those files are no longer available
  (you move directories, you delete them, etc), you should set
  "force.write.files = TRUE" (pREC will let you know if you need to do
  so).
  
  delete.rewritten helps prevent cluttering the disk. Files with the
  Viterbi sequence will be written to disk, read by C, and then
  deleted. Again, no information is lost, since the sequences are stored
  as part of the RJaCGH object. Note, however, that if you run RJaCGH
  with ".\_\_DELETE\_GZIPPED <- FALSE" this option has no effect, because
  it is implicit that you wanted, from the start, to preserve the
  files. In other words, "delete.rewritten" only has any effect if you
  either used "force.write.files = TRUE" or if you originally run RJaCGH
  without preserving the files in disk.
  
  
%   \code{pREC_A} calls the function \code{\link{getSequence}} that
%   creates temporal files in the working directory containing the
%   sequence of hidden states for every MCMC sample, so there should be
%   writing permisssions in that directory.  Then it calls repeatedly
%   \code{\link{prob.seq}} to compute the joint probability of sets of
%   probes over the MCMC samples.
}
\value{
  
  An object of class \code{pREC_A.none} or \code{pREC_A.Chromosomes},
  depending on whether or not the original RJaCGH had, or not, a
  Chromosomes component (the later only when the original object was of
  neither "RJaCGH.Chrom" or "RJaCGH.Genome").

  They are lists with a sublist for every region encountered and
  elements:

  \item{start}{Start position of the region.}
  \item{indexStart}{index position of the start of the region.}
  \item{indexEnd}{index position of the end of the region.}
  \item{end}{End position of the region.}
  \item{genes}{Number of genes in the region.}
  \item{prob}{Joint probability of gain/loss of the region.}

  If there are chromosome information (that is, the object inputed is
  of class \code{RJaCGH.Chrom}, \code{RJaCGH.Genome} or
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
  There have been major changes in how pREC is implemented. For details,
  see "Implementing\_pREC\_in\_C.pdf".
}

\seealso{\code{\link{RJaCGH}},
  \code{\link{states}}, \code{\link{modelAveraging}},
   \code{\link{print.pREC_A}}
%   \code{\link{getSequence}}
%   \code{\link{prob.seq}}
  \code{\link{pREC_S}}
}
\examples{
% \dontrun{
## MCR for a single array:
y <- c(rnorm(100, 0, 1), rnorm(10, -3, 1), rnorm(20, 3, 1),
       rnorm(100,0, 1)) 
Pos <- sample(x=1:500, size=230, replace=TRUE)
Pos <- cumsum(Pos)
Chrom <- rep(1:23, rep(10, 23))

jp <- list(sigma.tau.mu=rep(0.05, 4), sigma.tau.sigma.2=rep(0.03, 4),
           sigma.tau.beta=rep(0.07, 4), tau.split.mu=0.1, tau.split.beta=0.1)

fit.genome <- RJaCGH(y=y, Pos=Pos, Chrom=Chrom, model="Genome",
burnin=1000, TOT=1000, jump.parameters=jp, k.max = 4)
pREC_A(fit.genome, p=0.8, alteration="Gain")
pREC_A(fit.genome, p=0.8, alteration="Loss")

##MCR for two arrays:
z <- c(rnorm(110, 0, 1), rnorm(20, 3, 1),
       rnorm(100,0, 1)) 
fit.array.genome <- RJaCGH(y=cbind(y,z), Pos=Pos, Chrom=Chrom, model="Genome",
burnin=1000, TOT=1000, jump.parameters=jp, k.max = 4)
pREC_A(fit.array.genome, p=0.4, alteration="Gain")
pREC_A(fit.array.genome, p=0.4, alteration="Loss")
% }
}
\keyword{models}% at least one, from doc/KEYWORDS

