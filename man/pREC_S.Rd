\name{pREC_S}
\alias{pREC_S}
% \alias{pREC_S.RJaCGH}
% \alias{pREC_S.RJaCGH.array}
% \alias{pREC_S.RJaCGH.array.Chrom}
% \alias{pREC_S.RJaCGH.array.genome}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Subgroups of arrays that share common alterations }
\description{
  An algorithm to find regions of gain/lost copy number
  shared by a given proportion of arrays over a 
  probability threshold.
}
\usage{
pREC_S(obj, p, freq.array, alteration = "Gain",
       force.write.files = FALSE, verbose = FALSE,
       delete.rewritten = TRUE)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{An object of class 'RJaCGH.array'. }
  \item{p}{ Threshold for the minimum joint probability
    of the region on every array.}
  \item{freq.array}{ Minimum number of arrays that share every
    region.}
  \item{alteration}{ Either 'Gain' or 'Loss'. }
    \item{force.write.files}{If TRUE, the gzipped files with the Viterbi
    sequence are written to disk even if they are, supposedly, already
    there. }
  \item{verbose}{If TRUE provide more details on what is being done, including
    intermediate output from the C functions themselves. Only helpful
    for debugging or if you are bored; this will often write more output than
    you want.}
  \item{delete.rewritten}{If TRUE, delete the files with the Viterbi
    sequence from disk.}
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
  

}
\value{
  An object of class \code{pREC_S.none}, or \code{pREC_S.Chromosomes}
  depending on whether or not the original RJaCGH had, or not, a
  Chromosomes component (the later only when the original object was of
  neither "RJaCGH.Chrom" or "RJaCGH.Genome").
  
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
  There have been major changes in how pREC is implemented. For details,
  see "Implementing\_pREC\_in\_C.pdf".
}
\seealso{\code{\link{RJaCGH}},
  \code{\link{states}}, \code{\link{model.averaging}},
   \code{\link{print.pREC_S}}
  \code{\link{plot.pREC_S}}
  \code{\link{pREC_A}}
}
\examples{
% \dontrun{
## MCR for a single array:
y <- c(rnorm(100, 0, 1), rnorm(10, -3, 1), rnorm(20, 3, 1),
       rnorm(100,0, 1)) 
z <- c(rnorm(110, 0, 1), rnorm(20, 3, 1),
       rnorm(100,0, 1)) 
zz <- c(rnorm(90, 0, 1), rnorm(40, 3, 1),
       rnorm(100,0, 1)) 

Pos <- sample(x=1:500, size=230, replace=TRUE)
Pos <- cumsum(Pos)
Chrom <- rep(1:23, rep(10, 23))

jp <- list(sigma.tau.mu=rep(0.05, 4), sigma.tau.sigma.2=rep(0.03, 4),
           sigma.tau.beta=rep(0.07, 4), tau.split.mu=0.1)

fit.array.genome <- RJaCGH(y=cbind(y,z,zz),
Pos=Pos, Chrom=Chrom, model="Genome",
burnin=1000, TOT=1000, jump.parameters=jp, k.max = 4)
pREC_S(fit.array.genome, p=0.4, freq.array=2,
alteration="Gain")
pREC_S(fit.array.genome, p=0.4, freq.array=2, alteration="Loss")
% }
}
\keyword{models}% at least one, from doc/KEYWORDS

