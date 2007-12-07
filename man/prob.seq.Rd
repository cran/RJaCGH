\name{prob.seq}
\alias{prob.seq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Joint posterior probabilities of alteration}
\description{
  This function computes posterior joint probabilities of
  alteration for arbitrary (contiguous) sets of probes.
  They are computed counting the times the  particular
  sequence of alterations have ocurred in the Viterbi paths
  of every MCMC sample.
}
\usage{
prob.seq(obj, from, to, filename, alteration = "Gain")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{An object of class 'RJaCGH', 'RJaCGH.genome'
    or 'RJaCGH.Chrom'}
  \item{from}{Index of the starting probe to compute joint
    probability.}
  \item{to}{Index of the ending probe to compute joint
    probaility.}
  \item{filename}{name of the file with the Viterbi paths
    (created with \code{getSequence}.}
  \item{alteration}{either 'Gain' or 'Loss'. }
}
\details{
  Before using this function, a call to \code{\link{getSequence}}
  must be made in order to create the Viterbi paths.
  Note that the indexes for \code{from} and \code{to} are
  absolute, not relative to the chromosome.
}
\value{
  A joint probability of 'alteration' for the probes between
 'from' and 'to'.
}
\references{
    Rueda OM, Diaz-Uriarte R.
  Flexible and Accurate Detection of Genomic Copy-Number Changes from
  aCGH.
  PLoS Comput Biol. 2007;3(6):e122
}
\author{Oscar M. Rueda and Ramon Diaz Uriarte}

\seealso{
  \code{\link{RJaCGH}}, \code{\link{pREC_A}}, \code{link{\getSequence}},
  \code{\link{pREC_S}}
}
\examples{
\dontrun{
y <- c(rnorm(100, 0, 1), rnorm(10, -3, 1), rnorm(20, 3, 1),
       rnorm(100,0, 1)) 
Pos <- sample(x=1:500, size=230, replace=TRUE)
Pos <- cumsum(Pos)
Chrom <- rep(1:23, rep(10, 23))

jp <- list(sigma.tau.mu=rep(0.5, 4), sigma.tau.sigma.2=rep(0.3, 4),
           sigma.tau.beta=rep(0.7, 4), tau.split.mu=0.5, tau.split.beta=0.5)

fit.genome <- RJaCGH(y=y, Pos=Pos, Chrom=Chrom, model="genome",
                    burnin=10, TOT=1000, k.max = 4,
                    jump.parameters=jp)
getSequence(fit.genome, 'sequence', 'Gain')
prob.seq(fit.genome, 1, 4, 'sequence', 'Gain')
}
}
\keyword{models}% at least one, from doc/KEYWORDS

