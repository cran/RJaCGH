\name{getSequence}
\alias{getSequence}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Viterbi path from every MCMC sample}
\description{
  This function creates files with the
  maximum a posteriory (Viterbi) sequence of hidden states
  for every MCMC sample.
}
\usage{
getSequence(obj, filename, alteration)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{An object of class 'RJaCGH', 'RJaCGH.genome'
  or 'RJaCGH.Chrom'}
  \item{filename}{A string with a name for the file(s) to
    be created. See details.}
  \item{alteration}{Either 'Gain' or 'Loss'}
}
\details{
  'k' files are created (where 'k' is the number of different
  models fitted in 'obj') in the working directory.
  Each file contains as many rows as MCMC samples for that
  model. Each row describes the hidden state sequence, in the
  format 'state', '\t', 'breakpoint position', etc. until a
  '-1' showing the end of the sequence and a final number,
  the times that that particular path has occurred in the
  MCMC.
  
}
\value{
  'k' files are created with names 'filename' plus the number
  of the model (in fact, the number of hidden states in that
  model.
}
\references{
  Rueda OM, Diaz-Uriarte R.
  Flexible and Accurate Detection of Genomic Copy-Number Changes from
  aCGH.
  PLoS Comput Biol. 2007;3(6):e122
}
\author{Oscar M. Rueda and Ramon Diaz Uriarte}

\seealso{
  \code{\link{RJaCGH}}, \code{link{pMCR}}, \code{link{prob.seq}},
  \code{\link{pREC_A}}, \code{\link{pREC_S}}
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
}
}
\keyword{models}% at least one, from doc/KEYWORDS

