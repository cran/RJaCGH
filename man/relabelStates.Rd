\name{relabelStates}
\alias{relabelStates}
\alias{relabelStates.RJaCGH}
\alias{relabelStates.RJaCGH.Chrom}
\alias{relabelStates.RJaCGH.Genome}
\alias{relabelStates.RJaCGH.array}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Relabelling of hidden states to biological states of alteration. }
\description{
  For every model, each hidden state is assigned to a state of copy number
  alteration ('normal', 'loss1', 'loss2', 'gain1', 'gain2'...)
}
\usage{
relabelStates(obj, normal.reference = 0,
              window = 1, singleState = FALSE)
\method{relabelStates}{RJaCGH}(obj, normal.reference = 0,
                               window = 1, singleState = FALSE)
\method{relabelStates}{RJaCGH.Chrom}(obj, normal.reference = 0,
                              window = 1, singleState = FALSE)
\method{relabelStates}{RJaCGH.Genome}(obj, normal.reference = 0,
                       window = 1, singleState = FALSE)
\method{relabelStates}{RJaCGH.array}(obj, normal.reference = 0,
                       window = 1, singleState = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{An object returned form RJaCGH of class 'RJaCGH',
  'RJaCGH.Chrom', 'RJaCGH.Genome', 'RJaCGH.array'.}
    \item{normal.reference}{The value considered as the mean of the normal
    state. See details. By default is \code{0}.}
  \item{window}{Multiplier of the standard deviation of the data to
    dtermine the width of the normal state. See details. Default is 1.}
    \item{singleState}{If TRUE, each state is assigned a probability of
      1 of being Gained or Lost or Normal, and 0 of being
      anyting else. This mimics the behavior of previous versions of
      RJaCGH. See details.}

}
\details{
  A relabelling of hidden states is performed to match biological
  states. We first define an upper and lower limit as: normal.reference
  +/- window. Next, we compute as the probability that a given state is
  gained the area beyond the upper limit, and the area up to the lower
  limit as the probability that a given state is lost. The areas are
  obtained from a Normal with posterior estimates of mean and
  standard deviation for each hidden state. The probability of being in
  a normal state is defined as 1 - prob.gained - prob.loss.
  
  Labels are assigned based on the state with highest probability. If we
  set singleState = TRUE, all the probability is assigned to a single
  condition, the one with largest probability.
}
\value{
An object of the same class as \code{obj} with hidden states relabelled.
}
\references{
  Rueda OM, Diaz-Uriarte R.
  Flexible and Accurate Detection of Genomic Copy-Number Changes from
  aCGH.
  PLoS Comput Biol. 2007;3(6):e122
}
\author{Oscar M. Rueda and Ramon Diaz Uriarte}
\seealso{\code{\link{RJaCGH}}}

\examples{
y <- c(rnorm(100, 0, 1), rnorm(10, -3, 1), rnorm(20, 3, 1),
       rnorm(100,0, 1)) 
Pos <- sample(x=1:500, size=230, replace=TRUE)
Pos <- cumsum(Pos)
Chrom <- rep(1:23, rep(10, 23))

jp <- list(sigma.tau.mu=rep(0.05, 4), sigma.tau.sigma.2=rep(0.03, 4),
           sigma.tau.beta=rep(0.07, 4), tau.split.mu=0.1, tau.split.beta=0.1)

fit.chrom <- RJaCGH(y=y, Pos=Pos, Chrom=Chrom, model="Chrom",
                    burnin=10, TOT=1000, k.max = 4,
                    jump.parameters=jp)
plot(fit.chrom)
fit.chrom.2 <- relabelStates(fit.chrom, normal.reference=3)
plot(fit.chrom.2)

}

\keyword{ models }% at least one, from doc/KEYWORDS

