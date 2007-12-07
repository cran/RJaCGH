\name{relabelStates}
\alias{relabelStates}
\alias{relabelStates.RJaCGH}
\alias{relabelStates.RJaCGH.Chrom}
\alias{relabelStates.RJaCGH.genome}
\alias{relabelStates.RJaCGH.array}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Relabelling of hidden states to biological states of alteration. }
\description{
  For every model, each hidden state is assigned to a state of copy number
  alteration ('normal', 'loss1', 'loss2', 'gain1', 'gain2'...)
}
\usage{
relabelStates(obj, normal.reference = 0, normal.ref.percentile = 0.95, auto.label = NULL)
\method{relabelStates}{RJaCGH}(obj, normal.reference = 0,
normal.ref.percentile = 0.95, auto.label = NULL)
\method{relabelStates}{RJaCGH.Chrom}(obj, normal.reference = 0,
normal.ref.percentile = 0.95, auto.label = NULL)
\method{relabelStates}{RJaCGH.genome}(obj, normal.reference = 0,
normal.ref.percentile = 0.95, auto.label = NULL)
\method{relabelStates}{RJaCGH.array}(obj, normal.reference = 0,
normal.ref.percentile = 0.95, auto.label = NULL)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{An object returned form RJaCGH of class 'RJaCGH',
  'RJaCGH.Chrom', 'RJaCGH.genome', 'RJaCGH.array'.}
    \item{normal.reference}{The value considered as the mean of the normal
    state. See details. By default is \code{0}.}
  \item{normal.ref.percentile}{Percentage for the relabelling of
    states. See details. by default is 0.95.}
    \item{auto.label}{If not NULL, should be the minimum proportion of
  observations labeled as 'Normal'. See details.}

}
\details{
  A relabelling of hidden states is performed to match biological
  states. The states that have the \code{normal.reference} value
  inside a \code{normal.ref.percentile}\% probability interval
  based on a normal distribution with means the median of \code{mu}
  and sd the square root of the median of \code{sigma.2} are labelled as
  'Normal'. If no state is close enough to \code{normal.reference} then
  there will not be a normal state. Bear this in mind for
  normalization issues.
  If \code{auto.label} is not null, closest states to 'Normal' are also
  labelled as 'Normal' until a proportion of \code{auto.label} is
  reached. Please note that the default value is 0.60, so at least the
  60\% of the observations will be labelled as 'Normal'.
  If this laeblling is not satisfactory, you can relabel manually. See 
  the example.
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

## Manual labelling
fit.chrom.2[[1]][[2]]$state.labels <- c("Normal", "Normal")
}

\keyword{ models }% at least one, from doc/KEYWORDS

