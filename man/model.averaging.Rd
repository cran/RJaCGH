\name{model.averaging}
\alias{model.averaging}
\alias{model.averaging.RJaCGH}
\alias{model.averaging.RJaCGH.Chrom}
\alias{model.averaging.RJaCGH.Genome}
\alias{model.averaging.RJaCGH.array}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Method for model averaging for RJaCGH objects.}
\description{
Bayesian model averaging for the estimation of hidden state sequence.
}
\usage{
model.averaging(obj)
\method{model.averaging}{RJaCGH}(obj)
\method{model.averaging}{RJaCGH.Chrom}(obj)
\method{model.averaging}{RJaCGH.Genome}(obj)
\method{model.averaging}{RJaCGH.array}(obj)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{An object of corresponding class}
}
\details{
  With the posterior distribution of the number of hidden states,
  bayesian model averaging is performed on every model using
  \code{\link{states}} method. \newline
  As the other methods, it may return a list with sublists according to
  the hierarchy of RJaCGH objects.
}
\value{
  \item{states}{Factor with the hidden state sequence}
  \item{prob.states}{Matrix with the probabilities associated to every
    states for every observation.}
%   \item{prob.joint.states}{Matrix with the joint probabilities
%     associated
%     to each pair of states for every observation.}

}
\references{
    Rueda OM, Diaz-Uriarte R.
  Flexible and Accurate Detection of Genomic Copy-Number Changes from
  aCGH.
  PLoS Comput Biol. 2007;3(6):e122
}
\author{Oscar M. Rueda and Ramon Diaz Uriarte}

\seealso{\code{\link{RJaCGH}},
  \code{\link{summary.RJaCGH}}, \code{\link{states}},
  \code{\link{plot.RJaCGH}}, \code{\link{trace.plot}},
  \code{\link{gelman.rubin.plot}}, \code{\link{collapseChain}}}

\examples{
\dontrun{y <- c(rnorm(100, 0, 1), rnorm(10, -3, 1), rnorm(20, 3, 1),
       rnorm(100,0, 1)) 
Pos <- sample(x=1:500, size=230, replace=TRUE)
Pos <- cumsum(Pos)
Chrom <- rep(1:23, rep(10, 23))

jp <- list(sigma.tau.mu=rep(0.5, 5), sigma.tau.sigma.2=rep(0.3, 5),
sigma.tau.beta=rep(0.7, 5), tau.split.mu=0.5, tau.split.beta=0.5)
fit.genome <- RJaCGH(y=y, Pos=Pos, Chrom=Chrom, model="Genome",
burnin=1000, TOT=10000, jump.parameters=jp, max.k=5)
mo <- model.averaging(fit.genome)
print(mo)}

}

\keyword{models}% at least one, from doc/KEYWORDS

