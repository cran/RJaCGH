\name{smoothMeans}
\alias{smoothMeans}
\alias{smoothMeans.RJaCGH}
\alias{smoothMeans.RJaCGH.Chrom}
\alias{smoothMeans.RJaCGH.Genome}
\alias{smoothMeans.RJaCGH.array}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Smoothed posterior mean}
\description{
Smoothed posterior mean for every probe after fitting a RJaCGH model.
}
\usage{
smoothMeans(obj, k = NULL)
\method{smoothMeans}{RJaCGH}(obj, k=NULL)
\method{smoothMeans}{RJaCGH.Chrom}(obj, k=NULL)
\method{smoothMeans}{RJaCGH.Genome}(obj, k=NULL)
\method{smoothMeans}{RJaCGH.array}(obj, k=NULL)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{An \code{\link{RJaCGH}} object, of class 'RJaCGH',
    'RJaCGH.Chrom', 'RJaCGH.Genome' or 'RJaCGH.array'.
  }
  \item{k}{Number of states (or model) to get the smoothed means from.
  If NULL, Bayesian Model Averaging is used.}
}
\details{
  For a model with \code{k} hidden states, the mean from the MCMC samples
  from \code{mu} is computed for every hidden state.
  Then, for every probe these means are averaged by its posterior
  probability of belonging to every hidden state.
  If \code{k} is NULL, then this smoothed means are computed for every
  model and averaged by the posterior probability of each model.
  
}
\value{
  For class 'RJaCGH', 'RJaCGH.Chrom' and 'RJaCGH.Genome' a vector with
  the smoothed means for every probe. For class 'RJaCGH.array' a list
  with as many elements as arrays, each one a vector with the smoothed
  means for that array.
}
\references{
  Rueda OM, Diaz-Uriarte R.
  Flexible and Accurate Detection of Genomic Copy-Number Changes from
  aCGH.
  PLoS Comput Biol. 2007;3(6):e122
}
\author{Oscar M. Rueda and Ramon Diaz Uriarte}

\seealso{
  \code{\link{RJaCGH}}, \code{\link{plot.RJaCGH}}
  }
  \examples{
y <- c(rnorm(100, 0, 1), rnorm(10, -3, 1), rnorm(20, 3, 1),
       rnorm(100,0, 1)) 
Pos <- sample(x=1:500, size=230, replace=TRUE)
Pos <- cumsum(Pos)
Chrom <- rep(1:23, rep(10, 23))

jp <- list(sigma.tau.mu=rep(0.5, 4), sigma.tau.sigma.2=rep(0.3, 4),
           sigma.tau.beta=rep(0.7, 4), tau.split.mu=0.5)

fit.genome <- RJaCGH(y=y, Pos=Pos, Chrom=Chrom, model="Genome",
                    burnin=10, TOT=1000, k.max = 4,
                    jump.parameters=jp)
plot(y~Pos)
lines(smoothMeans(fit.genome) ~ Pos)
}
\keyword{models}% at least one, from doc/KEYWORDS
