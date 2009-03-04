\name{plot.RJaCGH}
\alias{plot.RJaCGH}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{'plot' method for RJaCGH objects}
\description{
A plot is drawed with information from the fit of a RJaCGH object.
}
\usage{
\method{plot}{RJaCGH}(x, array=NULL, k = NULL,
Chrom=NULL, show="average", weights=NULL,
model.averaging = TRUE, cex=1,
smoother=FALSE, \ldots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{any of RJaCGH, RJaCGH.Chrom, RJaCGH.Genome, RJaCGH.array objects}
  \item{array}{Name of the array to be plotted. If NULL,
    the weigthed average of the arrays is plotted.}
  \item{k}{Model to plot (i.e., number of hidden states). If NULL,
    the most visited is taken (only if a single array is plotted
    and a single model (genome/chromosome).}
  \item{Chrom}{Chromosome to be plotted (only if a single array
    is plotted and a different model for each chromosome has
    been fitted.}
  \item{show}{one of "average" or "frequency". See details.}
  \item{weights}{vector of weights for each array. Must have the length
    of the number of arrays. If NULL, the weights are uniform.}
  \item{model.averaging}{if TRUE, \code{\link{modelAveraging}} is
    performed. If FALSE, a call to \code{\link{states}} is made to get
    hidden state sequence.}
  \item{cex}{A numerical value giving the amount by which plotting text
          and symbols should be scaled relative to the default.}.
  \item{smoother}{Logical. Smoothed means by model averaging.}
  \item{...}{additional arguments passed to plot.}
}
\details{
  Depending on the object and the parameters passed,
  a different plot is drawed:
  If \code{array} in the case of a single model to all genome or
  \code{array} and \code{Chrom} in the case of a different model to
  each chromosome are passed, a panel with 5 subplots is
  returned. The first one is a barplot with the posterior distribution
  of the number of hidden states. The second and third are a density
  plot of the posterior distribution of means and variances. The four
  one is the probability of staying in the same hidden state, as
  returned by \code{\link{plot.Q.NH}}, and the last one shows the
  original observations colored by thir hidden state and the probability
  of being in that hidden state.

  On every plot, the 'Normal' state is coloured black. The 'Gain' states
  are red and the 'Loss' ones green.

  If \code{array} is NULL and \code{show} is 'average', the last
  one of the plots is drawn, but the hidden state sequence and its
  probability is computed averaging on all the arrays with weights
  according to \code{weights} vector. If \code{show} is 'frequency',
  again the  last plot is drawn, but the percentage of arrays in
  which every gene  is Gain/Lost is shown, weighted by the
  \code{weights} vector.

  If \code{smoother} is TRUE, the smoothed mean is drawn. See
  \code{\link{smoothMeans}}, except when
  \code{show} is 'frequency'.
}
\value{
  A plot.
}
\references{
    Rueda OM, Diaz-Uriarte R.
  Flexible and Accurate Detection of Genomic Copy-Number Changes from
  aCGH.
  PLoS Comput Biol. 2007;3(6):e122
}
\author{Oscar M. Rueda and Ramon Diaz Uriarte}
\seealso{\code{\link{RJaCGH}},
  \code{\link{smoothMeans}},
  \code{\link{summary.RJaCGH}}, \code{\link{modelAveraging}},
  \code{\link{states}}, \code{\link{trace.plot}}
%   \code{\link{gelman.rubin.plot}}, \code{\link{collapseChain}}
}
\examples{

y <- c(rnorm(100, 0, 1), rnorm(10, -3, 1), rnorm(20, 3, 1), rnorm(100,
                                                                  0, 1))
Pos <- round(runif(230))
Pos <- cumsum(Pos)
Chrom <- rep(1:23, rep(10, 23))
jp <- list(sigma.tau.mu=rep(0.5, 4), sigma.tau.sigma.2=rep(0.3, 4),
           sigma.tau.beta=rep(0.7, 4), tau.split.mu=0.5, tau.split.beta=0.5)
fit.Chrom <- RJaCGH(y=y, Pos=Pos, Chrom=Chrom, model="Chrom",
                    burnin=100, TOT=1000, jump.parameters=jp, k.max=4)
fit.Genom <- RJaCGH(y=y, Pos=Pos, Chrom=Chrom, model="Genome", burnin=100, TOT=1000, jump.parameters=jp, k.max=4)
fit.none <- RJaCGH(y=y, Pos=Pos, Chrom=NULL, model="None",
                   burnin=100, TOT=1000, jump.parameters=jp, k.max=4)

plot(fit.Chrom)
plot(fit.Chrom, array="array1")
plot(fit.Genom)
plot(fit.none)


y2 <- c(rnorm(100, 0, 1), rnorm(10, -3, 1), rnorm(20, 3, 1),
        rnorm(100, 0, 1))

ya <- cbind(y, y2)

fit.Chrom.array <- RJaCGH(y=ya, Pos=Pos, Chrom=Chrom, model="Chrom",
                          burnin=100, TOT=1000, jump.parameters=jp, k.max=4)
fit.Genom.array <- RJaCGH(y=ya, Pos=Pos, Chrom=Chrom, model="Genome", burnin=100, TOT=1000, jump.parameters=jp, k.max=4)
fit.none.array <- RJaCGH(y=ya, Pos=Pos, Chrom=NULL, model="None",
                         burnin=100, TOT=1000, jump.parameters=jp, k.max=4)

plot(fit.Chrom.array)
plot(fit.Genom.array)
plot(fit.none.array)

}
\keyword{models}% at least one, from doc/KEYWORDS

