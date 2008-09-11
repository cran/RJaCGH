\name{plot.RJaCGH}
\alias{plot.RJaCGH}
\alias{plot.RJaCGH.Chrom}
\alias{plot.RJaCGH.Genome}
\alias{plot.RJaCGH.array}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{'plot' method for RJaCGH objects}
\description{
A plot is drawed with information from the fit of a RJaCGH object.
}
\usage{
\method{plot}{RJaCGH}(x, k = NULL, model.averaging = TRUE, cex=1,
smoother=FALSE, \ldots)
\method{plot}{RJaCGH.Chrom}(x, Chrom="Genome", model.averaging=TRUE,
cex=1, k=NULL, smoother=FALSE, \ldots)
\method{plot}{RJaCGH.Genome}(x, k=NULL, model.averaging=TRUE, cex=1,
smoother=FALSE, \ldots)
\method{plot}{RJaCGH.array}(x, show="frequency", weights=NULL, cex=1,
smoother=FALSE, \ldots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{any of RJaCGH, RJaCGH.Chrom, RJaCGH.Genome, RJaCGH.array objects}
  \item{k}{Model to plot (i.e., number of hidden states). If NULL,
    the most visited is taken.}
  \item{cex}{A numerical value giving the amount by which plotting text
          and symbols should be scaled relative to the default.}.
  \item{model.averaging}{if TRUE, \code{\link{model.averaging}} is
    performed. If FALSE, a call to \code{link{states}} is made to get
    hidden state sequence.}
  \item{Chrom}{For RJaCGH.chrom objects. If "Genome", a whole plot of
    the genome is shown. If a chromosome number, olnly the model
    corresponding to that chromosome. See details.}
  \item{show}{one of "average" or "frequency"}
  \item{weights}{vector of weights for each array. Must have the length
    of the number of arrays. If NULL, the weights are uniform.}
  \item{smoother}{Logical. Smoothed means by model averaging.}
  \item{...}{additional arguments passed to plot.}
}
\details{
  Depending on the object, a different plot is drawed:
  If the object is RJaCGH or RJaCGH.Genome, a panel with 5 subplots is
  returned. The first one is a barplot with the posterior distribution
  of the number of hidden states. The second and third are a density
  plot of the posterior distribution of means and variances. The four
  one is the probability of staying in the same hidden state, as
  returned by \code{\link{plot.Q.NH}}, and the last one shows the
  original observations colored by thir hidden state and the probability
  of being in that hidden state.

  
  On every plot, the 'Normal' state is coloured black. The 'Gain' states
  are red and the 'Loss' ones green.


  
  If the object is 'RJaCGH.Chrom', only the last one of the plots are
  shown. For the complete panel, call the method on every chromosome.

  If the object is 'RJaCGH.array' and \code{show} is 'average', the last
  one of the plots is drawn, but the hidden state sequence and its
  probability is computed averaging on all the arrays with weights
  according to \code{weights} vector. If \code{show} is 'frequency',
  again the  last plot is drawn, but the percentage of arrays in
  which every gene  is Gain/Lost is shown, weighted by the
  \code{weights} vector.

  If \code{smoother} is TRUE, the smoothed mean is drawn. See
  \code{\link{smoothMeans}} (except in \code{plot.RJaCGH.array} when
  method is not 'averaging').
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
  \code{\link{summary.RJaCGH}}, \code{\link{model.averaging}},
  \code{\link{states}}, \code{\link{trace.plot}},
  \code{\link{gelman.rubin.plot}}, \code{\link{collapseChain}}}

\examples{

y <- c(rnorm(100, 0, 1), rnorm(10, -3, 1), rnorm(20, 3, 1), rnorm(100,
                                                                  0, 1))
% this does not work Pos <- round(runif(230))
Pos <- sample(1:10, 230, replace=TRUE)
Pos <- cumsum(Pos)
Chrom <- rep(1:23, rep(10, 23))
jp <- list(sigma.tau.mu=rep(0.5, 4), sigma.tau.sigma.2=rep(0.3, 4),
           sigma.tau.beta=rep(0.7, 4), tau.split.mu=0.5)
fit.Chrom <- RJaCGH(y=y, Pos=Pos, Chrom=Chrom, model="Chrom",
                    burnin=100, TOT=1000, jump.parameters=jp, k.max=4)
fit.Genom <- RJaCGH(y=y, Pos=Pos, Chrom=Chrom, model="Genome", burnin=100, TOT=1000, jump.parameters=jp, k.max=4)
fit.none <- RJaCGH(y=y, Pos=Pos, Chrom=NULL, model="None",
                   burnin=100, TOT=1000, jump.parameters=jp, k.max=4)

plot(fit.Chrom)
plot(fit.Chrom[[1]])
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

