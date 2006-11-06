\name{plot.RJaCGH}
\alias{plot.RJaCGH}
\alias{plot.RJaCGH.Chrom}
\alias{plot.RJaCGH.genome}
\alias{plot.RJaCGH.array}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{'plot' method for RJaCGH objects}
\description{
A plot is drawed with information from the fit of a RJaCGH object.
}
\usage{
plot.RJaCGH(x, k = NULL, model.averaging = TRUE, cex=1, ...)
plot.RJaCGH.Chrom(x, Chrom="genome", model.averaging=TRUE, cex=1, ...)
plot.RJaCGH.genome(x, k=NULL, model.averaging=TRUE, cex=1, ...)
plot.RJaCGH.array(x, method="averaging", weights=NULL, cex=1, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{any of RJaCGH, RJaCGH.Chrom, RJaCGH.genome, RJaCGH.array objects}
  \item{k}{Model to plot (i.e., number of hidden states). If NULL,
    the most visited is taken.}
  \item{cex}{A numerical value giving the amount by which plotting text
          and symbols should be scaled relative to the default.}.
  \item{model.averaging}{if TRUE, \code{\link{model.averaging}} is
    performed. If FALSE, a call to \code{link{states}} is made to get
    hidden state sequence.}
  \item{Chrom}{For RJaCGH.chrom objects. If "genome", a whole plot of
    the genome is shown. If a chromosome number, olnly the model
    corresponding to that chromosome. See details.}
  \item{method}{one of "averaging" or "other"}
  \item{weights}{vector of weights for each array. Must have the length
    of the number of arrays. If NULL, the weights are uniform.}
  \item{...}{additional arguments passed to plot.}
}
\details{
  Depending on the object, a different plot is drawed:
  If the object is RJaCGH or RJaCGH.genome, a panel with 5 subplots is
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

  If the object is 'RJaCGH.array' and method is 'averaging', the last
  one of the plots is drawn, but the hidden state sequence and its
  probability is computed averaging on all the arrays with weights
  according to \code{weights} vector. If method is 'other', again the
  last plot is drawn, but the percentage of arrays in which every gene
  is Gain/Lost is shown, weighted by the \code{weights} vector.
}
\value{
  A plot.
}
\references{Oscar Rueda and Ramon Diaz Uriarte, in prep.}
\author{Oscar Rueda and Ramon Diaz Uriarte}
 \seealso{\code{\link{RJaCGH}},
  \code{\link{summary.RJaCGH}}, \code{\link{model.averaging}},
  \code{\link{states}}, \code{\link{trace.plot}},
  \code{\link{gelman.brooks.plot}}, \code{\link{collapseChain}}}

\examples{
\dontrun{y <- c(rnorm(100, 0, 1), rnorm(10, -3, 1), rnorm(20, 3, 1), rnorm(100,
0, 1))
Pos <- runif(230)
Pos <- cumsum(Pos)
Chrom <- rep(1:23, rep(10, 23))
jp <- list(sigma.tau.mu=rep(0.5, 4), sigma.tau.sigma.2=rep(0.3, 4),
           sigma.tau.beta=rep(0.7, 4), tau.split.mu=0.5, tau.split.beta=0.5)
fit.Chrom <- RJaCGH(y=y, Pos=Pos, Chrom=Chrom, model="Chrom",
burnin=100, TOT=1000, jump.parameters=jp, k.max=4)
plot(fit.Chrom)
plot(fit.Chrom[[1]])}
}
\keyword{models}% at least one, from doc/KEYWORDS

