\name{genomePlot}
\alias{genomePlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Plot of the genome with probabilities of alteration. }
\description{
Plot of the genome showing, with a color key, the marginal probability
of every gene of alteration.
}
\usage{
genomePlot(obj, array=NULL, weights=NULL, 
col = NULL, breakpoints = NULL, legend.pos=NULL,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{An object of class RJaCGH.Chrom, RJaCGH.Genome or RJaCGH.array. }
  \item{array}{Name of the array to be plotted. If NULL, the
    weigthed average of all is computed.}
  \item{weights}{vector of weights for each array. Must have the length
    of the number of arrays. If NULL, the weights are uniform.}
  \item{col}{ A vector of length \emph{k} for the color of every range
    of probabilities of alteration, starting from loss to gain. }
  \item{breakpoints}{ A vector of length \emph{k-1} for the
    breakpoints of the color key. The corresponding to losses must be
    negative. See example for details.}
  \item{legend.pos}{Position of the legend. Must be a vector with two
    elements; the position of the x and y coordinates. If \code{NULL},
    the legend is placed at the right.}
  \item{...}{Aditional parameters passed to plot.}
}
\details{
If \code{col} and \code{breakpoints} are \code{NULL}, a default
color key is drawn.
}
\value{

  A plot is drawn.
}
\references{
  Rueda OM, Diaz-Uriarte R.
  Flexible and Accurate Detection of Genomic Copy-Number Changes from
  aCGH.
  PLoS Comput Biol. 2007;3(6):e122
}
\author{Oscar M. Rueda and Ramon Diaz Uriarte}

\note{The positions of the genes should be relative to the chromosome
  for the plot to make sense.}

\examples{
data(snijders)
y <- gm13330$LogRatio[!is.na(gm13330$LogRatio)]
Pos <- gm13330$PosBase[!is.na(gm13330$LogRatio)]
Chrom <- gm13330$Chromosome[!is.na(gm13330$LogRatio)]

## Sort positions

for (i in unique(Chrom)) {
if(any(diff(Pos[Chrom==i]) < 0)) {
id <- order(Pos[Chrom==i])
y[Chrom==i] <- y[Chrom==i][id]
Pos[Chrom==i] <- Pos[Chrom==i][id]
}
}

jp <- list(sigma.tau.mu=rep(0.05, 4), sigma.tau.sigma.2=rep(0.03, 4),
           sigma.tau.beta=rep(0.07, 4), tau.split.mu=0.1, tau.split.beta=0.1)
fit.genome <- RJaCGH(y=y, Pos=Pos, Chrom=Chrom, model="Genome",
burnin=1000, TOT=1000, jump.parameters=jp, k.max = 4)
genomePlot(fit.genome)
genomePlot(fit.genome, col=c(3, 1, 2), breakpoints=c(-0.5, 0.5))
}
\keyword{ models }% __ONLY ONE__ keyword per line
