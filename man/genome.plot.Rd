\name{genome.plot}
\alias{genome.plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Plot of the genome with probabilities of alteration. }
\description{
Plot of the genome showing, with a color key, the marginal probability
of every gene of alteration.
}
\usage{
genome.plot(obj, col = NULL, breakpoints = NULL, legend.pos=NULL,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{ An object of class RJaCGH.Chrom, RJaCGH.genome or RJaCGH.array. }
  \item{col}{ A vector of length \texttt{k} for the color of every range
    of probabilities of alteration, starting from loss to gain. }
  \item{breakpoints}{ A vector of length \texttt{k-1} for the
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
  Oscar M. Rueda and Ramon Diaz Uriarte. A flexible, accurate and
  extensible statistical method for detecting genomic copy-number
  changes. \url{http://biostats.bepress.com/cobra/ps/art9/}.
  {http://biostats.bepress.com/cobra/ps/art9/}.
}
\author{Oscar M. Rueda and Ramon Diaz Uriarte}

\note{The positions of the genes should be relative to the chromosome
  for the plot to make sense.}

\examples{
data(snijders)
y <- gm13330$LogRatio[!is.na(gm13330$LogRatio)]
Pos <- gm13330$PosBase[!is.na(gm13330$LogRatio)]
Chrom <- gm13330$Chromosome[!is.na(gm13330$LogRatio)]

jp <- list(sigma.tau.mu=rep(0.05, 4), sigma.tau.sigma.2=rep(0.03, 4),
           sigma.tau.beta=rep(0.07, 4), tau.split.mu=0.1, tau.split.beta=0.1)
fit.genome <- RJaCGH(y=y, Pos=Pos, Chrom=Chrom, model="genome",
burnin=1000, TOT=1000, jump.parameters=jp, k.max = 4)
genome.plot(fit.genome)
genome.plot(fit.genome, col=c(3, 1, 2), breakpoints=c(-0.5, 0.5))
}
\keyword{ models }% __ONLY ONE__ keyword per line
