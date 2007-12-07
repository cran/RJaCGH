\name{plot.pREC_S.RJaCGH.array}
\alias{plot.pREC_S.RJaCGH.array}
\alias{plot.pREC_S.RJaCGH.array.Chrom}
\alias{plot.pREC_S.RJaCGH.array.genome}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Plot number of probes shared by pairs of arrays }
\description{
  An image plot showing the results of \code{\link{pREC_S}} on a
  group of arrays. 
}
\usage{
\method{plot}{pREC_S.RJaCGH.array}(x, array.labels = NULL, col = NULL,
breaks = NULL, stats = TRUE, dend=TRUE, \ldots)
\method{plot}{pREC_S.RJaCGH.array.Chrom}(x, array.labels = NULL,
Chrom = NULL, stats = TRUE, col = NULL,
breaks = NULL, dend=TRUE, \ldots)
\method{plot}{pREC_S.RJaCGH.array.genome}(x, array.labels = NULL,
Chrom=NULL, stats = TRUE, col = NULL,
breaks = NULL, dend=TRUE, \ldots)


}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{An object of class \code{\link{pREC_S.RJaCGH.array}},
    \code{\link{pREC_S.RJaCGH.array.Chrom}} or
    \code{\link{pREC_S.RJaCGH.array.genome}}}
  \item{array.labels}{A vector for alternative labels for the arrays.}
  \item{Chrom}{Chromosome to plot. If \code{NULL}, all chromosomes
    are plotted.}
  \item{stats}{Logical. If \code{TRUE}, over every cell the number of common probes and
    the mean length is printed.}
  \item{col}{A vector of color codes for the image plot.}
  \item{breaks}{Breakpoints for the code color. Must be a vector of
    length \code{length(col) + 1}}
  \item{dend}{Logical. If \code{TRUE}, a clustering of arrays is performed with
    \code{\link{hclust}} and arrays reordered.}
  \item{\dots}{Additional arguments passed to \code{\link{image}}}
}
\details{
  First, the number of probes shared by every pair of arrays and their
  mean length is computed.
  The plot consists of a square with as many rows and as many
  columns as the number of arrays are. The more altered probes
  two arrays share the brighter the color is. The diagonals
  are turned off to improve the visibility of the groups. If
  \code{dend} is \code{TRUE}, a hierarchical clustering on arrays
  is performed based on the dissimilarity measure defined as:
  $1 - (inc.mat / max(inc.mat)$ where \code{inc.mat} is the matrix with
  the number of arrays shared by every pair of arrays. Then a dendrogram
  is plotted and the arrays are reordered.
  The diagonals of the plot are turned off to improve the perception of
  the relationships between arrays.
  Note that the number of probes shared depends on the parameters passed
  to \code{\link{pREC_S}}, such as the probability threshold
  \code{p} and the minimum number of arrays requiered to form a
  region \code{freq.array}.
}
\value{
  A list with elements
  \item{probes}{Matrix with the number of probes shared by every pair of
  arrays.}
  \item{length}{Matrix with the mean length of probe shared by every
    pair of arrays.}
}

\references{
  Rueda OM, Diaz-Uriarte R.
  Flexible and Accurate Detection of Genomic Copy-Number Changes from
  aCGH.
  PLoS Comput Biol. 2007;3(6):e122
}

\author{Oscar M. Rueda and Ramon Diaz Uriarte}

\seealso{\code{\link{pREC_S}}, \code{\link{hclust}}, \code{\link{image}}}
\examples{
y <- c(rnorm(100, 0, 1), rnorm(10, -3, 1), rnorm(20, 3, 1),
       rnorm(100,0, 1)) 
Pos <- sample(x=1:500, size=230, replace=TRUE)
Pos <- cumsum(Pos)
Chrom <- rep(1:23, rep(10, 23))

jp <- list(sigma.tau.mu=rep(0.05, 4), sigma.tau.sigma.2=rep(0.03, 4),
           sigma.tau.beta=rep(0.07, 4), tau.split.mu=0.1, tau.split.beta=0.1)

z <- c(rnorm(110, 0, 1), rnorm(20, 3, 1),
       rnorm(100,0, 1)) 
zz <- c(rnorm(90, 0, 1), rnorm(40, 3, 1),
       rnorm(100,0, 1)) 

fit.array.genome <- RJaCGH(y=cbind(y,z,zz),
Pos=Pos, Chrom=Chrom, model="genome",
burnin=1000, TOT=1000, jump.parameters=jp, k.max = 4)
Reg1 <- pREC_S(fit.array.genome, p=0.4, freq.array=2,
alteration="Gain")
plot(Reg1)
}
\keyword{models}% at least one, from doc/KEYWORDS


