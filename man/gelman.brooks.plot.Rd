\name{gelman.brooks.plot}
\alias{gelman.brooks.plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{gelman-brooks plot for 'RJaCGH' objects}
\description{
A plot to show the convergence of several parallel chains, as described
in Brooks and Gelman, 1998.
}
\usage{
gelman.brooks.plot(obj, bin = 1000, array = NULL, Chrom = NULL, k = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{a list containing several parallel chains; that is objects
    of any of RJaCGH, RJaCGH.Chrom, RJaCGH.genome, RJaCGH.array classes
    (obviously, all of the same class).}
  \item{bin}{Number of observations taken in every subchain}
  \item{array}{if \code{obj} is 'RJaCGH.array', the name of the array to
  plot must be given.}
  \item{Chrom}{if \code{obj} is 'RJaCGH.Chrom', the number of the chromosome to
  plot must be given.}
  \item{k}{Model to monitorize (i.e., number of hidden states). If NULL,
    the most visited is taken.}
}
\details{
  As described in the references, for every \code{bin} runs of the
  chain, the R value of \code{k}, \code{mu}, \code{sigma.2} and
  \code{beta} are computed.
  Please note that the input to this function must be a list of parallel
  chains. If the convergence diagnostics returned by it are
  satisfactory, then one can run \code{collapseChain} to join them all,
  but one can not collapse all the chains and later run
  \code{gelman.brooks.plot} to recover the diagnostics, because the list
  of chains has been merged into one object.
  
}
\value{
  A plot with four panels showing the R values is drawn.
  The first one shows the R value for the number of hidden states.
  The second one shows the R values for the means of the most visited
  model (there will be a R value for each mean).
  The third and the four show the R values for variances and betas of
  the most visited model.
  Besides, a list is returned with components
  \item{k}{R values for the number of states}
  \item{mu}{R values for the means of the states}
  \item{sigma.2}{R values for the variances of the states}
  \item{beta}{R values for the beta parameters}
}
\references{Brooks, S.P. and Gelman, A. (1998). General Methods for
  Monitoring convergence of iterative simulations. Journal of
  Computational and Graphical Statistics. p434-455.

  
  Rueda OM, Diaz-Uriarte R.
  Flexible and Accurate Detection of Genomic Copy-Number Changes from
  aCGH.
  PLoS Comput Biol. 2007;3(6):e122
}  
\author{Oscar M. Rueda and Ramon Diaz Uriarte}

\seealso{\code{\link{RJaCGH}},
  \code{\link{summary.RJaCGH}}, \code{\link{model.averaging}},
  \code{\link{plot.RJaCGH}}, \code{\link{states}},
  \code{\link{trace.plot}}, \code{\link{collapseChain}}}

\examples{
y <- c(rnorm(100, 0, 1), rnorm(10, -3, 1), rnorm(20, 3, 1),
       rnorm(100,0, 1)) 
Pos <- sample(x=1:500, size=230, replace=TRUE)
Pos <- cumsum(Pos)
Chrom <- rep(1:23, rep(10, 23))

jp <- list(sigma.tau.mu=rep(0.5, 4), sigma.tau.sigma.2=rep(0.3, 4),
           sigma.tau.beta=rep(0.7, 4), tau.split.mu=0.5, tau.split.beta=0.5)

fit.genome <- list()
for (i in 1:4) {
fit.genome[[i]] <- RJaCGH(y=y, Pos=Pos, Chrom=Chrom, model="genome",
burnin=10, TOT=1000, jump.parameters=jp, k.max = 4)
}

\dontrun{gelman.brooks.plot(fit.genome)}
}
\keyword{models}% at least one, from doc/KEYWORDS
 


