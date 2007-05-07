\name{RJaCGH}
\alias{RJaCGH}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Reversible Jump MCMC for the analysis of arrays of CGH }
\description{
This function fits a non-homogeneous hidden Markov model to CGH data
through bayesian methods and Reversible Jump Markov chain Montecarlo.
}
\usage{
RJaCGH(y, Chrom = NULL, Pos = NULL, Dist=NULL, model = "genome",
max.dist=NULL, normal.reference=0, normal.ref.percentile=0.95,
burnin = 10000, TOT =10000, k.max = 6, stat = NULL, mu.alfa = NULL,
mu.beta = NULL, ka = NULL,
g = NULL, prob.k = NULL, jump.parameters=list(), start.k = NULL,
RJ=TRUE, auto.label=0.60)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{Vector with Log Ratio observations.}
  \item{Chrom}{Vector with Chromosome indicator. Must be numeric (1 to 23)}
  \item{Pos}{Vector with Positions of every gene. They can be absolute
    to the genome or relative to the chromosome. They should be ordered
    within every chromosome. This is, the arrays must be ordered by their
    positions in the genome. They must be integers.}
  \item{Dist}{Optional vector of distances between genes. It should be a vector
    of length \code{length(y)-1}. Note that when \code{Chrom} is not NULL,
    every last value of every Chromosome is not used.}
  \item{model}{if \code{model}="genome", the same model is fitted for
    the whole genome. If \code{model}="Chrom", a different model is
    fitted for each chromosome.}
  \item{max.dist}{maximal distance between spots. When two spots have
    a distance between them as far or further than \code{max.dist}, they
    are considered independent. That is, the state of that spot
    does not affect the state of the other. If \code{NULL} (the default) the
    maximum \code{Dist} or maximum diference in \code{Pos} is taken.}
  \item{normal.reference}{The value considered as the mean of the normal
    state. See details. By default is \code{0}.}
  \item{normal.ref.percentile}{Percentage for the relabelling of
    states. See details. by default is 0.95.}
  \item{burnin}{Number of burn-in iterations in the Markov Chain}
  \item{TOT}{Number of iterations after the burn-in}
  \item{k.max}{Maximum number of hidden states to fit.}
  \item{stat}{Initial Distribution for the hidden states. Must be a
    vector of size 1 + 2 + ... +\code{k.max}. If NULL, it is assumed a
    uniform distribution for every model.}
  \item{mu.alfa}{Hyperparameter. See details}
  \item{mu.beta}{Hyperparameter. See details}
  \item{ka}{Hyperparameter. See details}
  \item{g}{Hyperparameter. See details}
  \item{prob.k}{Hyperparameter. See details}
  \item{jump.parameters}{List with the parameters for the MCMC
    jumps. See details.}
  \item{start.k}{Initial number of states. if NULL, a random draw from
    \code{prob.k} is chosen.}
  \item{RJ}{Logical. If TRUE, Reversible Jump is performed. If not, MCMC
    over a fixed number of hidden states.}
  \item{auto.label}{If not NULL, should be the minimum proportion of
  observations labeled as 'Normal'. See details.}
}
\details{
  RJaCGH fits the following bayesian model: There is a priori
  distribution for the number of hidden states (different copy numbers)
  as stated by \code{prob.k}. If NULL, a uniform distribution between 1
  and \code{k.max} is used.

  
  The hidden states follow a normal distribution which mean (\code{mu}) follows
  itself a normal distribution with mean
  \code{mu.alfa} and stdev \code{mu.beta}. If NULL, these are the
  median of the data and the range. The variance (\code{sigma.2})of the hidden states
  follows an Inverse Gamma distribution with parameters \code{ka} and
  \code{g}. If NULL, these are \code{2} and \code{diff(range(y))^2 / 50},
  respectively.

  
  The model for the transition matrix is based on a random matrix
  \code{beta} whose diagonal is zero. The transition matrix, \code{Q},
  has the form:
    Q[i,j] = exp(-beta[i,j] + beta[i,j]*x) / sum(i,.) {exp(-beta[i,.] +
      beta[i,.]*x}

    
  The prior distribution for \code{beta} is gamma with parameters 1, 1. 
  The \code{x} are the distances between positions, normalized to lay
  between zero and 1 (\code{x=diff(Pos) / max(diff(Pos))})

  
  RJaCGH performs Markov Chain MonteCarlo with Reversible Jump to sample
  for the posterior distribution. Every sweep has 3 steps:

  
  1.- A Metropolis-Hastings move is used to update, for a fixed number
  of hidden states, \code{mu}, \code{sigma.2} and \code{beta}. A
  symmetric proposal with a normal distribution and standard deviation
  \code{sigma.tau.mu}, \code{sigma.tau.sigma.2} and
  \code{sigma.tau.beta} is sampled.

  
  2.- A transdimensional move is chosen, between birth (a new hidden
  state is sampled from the prior) or death (an existing hidden state is
  erased).

  
  3.- Another transdimensional move is performed; an split move (divide
  an existing state in two) or a combine move (join two adjacent
  states). The length of the split is sampled from a normal distribution
  with standard deviation \code{tau.split.mu} for the \code{mu} and
  \code{tau.split.beta} for \code{beta}.

  
  \code{jump.parameters} must be a list with the parameters for the
  moves. It must have components \code{sigma.tau.mu},
  \code{sigma.tau.sigma.2}, \code{sigma.tau.beta} These are vectors of
  length \code{k.max}. \code{tau.split.mu},  \code{tau.split.beta} are vectors of
  length 1. If any of them is NULL, a call to the internal function
  \code{get.jump()} is made to find 'good' values.

  A relabelling of hidden states is performed to match biological
  states. The states whom have the \code{normal.reference} value
  inside a \code{normal.ref.percentile}\% confidence
  based on a normal distribution with means the median of \code{mu}
  and sd the square root of the median of \code{sigma.2} are labelled as
  'Normal'. If no state is labelled as 'Normal', the closest state to
  \code{normal.reference} is chosen.
  If \code{auto.label} is not null, closest states to 'Normal' are also
  labelled as 'Normal' until a proportion of \code{auto.label} is
  reached. Please note that the default value is 0.60, so at least the
  60\% of the observations will be labelled as 'Normal'.
  
    
}
\value{
  The object returned follows a hierarchy:
  
  If y is a matrix or data.frame (i.e., several arrays), an object of
  class \code{RJaCGH.array} is returned, with components:
  \item{[[]]}{A list with an object of corresponding class (see below) for
    every array.}
  \item{array.names}{Vector with the names of the arrays.}
  If \code{model} is "genome", an object of class \code{RJaCGH.genome}
  is returned, with components:
  \item{[[]]}{a list with as many objects as k.max, with the fits.}
  \item{k}{sequence of number of hidden states sampled.}
  \item{prob.b}{Number of birth moves performed (Includes burn-in.}
  \item{prob.d}{Number of death moves performed (Includes burn-in.}
  \item{prob.s}{Number of split moves performed (Includes burn-in.}
  \item{prob.c}{Number of combine moves performed (Includes burn-in.}
  \item{y}{y vector.}
  \item{Pos}{Pos vector.}
  \item{model}{model.}
  \item{Chrom}{Chromosome vector.}
  \item{x}{x vector of distances between genes.}
  If \code{model} is "Chrom", an object of class \code{RJaCGH.Chrom} is
  returned, with the following components:
  \item{[[]]}{a list with as many components as chromosomes, of class
    \code{RJaCGH} (See below).}
  \item{Pos}{Pos vector.}
  \item{model}{model.}
  \item{Chrom}{Chromosome vector.}
  
  If no \code{model} was specified and no \code{Chrom} was given, an
  object of class \code{RJaCGH} is returned, with components \code{k},
  \code{prob.b}, \code{prob.d}, \code{prob.s}, \code{prob.c}, \code{y},
  \code{Pos}, \code{x} as described before, plus a list with as many
  components of number of max hidden states fitted.
  For every hidden markov model fitted, a list is returned with
  components:
  
  \item{mu}{a matrix with the means sampled}
  \item{sigma.2}{a matrix with the variances sampled}
  \item{beta}{an array of dimension 3 with beta values sampled}
  \item{stat}{vector of initial distribution}
  \item{loglik}{log likelihoods of every MCMC iteration}
  \item{prob.mu}{probability of aceptance of \code{mu} in the
    Metropolis-Hastings step.}
  \item{prob.sigma.2}{probability of aceptance of \code{sigma.2} in the
    Metropolis-Hastings step.}
  \item{prob.beta}{probability of aceptance of \code{beta} in the
    Metropolis-Hastings step.}
  \item{state.labels}{Labels of the biological states.}
  \item{prob.states}{Marginal posterior probabilities of belonging to
    every hidden state.}
}
\references{
  Oscar M. Rueda and Ramon Diaz Uriarte. A flexible, accurate and
  extensible statistical method for detecting genomic copy-number
  changes. \url{http://biostats.bepress.com/cobra/ps/art9/}.
  {http://biostats.bepress.com/cobra/ps/art9/}.
  
  Cappe, Moulines and Ryden, 2005. Inference in Hidden Markov
  Models. Springer.
  
  Green, P.J. (1995) Reversible Jump Markov Chain Monte Carlo
  computation and Bayesian model determination. Biometrika, 82,
  711-732.
}
\author{Oscar M. Rueda and Ramon Diaz Uriarte}

\seealso{\code{\link{summary.RJaCGH}},
  \code{\link{states}}, \code{\link{model.averaging}},
  \code{\link{plot.RJaCGH}}, \code{\link{trace.plot}},
  \code{\link{gelman.brooks.plot}}, \code{\link{collapseChain}}}

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
##RJ results for chromosome 5
table(fit.chrom[[5]]$k)
fit.genome <- RJaCGH(y=y, Pos=Pos, Chrom=Chrom, model="genome",
burnin=100, TOT=1000, jump.parameters=jp, k.max = 4)
## Results for the model with 3 states:
apply(fit.genome[[3]]$mu, 2, summary)
apply(fit.genome[[3]]$sigma.2, 2, summary)
apply(fit.genome[[3]]$beta, c(1,2), summary)
}
\keyword{models}% at least one, from doc/KEYWORDS

