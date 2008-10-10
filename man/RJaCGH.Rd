\name{RJaCGH}
\alias{RJaCGH}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Reversible Jump MCMC for the analysis of arrays of CGH }
\description{
This function fits a non-homogeneous hidden Markov model to CGH data
through bayesian methods and Reversible Jump Markov chain Montecarlo.
}
\usage{
RJaCGH(y, Chrom = NULL, Start=NULL, End=NULL, Pos = NULL,
       Dist=NULL, probe.names=NULL, maxVar=NULL, model = "Genome",
       var.equal=TRUE, max.dist=NULL, normal.reference=0,
       window = NULL, burnin = 10000, TOT =10000,
       k.max = 6, stat = NULL, mu.alfa = NULL, mu.beta = NULL,
       s1=NULL, s2=NULL, init.mu=NULL, init.sigma.2=NULL, init.beta=NULL,
       prob.k = NULL, jump.parameters=list(),
       start.k = NULL, RJ=TRUE, NC = 1, deltaT = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{Vector with Log Ratio observations.}
  \item{Chrom}{Vector with Chromosome indicator.}
  \item{Start}{Vector with start positions of the probes.}
  \item{End}{Vector with end positions of the probes.}
  \item{Pos}{Vector with Positions of every gene. They can be absolute
    to the genome or relative to the chromosome. They should be ordered
    within every chromosome. This is, the arrays must be ordered by their
    positions in the genome. They must be integers. Positions can be
    specified with \code{Start}, \code{End} or \code{Pos}, but only in
    one of the two ways.}
  \item{Dist}{Optional vector of distances between genes. It should be a vector
    of length \code{length(y)-1}. Note that when \code{Chrom} is not NULL,
    every last value of every Chromosome is not used.}
  \item{probe.names}{Character vector with the number of the probes.}
  \item{maxVar}{Maximum value for the variance of the states. If
    \code{NULL}, the range of the data is chosen.}
  \item{model}{if \code{model}="Genome", the same model is fitted for
    the whole genome. If \code{model}="Chrom", a different model is
    fitted for each chromosome.}
  \item{var.equal}{Logical. If \code{TRUE} the variances of the hidden
    states are restricted to be the same.}
  \item{max.dist}{maximal distance between spots. When two spots have
    a distance between them as far or further than \code{max.dist}, they
    are considered independent. That is, the state of that spot
    does not affect the state of the other. If \code{NULL} (the default) the
    maximum \code{Dist} or maximum difference in \code{Pos} is taken.}
  \item{normal.reference}{The value considered as the mean of the normal
    state. See details. By default is \code{0}.}
  \item{window}{Multiplier of the standard deviation of the data to
    determine the width of the normal state. See details. Default (window
    = NULL) is 1.}
%   \item{normal.ref.percentile}{Percentage for the relabelling of
%     states. See details. by default is 0.95.}
  \item{burnin}{Number of burn-in iterations in the Markov Chain}
  \item{TOT}{Number of iterations after the burn-in}
  \item{k.max}{Maximum number of hidden states to fit.}
  \item{stat}{Initial Distribution for the hidden states. Must be a
    vector of size 1 + 2 + ... +\code{k.max}. If \code{NULL}, it is assumed a
    uniform distribution for every model.}
  \item{mu.alfa}{Hyperparameter. See details}
  \item{mu.beta}{Hyperparameter. See details}
  \item{prob.k}{Hyperparameter. See details}
  \item{jump.parameters}{List with the parameters for the MCMC
    jumps. See details.}
  \item{start.k}{Initial number of states. if \code{NULL}, a random draw from
    \code{prob.k} is chosen.}
  \item{RJ}{Logical. If \code{TRUE}, Reversible Jump is performed.
    If not, MCMC
    over a fixed number of hidden states. Note that if \code{FALSE}, most
    of the methods for extracting information won't work.}
  \item{s1}{Standard deviation for the creation of a new mean in the
    birth move (first attempt). If \code{NULL}, it is the prior \code{mu.beta}.
  }
  \item{s2}{Standard deviation for the creation of a new mean in the
    birth move (second attempt). If \code{NULL}, it is the prior \code{mu.beta}.}
  \item{init.mu}{Starting values for \code{mu} for the chain. See details}
  \item{init.sigma.2}{Starting values for \code{sigma.2} for the chain. See details}
  \item{init.beta}{Starting values for \code{beta} for the chain. See details}
  \item{NC}{Number of coupled parallel chains. See details}
  \item{deltaT}{Heat parameter for tempering the parallel chains. See details}
 
%   \item{auto.label}{If not \code{NULL}, should be the minimum proportion of
%   observations labeled as 'Normal'. See details.}
}
\details{
  RJaCGH fits the following bayesian model: There is a priori
  distribution for the number of hidden states (different copy numbers)
  as stated by \code{prob.k}. If \code{NULL}, a uniform distribution between 1
  and \code{k.max} is used.

  
  The hidden states follow a normal distribution which mean (\code{mu})
  follows
  itself a normal distribution with mean
  \code{mu.alfa} and stdev \code{mu.beta}. If \code{NULL}, these are the
  median of the data and the range. The square
  root of the variance (\code{sigma.2})of the hidden states
  follows a uniform distribution between $0$ and \code{maxVar}.

  The model for the transition matrix is based on a random matrix
  \code{beta} whose diagonal is zero. The transition matrix, \code{Q},
  has the form:
    Q[i,j] = exp(-beta[i,j] + beta[i,j]*x) / sum(i,.) {exp(-beta[i,.] +
      beta[i,.]*x}

    
  The prior distribution for \code{beta} is gamma with parameters 1, 1. 
  The \code{x} are the distances between positions, normalized to lay
  between zero and 1 (\code{x=diff(Pos) / max(diff(Pos))})

  
  RJaCGH performs Markov Chain MonteCarlo with Reversible Jump to sample
  for the posterior distribution.
  \code{NC} sets the number of chains from we will sample in
  parallel. Each of them is tempered in order to escape from local
  maximum. The temper parameter is \code{deltaT}, and it can be
  a value greater than $0$. each chain is tempered according to:

  $ 1 / (1 + \code{deltaT} * \code{NC}) $
  
  Every sweep is performed for all chains and has 3 steps plus
  another one common for all of them:

  
  1.- A Metropolis-Hastings move is used to update, for a fixed number
  of hidden states, \code{mu}, \code{sigma.2} and \code{beta}. A
  symmetric proposal with a normal distribution and standard deviation
  \code{sigma.tau.mu}, \code{sigma.tau.sigma.2} and
  \code{sigma.tau.beta} is sampled.

  
  2.- A transdimensional move is chosen, between birth (a new hidden
  state is sampled from the prior) or death (an existing hidden state is
  erased). Both moves are tried using delayed rejection. That is, if
  the move is rejected, is given another try. The means for the new
  state are drawn for the priors, but the standard deviation can be
  set for the two stages with parameters \code{s1} and \code{s2}.

  
  3.- Another transdimensional move is performed; an split move (divide
  an existing state in two) or a combine move (join two adjacent
  states). The length of the split is sampled from a normal distribution
  with standard deviation \code{tau.split.mu} for the \code{mu} and
  \code{tau.split.beta} for \code{beta}.

  4.- If \code{NC} is greater than 1, a swap move is tried to exchange
  information from two of the coupled parallel chains.
  
  \code{jump.parameters} must be a list with the parameters for the
  moves. It must have components \code{sigma.tau.mu},
  \code{sigma.tau.sigma.2}, \code{sigma.tau.beta} These are vectors of
  length \code{k.max}. \code{tau.split.mu} is a vector of
  length 1. If any of them is NULL, a call to the internal function
  \code{get.jump()} is made to find 'good' values.

  A relabelling of hidden states is performed to match biological
  states. See details in \code{\link{relabelStates}}.

  
  The initial values of the chain are drawn from an overdispersed
  distribution. One can start the chain in a given point with
  the parameters \code{start.k} (model to start from $1, \ldots,
  \code{max.k}$) and the initial values \code{init.mu},
  \code{init.sigma.2} (vectors of dimension \code{start.k}) and
  \code{init.beta} (matrix of positive values with \code{start.k})
  rows and \code{start.k}) columns. The diagonal must be zero.
  
%   The states that have the \code{normal.reference} value
%   inside a \code{normal.ref.percentile}\% probability interval
%   based on a normal distribution with means the median of \code{mu}
%   and sd the square root of the median of \code{sigma.2} are labelled as
%   'Normal'. If no state is close enough to \code{normal.reference} then
%   there will not be a normal state. Bear this in mind for
%   normalization issues.
%   If \code{auto.label} is not null, closest states to 'Normal' are also
%   labelled as 'Normal' until a proportion of \code{auto.label} is
%   reached. Please note that the default value is 0.60, so at least the
%   60\% of the observations will be labelled as 'Normal'.
%   If this laeblling is not satisfactory, you can relabel with
%   \code{\link{relabelStates}}.
    
}
\value{
  The object returned follows a hierarchy:
  
  If y is a matrix or data.frame (i.e., several arrays), an object of
  class \code{RJaCGH.array} is returned, with components:
  \item{[[]]}{A list with an object of corresponding class (see below) for
    every array.}
  \item{array.names}{Vector with the names of the arrays.}
  If \code{model} is "Genome", an object of class \code{RJaCGH.Genome}
  is returned, with components:
  \item{[[]]}{a list with as many objects as k.max, with the fits.}
  \item{k}{sequence of number of hidden states sampled.}
  \item{prob.b}{Number of birth moves performed in first and second
    proposal.(Includes burn-in.}
  \item{prob.d}{Number of death moves performed in first and second
    proposal.(Includes burn-in.}
  \item{prob.s}{Number of split moves performed (Includes burn-in.}
  \item{prob.c}{Number of combine moves performed (Includes burn-in.}
  \item{prob.e}{Number of exchange (swap) moves performed between
    tempered chains and with cool (main) chain.}
  \item{y}{y vector.}
  \item{Pos}{Pos vector.}
  \item{model}{model.}
  \item{Chrom}{Chromosome vector.}
  \item{x}{x vector of distances between genes.}
  \item{viterbi}{A list with as many components as chromosomes. For each
  chromosome, a list, with at least three components: gzipped\_sequence:
  the compacted and gzipped sequence of states; num\_sequences: the
  number of sequences in that compacted set of sequences;
  sum\_mcmc\_iter: the number of times a viterbi sequence was
  obtained. All these are only of use for calls to the pREC funcitons
  (\code{\link{pREC_A}} and \code{\link{pREC_S}}).}

  If \code{model} is "Chrom", an object of class \code{RJaCGH.Chrom} is
  returned, with the following components: \item{[[]]}{a list with as
  many components as chromosomes, of class \code{RJaCGH} (See below).}
  \item{Pos}{Pos vector.}  \item{Start}{Start positions.}
  \item{End}{End positions.}  \item{probe.names}{Names of the probes.}
  \item{model}{model.}  \item{Chrom}{Chromosome vector.}
  \item{viterbi}{Identical structure as above.}  If no \code{model} was
  specified and no \code{Chrom} was given, an object of class
  \code{RJaCGH} is returned, with components \code{k}, \code{prob.b},
  \code{prob.d}, \code{prob.s}, \code{prob.c}, \code{y}, \code{Pos},
  \code{x}, \code{viterbi} (but as there are no chroms, viterbi has a
  one component list ---analogous to having a single chromosome), as
  described before, plus a list with as many components of number of max
  hidden states fitted.  The length of \code{k} equals aproximately $2$
  times \code{TOT}, because in every sweep of the algorithm there are
  two tries to jump between models, so two times to explore the
  probability of the number of hidden states.  For every hidden markov
  model fitted, a list is returned with components:
  
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
  The number of rows of components \code{mu}, \code{sigma.2} and
  \code{beta} is random, because it depends on the number of times
  a particular model is visited and on the number of moves between
  models, because when we visit a new model we also explore the space
  of its means, variances and parameters of its transition functions.
}

\note{
  The data must be ordered by chromosome and within chromosome by position.
  }
\references{
  Rueda OM, Diaz-Uriarte R.
  Flexible and Accurate Detection of Genomic Copy-Number Changes from
  aCGH.
  PLoS Comput Biol. 2007;3(6):e122
  
  Cappe, Moulines and Ryden, 2005. Inference in Hidden Markov
  Models. Springer.
  
  Green, P.J. (1995) Reversible Jump Markov Chain Monte Carlo
  computation and Bayesian model determination. Biometrika, 82,
  711-732.

  Green, P.J. and Antonietta, M. (2001) Delayed Rejection in Reversible
  Jump Metropolis Hastings. Biometrika, 88 (4), 1035-1053.

  Geyer, C. J. (1991). Markov Chain Monte Carlo Maximum Likelihood.
  Proceedings of the 23th Symposium on the Interface, 156-163.
  
}
\author{Oscar M. Rueda and Ramon Diaz Uriarte}

\seealso{\code{\link{summary.RJaCGH}},
  \code{\link{states}}, \code{\link{model.averaging}},
  \code{\link{plot.RJaCGH}}, \code{\link{trace.plot}},
  \code{\link{gelman.rubin.plot}}, \code{\link{collapseChain}},
  \code{\link{relabelStates}}, \code{\link{pREC_A}},
  \code{\link{pREC_S}}
}

\examples{
y <- c(rnorm(100, 0, 1), rnorm(10, -3, 1), rnorm(20, 3, 1),
       rnorm(100,0, 1)) 
Pos <- sample(x=1:500, size=230, replace=TRUE)
Pos <- cumsum(Pos)
Chrom <- rep(1:23, rep(10, 23))

jp <- list(sigma.tau.mu=rep(0.05, 4), sigma.tau.sigma.2=rep(0.03, 4),
           sigma.tau.beta=rep(0.07, 4), tau.split.mu=0.1)

fit.chrom <- RJaCGH(y=y, Pos=Pos, Chrom=Chrom, model="Chrom",
                    burnin=10, TOT=1000, k.max = 4,
                    jump.parameters=jp)
##RJ results for chromosome 5
table(fit.chrom[[5]]$k)
fit.genome <- RJaCGH(y=y, Pos=Pos, Chrom=Chrom, model="Genome",
burnin=100, TOT=1000, jump.parameters=jp, k.max = 4)
## Results for the model with 3 states:
apply(fit.genome[[3]]$mu, 2, summary)
apply(fit.genome[[3]]$sigma.2, 2, summary)
apply(fit.genome[[3]]$beta, c(1,2), summary)
}
\keyword{models}% at least one, from doc/KEYWORDS

