\name{simulateRJaCGH}
\alias{simulateRJaCGH}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Simulate observations form a hidden Markov model with
  non-homogeneous transition probabilities.}
\description{
This function simulates observations from a hidden Markov model with
normal distributed observations and non-homogeneous transition matrix.
}
\usage{
simulateRJaCGH(n, x = NULL, mu, sigma.2, beta, start)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{Number of observations to simulate}
  \item{x}{Distance to the next observation. Must be a vector of size
    n-1 and normalized between zero and one.
    If NULL, a vector of zeros is taken}
  \item{mu}{Vector of means for the hidden states}
  \item{sigma.2}{Vector of variances for the hidden states}
  \item{beta}{\code{beta} parameter of the transition matrix. Must be a square
  matrix with the same size as the number of hidden states.}
  \item{start}{Starting states of the sequence. Must be an integer from
    1 to the number of hidden states.}
}
\details{
Please note that in RJaCGH model, parameter \code{q} is taken as \code{-beta}
}
\value{
A list with components
  \item{states }{Sequence of hidden states}
  \item{y }{Observations}

}
\references{
    Rueda OM, Diaz-Uriarte R.
  Flexible and Accurate Detection of Genomic Copy-Number Changes from
  aCGH.
  PLoS Comput Biol. 2007;3(6):e122
}

\author{Oscar M. Rueda and Ramon Diaz-Uriarte}

\seealso{Q.NH, RJaCGH}
\examples{
beta <- matrix(c(0, 5, 1, 1,  0, 1, 3, 5, 0), 3)
obs <- simulateRJaCGH(n=200, x=rexp(199), mu=c(-3, 0, 3), sigma.2=c(1,1,1),
beta=beta, start=2)
plot(obs$y, col=obs$states)
}
\keyword{models}% at least one, from doc/KEYWORDS

