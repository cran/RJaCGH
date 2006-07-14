\name{Q.NH}
\alias{Q.NH}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Transition Matrix for non-homogeneous Hidden Markov Model}
\description{
This function returns the transition matrix for a given distance between genes.
}
\usage{
Q.NH(beta, x, q=-beta)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{beta}{beta parameter of transition matrix. Must be a square
    matrix with dimension equal to the number of hidden states.}
  \item{x}{Distance between genes to compute the transition matrix. Must
    be a scalar.}
    \item{q}{q parameter of transition matrix. Note that in RJaCGH q is
    always -beta (details below). Must be a square matrix with dimension equal to the
    number of hidden states.}
}
\details{
  RJaCGH assumes a non-homogeneous transition matrix with this form: 
  Q[i,j] = exp(-beta[i,j] + beta[i,j]*x) / sum(i,.) {exp(-beta[i,.] +
    beta[i,.]*x}

  
  All beta[i,i] are constrained to be zero, for the model to be
  identifiable.
  
  All beta[i,j] are positive.  This model is chosen for its simplicity
  and because it agrees with 
  biological assumptions, for it makes the
  probabilities of staying in the same state decreasing with the
  distance to the next gene, and reaches a probability of 1/number of
  hidden state when the distance to the next gene is the maximum. \newline
  To avoid overflow errors, RJaCGH normalizes distances to be between 0
  and 1.
}
\value{
  ~Describe the value returned
  A matrix with the transition probabilities for that distance.
}
\references{Oscar Rueda and Ramon Diaz-Uriarte, in prep.}
\author{Oscar Rueda and Ramon Diaz-Uriarte}

\seealso{plot.Q.NH}
\examples{
## Model with two hidden states
## Note that RJaCGH normalizes distances to be between 0 and 1
beta <- matrix(c(0, 1, 3, 0), 2, 2)
Q.NH(beta=beta, x=0.4)
}
\keyword{models}% at least one, from doc/KEYWORDS

