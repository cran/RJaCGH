\name{plot.Q.NH}
\alias{plot.Q.NH}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot transition probabilities}
\description{
Plot the probabilities of staying in the same state for a non-homogenous
hidden Markov model.
}
\usage{
plot.Q.NH(x, beta, q=-beta, col = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Vector of distances between observations}
  \item{beta}{\code{beta} parameter of the transition matrix. Must be a square
    matrix with the same size as the number of hidden states}
  \item{q}{\code{q} parameter of the transition matrix. Must be a square
    matrix with the same size as the number of hidden states}
  \item{col}{vector of colors for each state. Must be of the same size
    as the number of hidden states}
  \item{\dots}{aditional arguments passed to plot}
}
\details{
Please note that RJaCGH model imposes that \code{q} is \code{-beta}, and
distances are normalized to lay between 0 and 1.
}
\value{
A plot is produced showing the probability of staying in the same hidden
state versus distance between adjacent genes, for every state.
}
\references{Oscar   Oscar M. Rueda and Ramon Diaz Uriarte. A flexible, accurate and
  extensible statistical method for detecting genomic copy-number
  changes. \url{http://biostats.bepress.com/cobra/ps/art9/}.
{http://biostats.bepress.com/cobra/ps/art9/}.
}

\author{Oscar M. Rueda and Ramon Diaz-Uriarte}
\examples{
## Model with two hidden states
## Note that RJaCGH normalizes distances to be between 0 and 1
x <- rexp(99)
x <- x/ max(x)
beta <- matrix(c(0, 1, 3, 0), 2, 2)
plot.Q.NH(x=x, beta=beta, q=-beta, col=c(1,2))
}
\keyword{models}% at least one, from doc/KEYWORDS

