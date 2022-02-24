\name{PropHM}
\alias{PropHM}
\title{Bayesian Hierarchical Model for Information Borrowing for Proportions}
\usage{
PropHM(x, n, weights)
}
\arguments{
\item{x}{Numeric vector of events.}

\item{n}{Numberic vector of group sample sizes.}

\item{weights}{Numberic vector of weights (for the power prior method).}

}
\value{
\itemize{
\item `data`: data with estimates.
\item `alpha`: parameter of the beta distribution.
\item `beta`: parameter of the beta distribution.
}
}
\description{
To compute the proportions of the subgroups assuming the subgroups
follow the same binomial distribution with parameter \code{p}.
A beta prior with parameters \code{alpha} and \code{beta} are used
(\code{alpha+beta=2}).
}
\examples{
out <- PropHM(x=c(5,10,2),n=c(20,50,30))
}
