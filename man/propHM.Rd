\name{PropHM}
\alias{PropHM}
\title{Bayesian Hierarchical Model for Information Borrowing for Proportions}
\usage{
PropHM(x, n, kappa)
}
\arguments{
\item{x}{Numeric vector of events.}
\item{n}{Numberic vector of group sample sizes.}
\item{kappa}{\code{kappa=alpha+beta>1}. Must be given if the number of
  subgroups is 2.}
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
The approach on partial pooling by Bob Carpenter has been used --
"Hierarchical Partial Pooling for Repeated Binary Trials"
https://mc-stan.org/users/documentation/case-studies/pool-binary-trials.html
}
\examples{
out <- PropHM(x=c(5,10,2),n=c(20,50,30))
}
