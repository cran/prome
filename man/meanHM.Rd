\name{MeanHM}
\alias{MeanHM}
\title{Bayesian Hierarchical Model for Information Borrowing for Means}
\usage{
MeanHM(x,sigma)
}
\arguments{
\item{x}{Numeric vector of observations for the subgroups.}
\item{sigma}{hyper-parameter. to be estimated or can be given.}
}
\value{
\itemize{
\item `theta`: population mean.
\item `sigma`: population standard deviation.
}
}
\description{
To compute the mean values of subgroups based on a Bayesian hierarchical model.
}
\examples{
x1 <- rnorm(100,2,1)
x2 <- rnorm(100,3,1.5)
x3 <- rnorm(100,4,1.9)
x <- cbind(x1,x2,x3)
MeanHM(x,sigma=0.5)
}
