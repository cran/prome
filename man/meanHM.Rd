\name{MeanHM}
\alias{MeanHM}
\title{Bayesian Hierarchical Model for Information Borrowing for Means}
\usage{
MeanHM(x,group,sigma)
}
\arguments{
\item{x}{Numeric vector of observations for the subgroups.}
\item{group}{group ID.}
\item{sigma}{hyper-parameter. to be estimated or can be given.}
}
\value{
\itemize{
\item `data`: data with estimates.
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
x <- c(x1,x2,x3)
group <- c(rep("grp1",length(x1)),
   rep("grp2",length(x2)),
   rep("grp3",length(x3)))
MeanHM(x,group,sigma=0.5)
}
