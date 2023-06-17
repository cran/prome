\name{xover}
\alias{xover}
\alias{print.xover}

\title{Bayesian analysis of 2x2 crossover trial data}
\usage{
  xover(group,y1,y2,y0,...)
}
\arguments{
\item{y0,y1,y2}{vectors of data from baseline, period 1, and period 2,
  respectively.}
\item{group}{group or treatment sequence.}
\item{...}{other parameters, i.e. 'control' for model fitting.}
}
\value{
\itemize{
\item `stat`: summary statistics.
\item `best`: estimates using Bayesian analysis.
}
}
\description{
A Bayesian hierachical model to analysis data from 2x2 (AB/BA) crossover trials.
}
\examples{
 \donttest{
   xover(y0=rnorm(20,34,1.5),y1=rnorm(20,30,2),
         y2=rnorm(20,25,1.5),group=round(runif(20)<0.5))
 }
}