\name{blinding.lslogit}
\alias{blinding.lslogit}
\alias{FitLogitStan}

\title{Latent Shift Logistic Regression}
\usage{
  blinding.lslogit(x, group, guess, mu0 = 0, s0 = 1, ...)
}
\arguments{
\item{x,guess}{outcome variable and guess response from blinding survey}
\item{group}{group assignments. Current version support one or two
  groups only}
\item{mu0,s0}{initial mean and sd of the latent variable of having sham effects}
\item{...}{Parameters ("adapt_delta","stepsize","max_treedepth") to improve model
  fitting/convergence.}
}
\value{
\itemize{
\item `sig.sham`: sd of effect size of sham treatment.
}
}
\description{
To be updated.
}
\examples{
\donttest{
guess <- rbinom(100,1,0.5)
group <- rbinom(100,1,0.5)
x <- rnorm(100,3,1)
blinding.lslogit(x=x,group=group,guess=guess)
}
}
