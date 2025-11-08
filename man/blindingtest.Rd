\name{blinding.test}
\alias{blinding.test}
\alias{FitLogitStan}

\title{Latent Shift Logistic Regression}
\usage{
  blinding.test(x, group, guess, mu0 = 0, s0 = 1,...)
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
u1      = 5.5 # trt
u2      = 2.0 # ctrl
theta   = 3.2 # sham
sigma2  = 2.5   # v(rij)
ntreat  = 500      
nsham   = 500

beta0 = 1.0
beta1 = 2.0
beta2 = 1.0 # no contamination

Tind  = c(rep(1, ntreat), rep(0,nsham))  #treatment group indicator
u1v   = rep(u1,ntreat)
u2v   = rep(u2,nsham)
uv    = c(u1v,u2v)
tauv  = uv - rep(u2, ntreat+nsham)
r = rnorm(ntreat + nsham, mean = 0, sd = sqrt(sigma2))
q = 1/(1 + exp(-(beta0 + beta1*Tind + beta2*(tauv+r))))
bernGen = function(qq){rbinom(1,1,qq)}
I = sapply(q,bernGen)
x = uv + theta*I + r   # fixed sham effect
## I have concerns about the error term(s). x.sham~N(theta,sigma.sham)?
sigma.sham = 1.5
r2 = rnorm(ntreat + nsham, mean = 0, sd = sqrt(sigma.sham))
x = (uv + r) + theta*I #+ r2   # fixed sham effect

out1 <- blinding.test(x=x,group=Tind,guess=I);
out1

##data(bd012)
##blinding.test(x=bd012$y, group=bd012$group,guess=bd012$guess)
##data(bd011)
##blinding.test(x=bd011$y, group=bd011$group,guess=bd011$guess)
##data(bd010)
##blinding.test(x=bd010$y, group=bd010$group,guess=bd010$guess)
}
}
