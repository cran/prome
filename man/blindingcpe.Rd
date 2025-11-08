\name{blinding.cpe}
\alias{cpe}
\alias{blinding.cpe}
\title{Changing point estimator to correct bias due to unblinding in RCTs}
\usage{
  blinding.cpe(x,group,guess)
}
\arguments{
\item{x}{outcome variable. numeric vector}
\item{group}{group assignment. Coded as "0"=control and
  "1"=active/treatment. If 'group' is a factor, the first level
  will be treated as "control" arm. For example, if there are two
  values (ie. "ctrl" and "active"), "active" will be treated as
  the control arm. If the two levels are "control" and "treatment",
  "control" will be treated as the control arm.}
\item{guess}{responses to the blinding survey question. The response
  corresponding to positive sham effect needs to be coded as "1", and
  the rest as "0". If the possible responses are "Active", "Control" and
  "I don't know", code "guess" as "1" for "Active" and "0" for the
  others if a subject responded "active" is expect to have positive sham
  effect. Otherwise, if a subject responded "active" is expect to have
  negative sham effect, code "guess" as "0" for "Active" and "1" for the
  others.}
}
\value{
3 estimates, BI indices.
}
\description{
 estimate the sham effect and correct effect size
}
\details{
 To be added
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

out1 <- blinding.cpe(x=x,group=Tind,guess=I);
out1

##data(bd012)
##blinding.cpe(x=bd012$y, group=bd012$group,guess=bd012$guess)
##data(bd011)
##blinding.cpe(x=bd011$y, group=bd011$group,guess=bd011$guess)
##data(bd010)
##blinding.cpe(x=bd010$y, group=bd010$group,guess=bd010$guess)
}

}
