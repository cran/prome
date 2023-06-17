data {
  int<lower=0> n1;         // number of subjects in group1
  int<lower=0> n2;         // number of subjects in group2
  vector[n1] y10;        // group1, baseline
  vector[n1] y11;        // group1, period 1
  vector[n1] y12;        // group1, period 2
  vector[n2] y20;        // group2, baseline
  vector[n2] y21;        // group2, period 1
  vector[n2] y22;        // group2, period 2
}
parameters {
  vector[n1] mu1;
  vector[n2] mu2;
  real mu;
  real<lower=0> sig;
  real tau1;
  real<lower=0> stau1;
  real tau2;
  real<lower=0> stau2;
  real pi_d;
  real<lower=0> spid;
  real lambda_d;
  real<lower=0> slmdd;
}
transformed parameters {
  real<lower=0> sig11;
  real<lower=0> sig12;
  real<lower=0> sig21;
  real<lower=0> sig22;
  sig11 = sqrt(sig*sig + stau1*stau1); 
  sig21 = sqrt(sig*sig + stau2*stau2); 
  sig12 = sqrt(sig*sig + spid*spid + stau2*stau2); 
  sig22 = sqrt(sig*sig + spid*spid + stau1*stau1 + slmdd*slmdd); 
}
model {
  y10 ~ normal(mu, sig);
  y20 ~ normal(mu, sig);
  mu1 ~ normal(mu, sig);
  mu2 ~ normal(mu, sig);
  stau1 ~ inv_gamma(1.0, 1.0);
  stau2 ~ inv_gamma(1.0, 1.0);
  spid ~ inv_gamma(1.0, 1.0);
  slmdd ~ inv_gamma(1.0, 1.0);
  y11 ~ normal(mu1 + tau1, sig11);
  y21 ~ normal(mu2 + tau2, sig21);
  y12 ~ normal(mu1 + tau2 + pi_d, sig12);
  y22 ~ normal(mu2 + tau1 + pi_d + lambda_d, sig22);
}
generated quantities {
  real tau_d; 
  tau_d = tau2 - tau1;
}
