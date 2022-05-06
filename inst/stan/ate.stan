data {
  int<lower=1> N; 
  vector[N] w; 
  vector[N] arm; 
  vector[N] tm;
  int<lower=1> nsubj;
  real mu;
  real sigma;
}
parameters {
  real<lower=-1, upper=1> p0;
  real<lower=-1, upper=1> p1;
  vector[nsubj] mui; 
  real mu0;              // baseline mean
  real sigma0;           // baseline SD
  real<lower=0> sa;  // var of active ES
  real<lower=0> ss;    // var of placebo ES
  real<lower=0> tau;   //var of ME+within subj var 
}
transformed parameters {
  vector[N] muii;
  vector[N] sii;
  for(i in 1:N){
    muii[i] = mui[1+(i-1)%nsubj]*(1.0+p0*tm[i]+p1*tm[i]*arm[i]);
    sii[i] = sqrt(tau + ss*tm[i] + sa*tm[i]*arm[i]);
  }
}
model {
  mui ~ normal(mu0, sqrt(sigma0));
  mu0 ~ normal(mu, sigma);
  tau ~ inv_gamma(1.0, 1.0);
  sigma0 ~ inv_gamma(1.0, 1.0);
  sa ~ inv_gamma(1.0, 1.0);
  ss ~ inv_gamma(1.0, 1.0);
  w ~ normal(muii,sii);
}
generated quantities {
  real es_abs; 
  real es_rel; 
  real es0; 
  real es1; 
  es_rel = p1;
  es_abs = p1 * mu0;
  es0 = p0 * mu0;
  es1 = (p1+p0) * mu0;
}
