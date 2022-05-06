data {
  int<lower=1> r0c; 
  int<lower=1> c0c; 
  matrix[r0c,c0c] w0c;   //x[i] is a vector for the ith row
  int<lower=1> r1c; 
  int<lower=1> c1c; 
  matrix[r1c,c1c] w1c; 
  int<lower=1> r0t; 
  int<lower=1> c0t; 
  matrix[r0t,c0t] w0t;
  int<lower=1> r1t; 
  int<lower=1> c1t; 
  matrix[r1t,c1t] w1t;
  real mu;
  real<lower=0> sigma;
}
parameters {
  real<lower=-1, upper=1> p0;
  real<lower=-1, upper=1> p1;
  row_vector[c0c] muc; 
  row_vector[c0t] mut; 
  real mu0;              // baseline mean
  real sigma0;           // baseline SD
  real<lower=0> sa;  // VAR of active ES
  real<lower=0> ss;    // VAR of placebo ES
  real<lower=0> tau;   //VAR of ME+within subj var 
}
transformed parameters {
  real<lower=0> sigma_t0;
  real<lower=0> sigma_t1;
  sigma_t0 = sqrt(tau+ss); 
  sigma_t1 = sqrt(tau+sa);
}
model {
  muc ~ normal(mu0, sqrt(sigma0));
  mut ~ normal(mu0, sqrt(sigma0));
  mu0 ~ normal(mu, sigma);
  sigma0 ~ inv_gamma(1.0, 1.0);
  tau ~ inv_gamma(1.0, 1.0);
  sa ~ inv_gamma(1.0, 1.0);
  ss ~ inv_gamma(1.0, 1.0);
  for(i in 1:r0c)
    w0c[i] ~ normal(muc, sqrt(tau));
  for(i in 1:r0t)
    w0t[i] ~ normal(mut, sqrt(tau));
  for(i in 1:r1c)
      w1c[i] ~ normal(muc*(1.0 + p0),sigma_t0);
  for(i in 1:r1t)
      w1t[i] ~ normal(mut*(1.0 + p0 + p1),sigma_t1);
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
