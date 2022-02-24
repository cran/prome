data {
  int<lower=0> N;         // number of subjects
  vector[N] group;        // group: 0=control,1=treatment
  int<lower=0> r0;        // number of repeats@T0
  int<lower=0> r1;        // number of repeats@T1
  matrix[r0,N] w0;  
  matrix[r1,N] w1;  
}
parameters {
  vector[N] x0;		 // unknown true value
  real es0;              // effect size under sham treatment
  real es1;              // effect size under active treatment
  real mu0;              // prior location=baseline mean
  real sigma0;           // prior scale=baseline SD
  real mu1;              // prior location=baseline mean
  real sigma1;           // prior scale=baseline SD
  real<lower=0> tau;     // measurement noise
  real<lower=0> sig_active;  // SD of active ES
  real<lower=0> sig_sham;    // SD of placebo ES
}
transformed parameters {
  real<lower=0> sigma_t0;
  real<lower=0> sigma_t1;
  sigma_t0 = sqrt(tau*tau+sig_sham*sig_sham); 
  sigma_t1 = sqrt(tau*tau+sig_active*sig_active+sig_sham*sig_sham);
}
model {
  x0 ~ normal(mu0*(1-group)+mu1*group, sigma0*(1-group)+sigma1*group); 
  for(i in 1:r0)
    w0[i] ~ normal(x0, tau);    // measurement model
  for(i in 1:r1)
    w1[i] ~ normal(x0 + es0 + es1 * group,
    sigma_t0 * (1-group) + sigma_t1 * group);
}
generated quantities {
  real es_abs; 
  real es_rel; 
  es_rel = es1/mu1;
  es_abs = es1;
}
