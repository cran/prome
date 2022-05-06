data {
  int<lower=0> N;         // number of subjects
  vector[N] group;        // group: 0=control,1=treatment
  real<lower=0> r1;        // number of repeats@T1
  int<lower=0> rz;        // number of covariates
  vector[N] w0;  
  vector[N] w1;
  real<lower=0> tau;     // measurement noise
  matrix[rz,N] z;  
}
parameters {
  vector[rz] beta;	 // coefficients of covariates
  real es0;              // effect size under sham treatment
  real es1;              // effect size under active treatment
  real mu0;              // prior location=baseline mean
  real sigma0;           // prior scale=baseline SD
  real<lower=0> sig_active;  // SD of active ES
  real<lower=0> sig_sham;    // SD of placebo ES
}
transformed parameters {
  real<lower=0> sigma_t0;
  real<lower=0> sigma_t1;
  sigma_t0 = sqrt(tau*tau/r1+sig_sham*sig_sham); 
  sigma_t1 = sqrt(tau*tau/r1+sig_active*sig_active+sig_sham*sig_sham);
}
model {
  w0 ~ normal(mu0, sigma0);
  for(i in 1:N)
    w1[i] ~ normal(w0[i]*(es0 + es1*group[i]) + z[i] * beta,
            sigma_t0 * (1-group[i]) + sigma_t1 * group[i]);
}
generated quantities {
  real es_abs; 
  real es_rel; 
  es_rel = es1;
  es_abs = es1*mu0;
}
