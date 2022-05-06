data {
  int<lower=0> N;         // number of subjects
  vector[N] group;        // group: 0=control,1=treatment
  int<lower=0> r1;         // number of subjects
  int<lower=0> r0;         // number of subjects
  matrix[r0,N] w0;//x[i]= {x[i,] as in R}  
  matrix[r1,N] w1;//x[i]= {x[i,] as in R}  
  real mu;
  real<lower=0> sigma;
}
parameters {
  real<lower=-1, upper=1> p0;
  real<lower=-1, upper=1> p1;
  vector[N] mui;        // group: 0=control,1=treatment
  real mu0;              // prior location=baseline mean
  real<lower=0> sigma0;           // prior scale=baseline SD
  real<lower=0> sa;  // var of active ES
  real<lower=0> ss;    // var of placebo ES
  real<lower=0> tau;   //var of ME+within subj var 
}
transformed parameters {
  real<lower=0> sigma_t0;
  real<lower=0> sigma_t1;
  sigma_t0 = sqrt(tau+ss); 
  sigma_t1 = sqrt(tau+sa);
}
model {
  mui ~ normal(mu0, sqrt(sigma0));
  mu0 ~ normal(mu, sigma);
  sigma0 ~ inv_gamma(1.0, 1.0);
  sa ~ inv_gamma(1.0, 1.0);
  ss ~ inv_gamma(1.0, 1.0);
  tau ~ inv_gamma(1.0, 1.0);
  for(i in 1:r0)
    w0[i] ~ normal(mui, sqrt(tau));
  for(i in 1:r1){
      w1[i] ~ normal(mui .* (1.0 + p0 + p1*group),
        sigma_t0 * (1.0-group) + sigma_t1 * group);
  }
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
