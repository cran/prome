data {
  int<lower=0> N;
  real<lower=0> sigma;
  vector[N] x;   
  vector[N] sig;   
}
parameters {
  real theta;    # overall effect
  vector[N] mu;
}
model {
  mu ~ normal(theta, sigma);
  x ~ normal(mu,sig);
}
