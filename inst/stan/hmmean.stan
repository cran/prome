data {
  int<lower=0> N;
  int<lower=0> M;
  matrix[M,N] x;//x[i]= {x[i,] as in R}
  real<lower=0> sigma;
}
parameters {
  real theta;
  vector[M] mu;
  vector<lower=0>[M] s;
}
model {
  s ~ inv_gamma(1.0, 1.0);
  mu ~ normal(theta, sigma);
  for(i in 1:M){
    x[i] ~ normal(mu[i],sqrt(s[i]));
  }
}
