data {
  int<lower=0> N;
  int<lower=0> M;
  matrix[M,N] x;//x[i]= {x[i,] as in R}  
}
parameters {
  real theta;
  real<lower=0> sigma;
  vector[M] mu;
  vector<lower=0>[M] s2;
}
model {
  s2 ~ inv_gamma(1.0, 1.0);
  sigma ~ inv_gamma(1.0, 1.0);
  mu ~ normal(theta, sqrt(sigma));
  for(i in 1:M){
    x[i] ~ normal(mu[i],sqrt(s2[i]));
  }
}
