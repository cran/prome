data {
  int<lower=0> J;
  int y[J];
  int n[J];
  real<lower=1> kappa;
}
parameters {
  real<lower=0,upper=1> phi;
  vector<lower=0,upper=1>[J] theta;
}
model {
  theta ~ beta(phi * kappa, (1-phi) * kappa);
  y ~ binomial(n, theta);
}
