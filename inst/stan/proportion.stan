data {
  int<lower=0> J;
  array[J] int y;
  array[J] int n;
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
