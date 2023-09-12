data {
  int<lower=0> J;
  array[J] int y;
  array[J] int n;
}
parameters {
  real<lower=0,upper=1> phi;
  real<lower=1> kappa;
  vector<lower=0,upper=1>[J] theta;
}
model {
  kappa ~ pareto(1, 1.5);
  theta ~ beta(phi * kappa, (1-phi) * kappa);
  y ~ binomial(n, theta);
}
