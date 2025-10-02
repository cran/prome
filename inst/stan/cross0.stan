data {
  int<lower=0> n1;         // number of subjects in group1
  int<lower=0> n2;         // number of subjects in group2
  vector[n1] y10;        // group1, baseline
  vector[n1] y11;        // group1, period 1
  vector[n1] y12;        // group1, period 2
  vector[n2] y20;        // group2, baseline
  vector[n2] y21;        // group2, period 1
  vector[n2] y22;        // group2, period 2
}
parameters {
  vector<lower=0,upper=1>[n1] mu1;
  vector<lower=0,upper=1>[n2] mu2;
  real mu;
  real<lower=0> sig;
  real mu11;
  real<lower=0> s11;
  real mu12;
  real<lower=0> s12;
  real mu21;
  real<lower=0> s21;
  real mu22;
  real<lower=0> s22;
}
//transformed parameters {
//  real<lower=0> a0;
//  real<lower=0> b0;
//  a0 = square(mu/sig)*(1-mu) - mu;
//  b0 = mu*(1-mu)/square(sig) - 1 - a0;
//}
model {
  y10 ~ normal(mu, sig);
  y20 ~ normal(mu, sig);
  s11 ~ inv_gamma(1.0, 1.0);
  s12 ~ inv_gamma(1.0, 1.0);
  s21 ~ inv_gamma(1.0, 1.0);
  s22 ~ inv_gamma(1.0, 1.0);
  //mu1 ~ beta(a0, b0);
  //mu2 ~ beta(a0, b0);
  mu1 ~ normal(mu, sig);
  mu2 ~ normal(mu, sig);
  y11 ~ normal(mu1 + mu11, s11);
  y21 ~ normal(mu2 + mu21, s21);
  y12 ~ normal(mu1 + mu12, s21);
  y22 ~ normal(mu2 + mu22, s22);
}
generated quantities {
  real tau_d; 
  real lambda_d; 
  real pi_d;
  lambda_d = (mu21 + mu22) - (mu11 + mu12);
  tau_d = mu12 - mu11;
  pi_d = 0.5 * (mu12 + mu22 - mu21 - mu11);
}
