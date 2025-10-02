data {
  int<lower=1> N;
  int<lower=0,upper=1> guess[N];
  int<lower=0,upper=1> group[N];
  vector[N] y;
  vector[N] z;
  real mu0;                // prior mean for theta
  real<lower=0> s0;        // prior sd for theta
}

parameters {
  // separate intercepts
  real beta0;     // outcome baseline
  real beta1;     // shifted outcome
  real beta_group;
  real beta_y;

  // shift parameter
  real gamma0;     // latent class intercept
  real gamma1;     // latent class slope
  real theta;
}

model {
  // priors
  theta ~ normal(mu0, s0);

  beta0 ~ normal(0, 5);
  beta1 ~ normal(0, 5);
  beta_group ~ normal(0, 5);
  beta_y ~ normal(0, 5);
  gamma0 ~ normal(0, 5);
  gamma1 ~ normal(0, 5);

  for (i in 1:N) {
    // outcome logits
    real eta0 = beta0 + beta_group * group[i] + beta_y * y[i];
    real eta1 = beta1 + beta_group * group[i] + beta_y * (y[i] - theta);

    // latent class logit 
    real pz = inv_logit(gamma0 + beta_group * group[i] + gamma1 * z[i]);

    // mixture
    real pi = pz * inv_logit(eta1) + (1 - pz) * inv_logit(eta0);

    target += bernoulli_lpmf(guess[i] | pi);
  }
}
