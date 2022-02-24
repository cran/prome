data {
  int<lower=0> J;         // number of subgroups 
  int y[J];              // number of events
  int n[J]; 		  // sample sizes
  real<lower=0,upper=1> pwr[J];	// power prior - weights
}
parameters {
  real<lower=0,upper=2> alpha;  // prior parameter
  //real bta;		        // prior parameter
  vector<lower=0,upper=1>[J] theta;    //treatment effect
}
model {
  theta ~ beta(alpha, 2-alpha);       // prior log-density
  for(i in 1:J)
    target += binomial_lpmf(y[i] | n[i], theta[i])*pwr[i];      
    //y ~ binomial(n, theta); // log-likelihood
}
