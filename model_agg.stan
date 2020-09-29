data {
  int<lower=0> J; //number of patients
  int<lower=0> N; //number of observations
  real y[N]; //patient outcome
  real t[N]; //treatment choice
  int<lower=0> id[N]; //patient identifier {1,..,J}
}
parameters {
  vector[2] beta_mu; //population mean
  vector<lower=0>[2] beta_sigma; //population standard deviation between patients
  vector[J] eta1; //variable to be transformed into betas for non-centered parameterization
  vector[J] eta2;
  real<lower=0> y_sig; //standard deviation within patients
}
transformed parameters {
  vector[J] beta1; // individual regression coefficients
  vector[J] beta2;
  
  beta1 = beta_mu[1] + beta_sigma[1] * eta1;
  beta2 = beta_mu[2] + beta_sigma[2] * eta2;
  //computing the group-level coefficient, based on non-centered parameterization
}
model {
  //Priors
      y_sig ~ gamma(2,0.1);
      beta_mu ~ normal(0,50);
      beta_sigma ~ gamma(2,0.1);
      eta1 ~ normal(0,1);
      eta2 ~ normal(0,1);
  //Likelihood
  {
  vector[N] xbeta;
  for (n in 1:N)
    xbeta[n] = beta1[id[n]] + beta2[id[n]].*t[n];
  y ~ normal(xbeta, y_sig);
  }
}
