data {
  int<lower=0> J; //number of patients
  int<lower=0> N; //number of observations
  real y[N]; //patient outcome
  real t[N]; //treatment choice
  int<lower=0> id[N]; //patient identifier {1,..,J}
}
parameters {
  vector[J] beta1; // individual regression coefficients
  vector[J] beta2;
  real<lower=0> y_sig; //standard deviation within patients
}
model {
  //Priors
      y_sig ~ gamma(2,0.1);
      beta1 ~ normal(0,50);
      beta2 ~ normal(0,50);
  //Likelihood
  {
  vector[N] xbeta;
  for (n in 1:N)
    xbeta[n] = beta1[id[n]] + beta2[id[n]].*t[n];
  y ~ normal(xbeta, y_sig);
  }
}
