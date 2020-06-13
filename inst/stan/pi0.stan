data {
  int<lower=0> N; // number of grid points
  vector[N] pvalues; // pvalue histogram centers
  vector[N] weights; // proportion of data in each bin
  real alpha_prior[2];
  real beta_prior[2];
}
parameters {
  real<lower=0,upper=1> pi0;
  real<lower=0, upper=1> alpha; // constrain <1...
  real<lower=0> beta_minus_one; // constrain >1 to model non-null p-values
}
transformed parameters {
  real beta = beta_minus_one + 1.0 ;
}
model {
  alpha ~ beta(alpha_prior[1], alpha_prior[2]);
  beta_minus_one ~ gamma(beta_prior[1], beta_prior[2]);
  for(i in 1:N) {
    target += weights[i] * log_sum_exp(log(pi0) + beta_lpdf(pvalues[i] | 1, 1), log1m(pi0) + beta_lpdf(pvalues[i] | alpha, beta));
  }
}
