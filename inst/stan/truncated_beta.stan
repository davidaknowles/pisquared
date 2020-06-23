// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0> N;
  real sum_logp; // sum(log(p)) of significant p-values
  real sum_log1mp; // sum(log(1-p)) of significant p-values
  real upper_limit; // upper bound of considered pvalues
  real alpha_prior[2];
  real beta_prior[2];
  int<lower = 0, upper = 1> beta_fixed; // whether to fix betas to 1
}

parameters {
  real<lower=0, upper=1> alpha; // constrain <1...
  real<lower=0> beta_minus_one[beta_fixed ? 0 : 1]; // constrain >1 to model non-null p-values
}

transformed parameters {
  real beta = beta_fixed ? 1.0 : (beta_minus_one[1] + 1.0);
}

model {
  alpha ~ beta(alpha_prior[1], alpha_prior[2]);
  if (beta_prior[2] > 0) beta_minus_one ~ gamma(beta_prior[1], beta_prior[2]);
  target += (alpha - 1.) * sum_logp + (beta - 1.) * sum_log1mp;
  // to account for truncation:
  target += - N * (beta_lcdf(upper_limit | alpha, beta) + lbeta(alpha, beta));
}

