data {
  int<lower=0> N; // number of grid points on each axis
  vector[N] pvalues;  // pvalue grid positions
  matrix[N,N] weights; // proportion of data in each pairwise bin
  real alpha_prior[2];
  real beta_prior[2];
  int<lower = 0, upper = 1> beta_fixed; // whether to fix betas to 1
}
parameters {
  simplex[4] pi_array;
  real<lower=0, upper=1> alpha[2];
  real<lower=0> beta_minus_one[beta_fixed ? 0 : 2]; // constrain >1 to model non-null p-values
}
transformed parameters {
  real beta[2];
  vector[4] log_pi_array = log(pi_array);
  for (i in 1:2)
    beta[i] = beta_fixed ? 1.0 : (beta_minus_one[i] + 1.0);
}
model {
  alpha ~ beta(alpha_prior[1], alpha_prior[2]);
  if (beta_prior[2] > 0) beta_minus_one ~ gamma(beta_prior[1], beta_prior[2]);
  for(i in 1:N) {
    for(j in 1:N) {
      vector[4] likelihoods;
      real like_1 = beta_lpdf(pvalues[i] | alpha[1], beta[1]);
      real like_2 = beta_lpdf(pvalues[j] | alpha[2], beta[2]);
      likelihoods[1] = 0;
      likelihoods[2] = like_1;
      likelihoods[3] = like_2;
      likelihoods[4] = like_1 + like_2;
      target += weights[i,j] * log_sum_exp(likelihoods + log_pi_array);
    }
  }
}
