data {
  int<lower=0> N; // number of grid points
  int<lower=0> K; // number of beta() mixture components for non-null
  vector[N] pvalues; // pvalue histogram centers
  vector[N] weights; // proportion of data in each bin
  real alpha_prior[2];
  real beta_prior[2];
  vector[K] mixture_weights_prior;
}
parameters {
  real<lower=0,upper=1> pi0;
  simplex[K] mixture_weights;
  real<lower=0, upper=1> alpha[K]; // constrain <1...
  real<lower=0> beta_minus_one[K]; // constrain >1 to model non-null p-values
}
transformed parameters {
  real beta[K];
  for (k in 1:K) beta[k] = beta_minus_one[k] + 1.0;
}
model {
  alpha ~ beta(alpha_prior[1], alpha_prior[2]);
  if (beta_prior[2] > 0) beta_minus_one ~ gamma(beta_prior[1], beta_prior[2]);
  mixture_weights ~ dirichlet(mixture_weights_prior);
  for(i in 1:N) {
    real mix_ll[K];
    for (k in 1:K) {
      mix_ll[k] = log(mixture_weights[k]) + beta_lpdf(pvalues[i] | alpha[k], beta[k]);
    }
    target += weights[i] * log_sum_exp(log(pi0) + beta_lpdf(pvalues[i] | 1, 1), log1m(pi0) + log_sum_exp(mix_ll));
  }
}
