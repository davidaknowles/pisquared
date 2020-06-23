---
title: "pisquared-extended"
author: David A Knowles
output: 
  github_document:
    toc: true
    toc_depth: 2
    fig_width: 6
---




```r
library(pisquared)
require(tidyverse)
require(magrittr)
require(rstan)
require(doMC)
registerDoMC(detectCores() - 1) 
theme_set(theme_bw(base_size = 14))
```

## Model-based estimation of marginal pi0

Let's check how good a job we can do of estimating marginal pi0 across a range of parameters. We'll simulate some p-values with specified non-null rates. We'll use Beta(a=0.2,b=2) for the non-null pvalues: 

```r
true_alpha = 0.2
true_beta = 2
N = 1000
marg_results = foreach(true_pi1 = c(0.001,0.01,0.05,seq(.1,1,by=0.1)), .combine = bind_rows) %:%
  foreach(rep = 1:30, .combine = bind_rows) %dopar% {
    p = simulate_pvalues(N = N, pi1 = true_pi1, alpha = true_alpha, beta = true_beta)
    res = tryCatch(
      model_based_pi0(p, method = "optimizing") %>% 
        as_tibble() %>% 
        mutate(true_pi1 = true_pi1, 
               pi1 = 1 - pi0, 
               rep = rep, 
               qvalue_pi1 = 1 - pi0est_copy(p)$pi0),
      error = function(e) NULL
    )
}
marg_results
```

```
## # A tibble: 390 x 11
##      pi0 alpha beta_minus_one  beta alpha_sd beta_sd pi0_sd true_pi1      pi1   rep qvalue_pi1
##    <dbl> <dbl>          <dbl> <dbl>    <dbl>   <dbl>  <dbl>    <dbl>    <dbl> <int>      <dbl>
##  1 1.00  0.500          1.00   2.00    0.258    2.24 0.500     0.001 7.87e-11     1     0.0503
##  2 0.994 0.502          1.01   2.01    0.221    1.96 0.159     0.001 5.82e- 3     2    -0.106 
##  3 0.984 0.603          1.07   2.07    0.218    1.86 0.131     0.001 1.58e- 2     3    -0.0303
##  4 0.995 0.488          1.03   2.03    0.230    2.71 0.157     0.001 5.49e- 3     4    -0.244 
##  5 0.993 0.576          0.992  1.99    0.293    1.83 0.342     0.001 6.81e- 3     5     0.101 
##  6 0.990 0.408          0.995  1.99    0.191    2.00 0.0356    0.001 1.00e- 2     6    -0.0221
##  7 0.973 0.488          1.02   2.02    0.161    1.85 0.0405    0.001 2.65e- 2     7    -0.0593
##  8 0.955 0.845          0.919  1.92    0.230    1.15 0.196     0.001 4.51e- 2     8     0.101 
##  9 1.00  0.500          1.00   2.00    0.262    1.97 0.500     0.001 1.50e-10     9    -0.289 
## 10 0.994 0.303          0.984  1.98    0.159    2.18 0.0272    0.001 6.27e- 3    10     0.0833
## # … with 380 more rows
```

Both Storey's qvalue and our model based approach give approximately unbiased estimates of pi0 in this setting: 

```r
res_gathered = marg_results %>% select(true_pi1, model_based = pi1, qvalue = qvalue_pi1) %>% 
  gather(method, est_pi1, -true_pi1) 
res_gathered  %>% ggplot(aes(true_pi1, est_pi1, col = method, shape = method)) + geom_point(alpha = 0.5, size = 2, position = position_dodge(width = 0.02)) + 
  geom_abline(intercept = 0, slope = 1) + xlab("True pi1") + ylab("Estimated pi1") + theme(legend.position = c(0.25,.8))
```

![plot of chunk vary-true-pi0](figs/vary-true-pi0-1.png)

Errors are substantially lower for pisquared than Storey's pi0, as expected since the data is generated from the pisquared model. Additionally Storey's pi0 can be negative which is clearly unreasonable. 

```r
err_summary = res_gathered %>% 
  group_by(method, true_pi1) %>% 
  summarize(RMSE = sqrt(mean((true_pi1 - est_pi1)^2))) %>% 
  ungroup()
err_summary %>% ggplot(aes(true_pi1, RMSE, col = method, shape = method)) + geom_point(size = 2) + expand_limits(y = 0) + theme(legend.position = c(0.8,.8))
```

![plot of chunk rmse-comparison](figs/rmse-comparison-1.png)

alpha and beta are estimated quite robustly: 

```r
marg_results %>% 
  ggplot(aes(true_pi1, alpha)) + geom_point(alpha = 0.3, size = 3) + geom_hline(yintercept = true_alpha) + xlab("True pi1") + ylab("Estimated alpha")
```

![plot of chunk alpha-beta-estimation](figs/alpha-beta-estimation-1.png)

```r
marg_results %>% 
  ggplot(aes(true_pi1, beta)) + geom_point(alpha = 0.3, size = 3) + geom_hline(yintercept = true_beta) + xlab("True pi1") + ylab("Estimated beta")
```

![plot of chunk alpha-beta-estimation](figs/alpha-beta-estimation-2.png)

## Estimating sharing

We first consider independent p-values but with different non-null rates: 

```r
true_pi1s = c(0.1, 0.5)
p1 = simulate_pvalues(N = N, pi1 = true_pi1s[1], alpha = true_alpha, beta = true_beta)
p2 = simulate_pvalues(N = N, pi1 = true_pi1s[2], alpha = true_alpha, beta = true_beta)
hist(p1, 50)
```

![plot of chunk independent-p](figs/independent-p-1.png)

```r
hist(p2, 50)
```

![plot of chunk independent-p](figs/independent-p-2.png)

In our simulation p1 and p2 are independent so the "true" Jaccard index is small:

```r
true_jaccard = true_pi1s[1] * true_pi1s[2] / (1 - (1 - true_pi1s[1]) * (1 - true_pi1s[2]))
true_jaccard
```

```
## [1] 0.09090909
```

pisquared can use three methods provided by stan: maximum likelihood estimation "optimizing", variational Bayes ("vb") and Markov chain Monte Carlo ("sampling"). The fastest fitting approach is to do maximum likelihood estimation via LBFGS. 

```
## # A tibble: 6 x 7
##   method     bin_pvalues jaccard jaccard_sd pi0_1 pi0_2 cpu_time
##   <chr>      <lgl>         <dbl>      <dbl> <dbl> <dbl>    <dbl>
## 1 optimizing TRUE         0.0973     0.0371 0.869 0.622    0.747
## 2 optimizing FALSE        0.0825     0.0248 0.905 0.541    0.148
## 3 vb         TRUE         0.0826     0.0233 0.886 0.636   30.2  
## 4 vb         FALSE        0.0820     0.0221 0.902 0.537    3.51 
## 5 sampling   TRUE         0.108      0.0430 0.850 0.622  425.   
## 6 sampling   FALSE        0.0839     0.0253 0.903 0.536   39.1
```

The red line here shows the true Jaccard sharing index. binning pvalues does come at a cost in accuracy. 

```r
results %>% ggplot(aes(method, jaccard, 
                       ymin = jaccard - jaccard_sd, 
                       ymax = jaccard + jaccard_sd, 
                       fill = bin_pvalues)) + 
  geom_col(position = "dodge") + 
  geom_errorbar(position = "dodge") +
  geom_hline(yintercept = true_jaccard, 
             col = "red")
```

![plot of chunk compare-methods-accuracy](figs/compare-methods-accuracy-1.png)

The computation time for these methods is very different:

```r
results %>% mutate(method = factor(method, methods)) %>% 
  ggplot(aes(method, 1000 * cpu_time, fill = bin_pvalues)) + 
  geom_col(position = "dodge") + scale_y_log10() + ylab("Computation time (ms)")
```

![plot of chunk compare-methods-time](figs/compare-methods-time-1.png)

Optimization performs as well as VB or sampling and is an order of magnitude (or two) faster so we'll use it for the rest of the vignette:

```r
pi2_result = pi2_estimator(p1,
                           p2,
                           method = "optimizing")
marginal_pi1 = c(rowSums(pi2_result$pi)[2], colSums(pi2_result$pi)[2])
marginal_pi1
```

```
## [1] 0.09548413 0.45888194
```

Plot fitted marginal p-value distributions (estimated vs true)

```r
dpi0 = function(p, pi1, alpha, beta) {
  (1-pi1)  + pi1 * dbeta(p, alpha, beta)
}
for (i in 1:2) { 
  plot(function(g) dpi0(g, true_pi1s[i], true_alpha, true_beta), 
       ylab = "non-null p-val dist", 
       ylim = c(0,2),
       col = 1)
  abline(h = 1-true_pi1s[i], col = 1, lty = 2)
  plot(function(g) dpi0(g, marginal_pi1[i], pi2_result$fit$par$alpha[i], pi2_result$fit$par$beta[i]), 
       col = 2,
       add = T)
  abline(h = 1-marginal_pi1[i], col = 2, lty = 2)
  legend("topright", legend=c("True", "Estimated"), col=1:2, lty = 1)
  grid()
}
```

![plot of chunk check-marginal-distribution](figs/check-marginal-distribution-1.png)![plot of chunk check-marginal-distribution](figs/check-marginal-distribution-2.png)

## Estimating sharing across a range of underlying sharing levels

Here we simulation two arrays of p-values with specific sharing ("true_jaccard") and 

```r
N = 10000
jacc_results = foreach(true_jaccard = seq(0,1,by=.1), .combine = bind_rows) %:% 
  foreach(rep = 1:10, .combine = bind_rows) %dopar% {  # .combine = bind_rows, 
    set.seed(rep)
    dat = simulate_joint_pvalues(N = N, 
                                pi_1 = 0.1, 
                                pi_2 = 0.1,
                                true_jaccard = true_jaccard, 
                                alpha = true_alpha, 
                                beta = true_beta)
    if (is.null(dat)) return(NULL)
    counts = dat %>% select(hit1, hit2) %>% table()
    pi2_result <- pi2_estimator(dat$p1, dat$p2) # , method = "vb")
    tibble(true_jaccard = true_jaccard, 
           actual_jacc = counts["TRUE", "TRUE"] / (counts["FALSE", "TRUE"] +  counts["TRUE", "FALSE"] + counts["TRUE", "TRUE"]), 
           est_pi1_1 = rowSums(pi2_result$pi)[1], 
           est_pi1_2 = colSums(pi2_result$pi)[1],
           est_jaccard = pi2_result$jaccard, 
           jaccard_sd = pi2_result$jaccard_sd)
}
jacc_results
```

```
## # A tibble: 110 x 6
##    true_jaccard actual_jacc est_pi1_1 est_pi1_2 est_jaccard jaccard_sd
##           <dbl>       <dbl>     <dbl>     <dbl>       <dbl>      <dbl>
##  1            0           0     0.901     0.899    2.57e- 7     0.253 
##  2            0           0     0.895     0.897    2.71e- 3     0.0581
##  3            0           0     0.898     0.901    4.83e- 8     0.247 
##  4            0           0     0.896     0.904    1.13e- 7     0.240 
##  5            0           0     0.905     0.893    2.45e- 6     0.253 
##  6            0           0     0.905     0.901    1.13e- 5     0.235 
##  7            0           0     0.899     0.891    3.61e- 9     0.262 
##  8            0           0     0.909     0.905    1.50e- 7     0.256 
##  9            0           0     0.897     0.905    6.21e- 6     0.225 
## 10            0           0     0.899     0.902    8.69e-10     0.248 
## # … with 100 more rows
```

Plot the results: we see pisquared provides unbiased and accurate sharing estimates. 

```r
jacc_results %>% ggplot(aes(actual_jacc, est_jaccard)) + geom_point() + xlab("Actual sharing (Jaccard)") + ylab("Estimated sharing (Jaccard)") + geom_abline(intercept = 0, slope = 1)
```

![plot of chunk true-v-est-jaccard](figs/true-v-est-jaccard-1.png)

