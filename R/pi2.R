#' Copy of qvalue's pi0est default. Included here to avoid dependency on bioconductor.
#'
#' Assumes pi0.method="smoother" and smooth.log.pi0 = F (see qvalue docs).
#'
#' @param lambda The value of the tuning parameter to estimate pi_0. Must be in [0,1). Optional, see Storey (2002).
#' @param smooth.df	Degrees-of-freedom to use when estimating pi_0 with a smoother. Optional.
#' @export
my_pi0est = function (p,
          lambda = seq(0.05, 0.95, 0.05),
          smooth.df = 3)
{
  if (length(p)==0) return(list(pi0=NA, pi0.lambda = NA*lambda, lambda = lambda))
  p <- p[!is.na(p)]
  lambda <- sort(lambda)
  ll <- length(lambda)
  pi0.lambda <- sapply(lambda, function(l) mean(p >= l)/(1 - l))
  spi0 <- smooth.spline(lambda, pi0.lambda, df = smooth.df)
  pi0Smooth <- predict(spi0, x = lambda)$y
  pi0 <- pi0Smooth[length(pi0Smooth)]
  return(list(pi0 = pi0, pi0.lambda = pi0.lambda, lambda = lambda))
}

#' Fit a beta distribution truncated to [0,upperlimit]
#'
#' @export
#' @param p Samples (p-values for pisquared's application)
#' @param fdr FDR threshold to use to define a set of "significant" p-values to fit the truncated beta to.
#' @param minimum_sig Require at least this many "significant" p-values for fitting. Otherwise return reasonable default.
#' @param alpha_prior Shape1 and shape2 of beta() prior on alpha. Shape1,shape2 > 1 prevents alpha collapsing to 0 or 1.
#' @param beta_prior Shape and rate of prior on beta_minus_one. Shape>1 enforces beta > 1.
#' @param method Can be "optimizing", "vb" or "sampling", determines approach for fitting stan model.
#' @param beta_fixed Whether to fix beta = 1.
#' @param chains Only used if method=="sampling": sets the number of MCMC chains to run.
#' @param refresh Only used if method=="vb": sets how frequently to output fitting info to stdout.
#' @param ... passed to rstan fitting function choosen by `method`
#' @import rstan
#' @return Fitted alpha and beta.
#'
fit_truncated_beta = function(p,
                              fdr = 0.05,
                              minimum_sig = 3,
                              alpha_prior = c(2,2),
                              beta_prior = c(2,1),
                              method = "optimizing",
                              beta_fixed = F,
                              chains = 1,
                              refresh = 0,
                              ...) {
  q = p.adjust(p, method = "BH")
  p_sig = p[q < fdr]
  if (length(p_sig) < minimum_sig) {
    warning(paste("Warning: only",length(p_sig),"significant p-values at 5% FDR"))
    return(list(alpha = 0.2, beta = 2)) # reasonable defaults
  }
  p_nonsig = p[q >= fdr]
  upper_limit = .5 * (max(p_sig) + min(p_nonsig))
  dat = list(N = length(p_sig),
             sum_logp = sum(log(p_sig)),
             sum_log1mp = sum(log(1-p_sig)),
             upper_limit = upper_limit,
             beta_fixed = beta_fixed,
             alpha_prior = alpha_prior,
             beta_prior = beta_prior)
  if (method == "optimizing") {
    tbeta_fit = optimizing( stanmodels$truncated_beta,
                                   data = dat,
                                   as_vector = F,
                                   ...)
    return(list(alpha = tbeta_fit$par$alpha,
                beta = tbeta_fit$par$beta))
  } else {

    inputs = list(...)
    if (method == "sampling") {
      inputs$chains = chains
    }
    inputs$refresh = refresh # control verbosity (see rstan::stan)
    inputs$object = stanmodels$truncated_beta
    inputs$data = dat
    tbeta_fit = do.call( if (method=="vb") vb else sampling, inputs )

    samples = extract(tbeta_fit, pars = c("alpha","beta"))
    return(list(alpha = mean(samples$alpha),
                beta = mean(samples$beta)))
  }
}

#' MAP fit of pi0 assuming non-null come from a beta distribution.
#'
#' @export
#' @param p Samples (p-values for pisquared's application)
#' @param bin_pvalues Whether to bin pvalues rather than using raw p-values. Can be more scalable.
#' @param breaks Grid for estimating p-value histogram (only used if bin_pvalues is True).
#' @param alpha_prior Shape1 and shape2 of beta() prior on alpha. Shape1,shape2 > 1 prevents alpha collapsing to 0 or 1.
#' @param beta_prior Shape and rate of prior on beta_minus_one. Shape>1 enforces beta > 1.
#' @param beta_fixed Whether to fix beta = 1.
#' @param truncated_beta_init Whether to initalize the beta distribution parameters using a truncated beta on significant p-values.
#' @param storey_init Whether to initialize pi0 using Storey's qvalue method.
#' @param min_pi0_init Minimum allowed pi0 initialization from Storey's method.
#' @param max_pi0_init Maximum allowed pi0 initialization from Storey's method.
#' @param method Can be "optimizing", "vb" or "sampling", determines approach for fitting stan model.
#' @param draws Only used if method=="optimizing": draw this many samples using Laplace's method (uses the Hessian of the loglikelihood at the optimum). These are used to get approximate standard errors on each parameter.
#' @param chains Only used if method=="sampling": sets the number of MCMC chains to run.
#' @param refresh Only used if method=="vb": sets how frequently to output fitting info to stdout.
#' @param ... passed to rstan fitting function choosen by `method`
#' @import rstan
#' @return Fitted model parameters. .
model_based_pi0=function(p,
                         bin_pvalues = F,
                         breaks=seq(0,1,by=0.01),
                         alpha_prior = c(2,2),
                         beta_prior = c(2,1),
                         beta_fixed = F,
                         truncated_beta_init = T,
                         storey_init = T,
                         min_pi0_init = 0.001,
                         max_pi0_init = .999,
                         method = "optimizing",
                         draws = 1000,
                         chains = 1,
                         refresh = 0,
                         ...) {
  init = list()
  if (truncated_beta_init) {
    tb_fit = fit_truncated_beta(p,
                               alpha_prior = alpha_prior,
                               beta_prior = beta_prior,
                               beta_fixed = beta_fixed,
                               method = method,
                               refresh = refresh)
    init$alpha = tb_fit$alpha
    if (!beta_fixed) init$beta_minus_one = array(tb_fit$beta - 1)
  }
  if (storey_init) {
    init$pi0 = my_pi0est(p)$pi0
    if (init$pi0 > max_pi0_init | init$pi0 < min_pi0_init)
      warning(paste("Warning: Storey's pi0 is",format(init$pi0,digits=5),"\n"))
    init$pi0 = pmin(pmax(init$pi0,min_pi0_init),max_pi0_init)
  }
  if (!bin_pvalues) {
    dat = list(N = length(p),
               weights = rep(1, length(p)),
               pvalues = p)
  } else {
    h = hist(p, breaks=breaks, plot=F)
    dat = list(N=length(h$density),
               weights=h$counts, # normalize? i.e. use density
               pvalues=h$mids)
  }
  dat$alpha_prior = alpha_prior
  dat$beta_prior = beta_prior
  dat$beta_fixed = as.integer(beta_fixed)
  if (method == "optimizing") {
    fit = optimizing( stanmodels$pi0,
                      data = dat,
                      init = init,
                      as_vector = F,
                      draws = draws,
                      ...)
    return(c(fit$par,
             alpha_sd = sd(fit$theta_tilde[,"alpha"]),
             beta_sd = sd(fit$theta_tilde[,"beta"]),
             pi0_sd = sd(fit$theta_tilde[,"pi0"])))
  } else {
    inputs = list(...)
    if (method == "sampling") {
      init = foreach(i=seq_len(chains)) %do% { init } # replicate
      inputs$chains = chains
    }
    inputs$refresh = refresh
    inputs$object = stanmodels$pi0
    inputs$data = dat
    inputs$init = init
    fit = do.call( if (method=="vb") vb else sampling, inputs )
    samples = extract(fit, pars = c("alpha","beta","pi0"))
    return(list(alpha = mean(samples$alpha),
                alpha_sd = sd(samples$alpha),
                beta = if(beta_fixed) 1.0 else mean(samples$beta),
                beta_sd = if(beta_fixed) 0.0 else sd(samples$beta),
                pi0 = mean(samples$pi0),
                pi0_sd = sd(samples$pi0)))
  }
}

#' Estimate sharing between two lists of p-values.
#'
#' Note the result is symmetric in p1 and p2. This will not work well if for either list all the p-values are uniform
#' (i.e. from the null) or all the p-values are significant. Also note this method (like Storey's pi0) assumes
#'
#' To run restarts in parallel registerDoMC needs to have been run beforehand.
#'
#' @export
#' @param p1 First list of p-values
#' @param p2 Second list of p-values
#' @param bin_pvalues Whether to bin pvalues rather than using raw p-values. Can be more scalable.
#' @param breaks Grid for estimating p-value histogram. Set to NULL for no binning (no pvalues can be 0!)
#' @param alpha_prior Shape1 and shape2 of beta() prior on alpha. Shape1,shape2 > 1 prevents alpha collapsing to 0 or 1.
#' @param beta_prior Shape and rate of prior on beta_minus_one. Shape>1 enforces beta > 1.
#' @param beta_fixed Whether to fix beta = 1.
#' @param truncated_beta_init Whether to initalize the beta distribution parameters using a truncated beta on significant p-values.
#' @param storey_init Whether to initialize pi0 using Storey's qvalue method.
#' @param method Can be "optimizing", "vb" or "sampling", determines approach for fitting stan model.
#' @param draws Only used if method=="optimizing": draw this many samples using Laplace's method (uses the Hessian of the loglikelihood at the optimum). These are used to get approximate standard errors on each parameter.
#' @param chains Only used if method=="sampling": sets the number of MCMC chains to run.
#' @param refresh Only used if method=="vb": sets how frequently to output fitting info to stdout.
#' @param ... passed to rstan fitting function choosen by `method`
#' @importFrom foreach "%do%" "%dopar%"
#' @import rstan
#' @return List of jaccard index, fitted model parameters, the stan fit object and the initiailization used.
pi2_estimator=function(p1,
                       p2,
                       bin_pvalues = F,
                       breaks=seq(0,1,by=0.01),
                       alpha_prior = c(2,2),
                       beta_prior = c(2,1),
                       beta_fixed = F,
                       truncated_beta_init = T,
                       storey_init = T,
                       method = "optimizing",
                       draws = 1000,
                       chains = 1,
                       refresh = 0,
                       ...) {
  pvalues = list(p1, p2)

  init = foreach(pv=pvalues) %do% {
    model_based_pi0(pv,
                    alpha_prior = alpha_prior,
                    beta_prior = beta_prior,
                    beta_fixed = beta_fixed,
                    bin_pvalues = bin_pvalues,
                    breaks = breaks,
                    truncated_beta_init = truncated_beta_init,
                    storey_init = storey_init,
                    method = method,
                    chains = chains,
                    refresh = refresh,
                    ...)
  }

  pi0_fit = c(init[[1]]$pi0, init[[2]]$pi0)
  pi_init = as.numeric(outer(c(pi0_fit[1], 1-pi0_fit[1]),
                              c(pi0_fit[2], 1-pi0_fit[2])))

  alpha=foreach(ini=init, .combine=c) %do% { ini$alpha }
  beta=foreach(ini=init, .combine=c) %do% { ini$beta }

  if (!bin_pvalues) {
    sm = stanmodels$pi2
    dat = list(N = length(p1),
               pvalues = rbind(p1,p2))
  } else {
    sm = stanmodels$pi2_hist
    h = gplots::hist2d(pvalues[[1]], pvalues[[2]], length(breaks), show = F)
    dat = list(N=length(breaks),
               pvalues=h$x,
               weights=h$counts) # normalize?
  }
  dat$alpha_prior = alpha_prior
  dat$beta_prior = beta_prior
  dat$beta_fixed = as.integer(beta_fixed)
  init=list(pi_array = pi_init,
            alpha = alpha)
  if (!beta_fixed) init$beta_minus_one = beta - 1
  if (method == "optimizing") {
    fit = optimizing(sm,
                    data = dat,
                    init = init,
                    draws = draws,
                    as_vector = F,
                    ...)
    # MAP estimate of jaccard is better than the pseudo-posterior mean for Jaccard index close
    # to 0 or 1. I suspect the Laplace approximation is poor in this regime.
    map_jacc = fit$par$pi_array[4] / sum(fit$par$pi_array[2:4])
    # these "samples" are generated by Stan using a Laplace approximation
    samples = list(pi_array = fit$theta_tilde[,paste0("pi_array[",1:4,"]"),drop=F],
                   alpha = fit$theta_tilde[,paste0("alpha[",1:2,"]"),drop=F],
                   beta = fit$theta_tilde[,paste0("beta[",1:2,"]"),drop=F] )
  } else {
    inputs = list(...)
    if (method == "sampling") {
      init = foreach(i=seq_len(chains)) %do% { init } # replicate
      inputs$chains = chains
    }
    inputs$refresh = refresh
    inputs$object = sm
    inputs$data = dat
    inputs$init = init
    fit = do.call( if (method=="vb") vb else sampling, inputs )
    samples = extract(fit, pars = c("alpha","beta","pi_array"))
  }
  jaccard = samples$pi_array[,4] / rowSums(samples$pi_array[,2:4,drop=F])
  list(fit = fit,
       init = init,
      jaccard = if (method=="optimizing") map_jacc else mean(jaccard),
      jaccard_sd = sd(jaccard),
      alpha = colMeans(samples$alpha),
      alpha_sd = apply(samples$alpha, 2, sd),
      beta = colMeans(samples$beta),
      beta_sd = apply(samples$beta, 2, sd),
      pi = matrix( if (method=="optimizing") fit$par$pi_array else colMeans(samples$pi_array),2),
      pi_sd = matrix(apply(samples$pi_array, 2, sd),2))
}


#' Simulate p-values for just one p-value array under the pisquared model.
#'
#' @export
#' @param N Sample size.
#' @param pi1 Marginal pi1 (proportion of non-null p-values).
#' @param alpha First parameter of beta distribution on non-null p-values. Should be <1.
#' @param beta Second parameter of beta distribution on non-null p-values. Should be >=1.
#' @return array of p-values
simulate_pvalues = function(N, pi1, alpha, beta) {
  ifelse( runif(N) > pi1, runif(N), rbeta(N, alpha, beta) )
}


#' Simulate p-values under the pisquared model with specified marginal pi0 and sharing rates.
#'
#' @export
#' @param N Sample size.
#' @param pi_1 Marginal pi1 for "discovery" data.
#' @param pi_2 Marginal pi1 for "replication" data.
#' @param true_jaccard Sharing level to aim for.
#' @param alpha First parameter of beta distribution on non-null p-values. Should be <1.
#' @param beta Second parameter of beta distribution on non-null p-values. Should be >=1.
#' @return Data.frame with which elements are hits in each cohort and corresponding p values.
simulate_joint_pvalues = function(N, pi_1, pi_2, true_jaccard, alpha, beta) {
  pi1_2_given_1not = ( pi_2 - true_jaccard * pi_1) / ((true_jaccard + 1) * (1 - pi_1))
  pi1_2_given_1hit = (pi_2 - (1 - pi_1) * pi1_2_given_1not) / pi_1
  if (pi1_2_given_1not < 0 | pi1_2_given_1not > 1 | pi1_2_given_1hit < 0 | pi1_2_given_1hit > 1) {
    warning("Warning: not possible to obtain the requested marginal pi1 and Jaccard index. Returning NULL.")
    return("NULL")
  }
  hit1 = runif(N) < pi_1
  p1 = ifelse( hit1,  rbeta(N, alpha, beta),  runif(N) )
  hit2 = ifelse(hit1, runif(N) < pi1_2_given_1hit, runif(N) < pi1_2_given_1not)
  p2 = ifelse( hit2,  rbeta(N, alpha, beta),  runif(N) )
  data.frame(hit1 = hit1, hit2 = hit2, p1 = p1, p2 = p2)
}
