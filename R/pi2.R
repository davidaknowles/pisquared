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
#' @param upper_limit Upper limit of truncation (lower limit is fixed to 0)
#' @param alpha_prior Shape1 and shape2 of beta() prior on alpha. Shape1,shape2 > 1 prevents alpha collapsing to 0 or 1.
#' @param beta_prior Shape and rate of prior on beta_minus_one. Shape>1 enforces beta > 1.
#' @param method Can be "optimizing", "vb" or "sampling", determines approach for fitting stan model.
#' @param ... passed to rstan fitting function choosen by `method`
#' @import rstan
#' @return Fitted alpha and beta.
#'
fit_truncated_beta = function(p,
                              upper_limit,
                              alpha_prior = c(2,2),
                              beta_prior = c(2,1),
                              method = "optimizing",
                              ...) {
  dat = list(N = length(p),
             sum_logp = sum(log(p)),
             sum_log1mp = sum(log(1-p)),
             upper_limit = upper_limit,
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
    fit_func = if (method=="vb") vb else sampling
    tbeta_fit = fit_func( stanmodels$truncated_beta,
                          data = dat,
                          ...)
    samples = extract(tbeta_fit, pars = c("alpha","beta"))
    return(list(alpha = mean(samples$alpha),
                beta = mean(samples$beta)))
  }
}

#' MAP fit of pi0 assuming non-null come from a beta distribution.
#'
#' @export
#' @param p Samples (p-values for pisquared's application)
#' @param alpha_prior Shape1 and shape2 of beta() prior on alpha. Shape1,shape2 > 1 prevents alpha collapsing to 0 or 1.
#' @param beta_prior Shape and rate of prior on beta_minus_one. Shape>1 enforces beta > 1.
#' @param ... passed to rstan fitting function choosen by `method`
#' @import rstan
#' @return Fitted model parameters. .
model_based_pi0=function(p,
                         breaks=seq(0,1,by=0.01),
                         alpha_prior = c(2,2),
                         beta_prior = c(2,1),
                         ...) {
  h = hist(p, breaks=breaks, plot=F)
  dat = list(N=length(h$density),
             weights=h$density,
             pvalues=h$mids,
             alpha_prior = alpha_prior,
             beta_prior = beta_prior)
  pi0_fit = optimizing( pisquared:::stanmodels$pi0, data = dat, as_vector = F, ...)
  pi0_fit
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
#' @param breaks Grid for estimating p-value histogram
#' @param restarts Number of random restarts to perform.
#' @param alpha_prior Shape1 and shape2 of beta() prior on alpha. Shape1,shape2 > 1 prevents alpha collapsing to 0 or 1.
#' @param beta_prior Shape and rate of prior on beta_minus_one. Shape>1 enforces beta > 1.
#' @param ... passed to rstan fitting function choosen by `method`
#' @importFrom foreach "%do%" "%dopar%"
#' @import rstan
#' @return List of jaccard index and fitted model parameters.
#'
pi2_estimator=function(p1,
                       p2,
                       breaks=seq(0,1,by=0.01),
                       alpha_prior = c(2,2),
                       beta_prior = c(2,1),
                       restarts=5,
                       ...) {
  pvalues = list(p1, p2)

  pvalues = foreach(pv=pvalues) %do% { pmin(pv, 1) }

  init = foreach(pv=pvalues) %do% {
    #pi0_init = my_pi0est(pv)$pi0 # use qvalue version instead?
    #if (pi0_init >= 1.) warning(paste("qvalue pi0 estimate is",pi0_init,"don't trust results\n"))
    #pi0_init = pmin(pmax(pi0_init, 1e-8), 1 - 1e-8)
    fits=foreach(i=seq_len(restarts)) %dopar% {
      model_based_pi0(pv, alpha_prior = alpha_prior, beta_prior = beta_prior) }
    fit_likelihoods = foreach(fit=fits, .combine = c) %do% { fit$value }
    fits[[which.max(fit_likelihoods)]]$par
  }
  pi0_fit = c(init[[1]]$pi0, init[[2]]$pi0)
  pi0_init=c(pi0_fit[1]*pi0_fit[2],
             (1-pi0_fit[1])*pi0_fit[2],
             pi0_fit[1]*(1-pi0_fit[2]),
             (1-pi0_fit[1])*(1-pi0_fit[2]))

  N=length(breaks)-1
  h=gplots::hist2d(pvalues[[1]], pvalues[[2]], N, show = F)

  alpha=foreach(ini=init, .combine=c) %do% { ini$alpha }
  beta=foreach(ini=init, .combine=c) %do% { ini$beta }

  data=list(N=h$N,
            pvalues=h$x,
            weights=h$counts/sum(h$counts),
            alpha_prior = alpha_prior,
            beta_prior = beta_prior)
  fits=foreach(i=seq_len(restarts)) %dopar% {
    rstan::optimizing(pisquared:::stanmodels$pi2,
               data = data,
               init=list(pi0=pi0_init,
                         alpha=alpha,
                         beta=beta),
               as_vector=F,
               ...)
  }
  fit_likelihoods = foreach(fit=fits) %do% { fit$value }
  o=fits[[ which.max(fit_likelihoods ) ]]

  #o$par$pi0[4] - (o$par$pi0[2]+o$par$pi0[4])*(o$par$pi0[3]+o$par$pi0[4]) # covariance
  list(jaccard = o$par$pi0[4] / sum(o$par$pi0[2:4]),  # jaccard like index
       init = init,
       fit = o$par)
}
