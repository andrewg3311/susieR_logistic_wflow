susie_logistic = function(Y, X, L = 10, V = 1, prior_weights = NULL, tol = 1e-3, maxit = 1000, intercept = TRUE) {
  
  p = ncol(X)
  n = nrow(X)
  
  # place to store posterior info for each l = 1, ..., L
  post_alpha = matrix(NA, nrow = p, ncol = L)
  post_mu = matrix(NA, nrow = p, ncol = L)
  post_sigma = matrix(NA, nrow = p, ncol = L)
  post_info = list(post_alpha = post_alpha, post_mu = post_mu, post_sigma = post_sigma)
  
  beta_post_init = matrix(Inf, nrow = p, ncol = L) # initialize
  beta_post_init2 = beta_post_init
  beta_post = matrix(0, nrow = p, ncol = L)
  
  fixed = rep(0, n) # fixed portion, estimated from l' != l other SER models
  
  iter = 0
  while((norm(beta_post - beta_post_init, "1") > tol) & (norm(beta_post - beta_post_init2, "1") > tol)) { # repeat until posterior means converge (ELBO not calculated here, so use this convergence criterion instead)
    beta_post_init2 = beta_post_init # store from 2 iterations ago
    beta_post_init = beta_post
    
    for (l in 1:L) {
      
      # below is old (inefficient) calculation of the fixed portion
      #fixed = rowSums(X %*% beta_post[, -l]) + int.coef # fixed portion from previous estimates (add intercept portion as well)
      fixed = fixed - (X %*% beta_post[, l]) # remove effect from previous iteration
      
      SER_logistric_l = single_effect_regression_logistic(Y, X, V, prior_weights, FALSE, fixed, intercept)
      # store
      post_info$post_alpha[, l] = SER_logistric_l$alpha
      post_info$post_mu[, l] = SER_logistric_l$mu
      post_info$post_sigma[, l] = SER_logistric_l$mu2 - SER_logistric_l$mu^2
      
      # update beta_post
      beta_post[, l] = SER_logistric_l$alpha * SER_logistric_l$mu
      
      fixed = fixed + (X %*% beta_post[, l]) # add back new fixed portion
      
    }
    
    iter = iter + 1
    if (iter > maxit) {
      stop("Maximum number of iterations reached")
    }
    
  }
  
  return(post_info)
  
}


# SER_logistic function
single_effect_regression_logistic = function(Y, X, V, prior_weights = NULL, optimize_V = FALSE, fixed = NULL, intercept = TRUE) {
  p = ncol(X)
  
  betahat = numeric(p)
  shat2 = numeric(p)
  
  if (is.null(fixed)) { # fixed is components from previous SER fits
    fixed = rep(0, length(Y))
  }
  
  # NOTE: could parallelize loop below if desired
  for (j in 1:p) { # logistic regression on each column of X separately
    if (intercept) {
      log.fit = glm(Y ~ X[, j] + 1 + offset(fixed), family = "binomial") # fit w/ intercept
    } else {
      log.fit = glm(Y ~ X[, j] - 1 + offset(fixed), family = "binomial") # fit w/out intercept
    }
    log.fit.coef = summary(log.fit)$coefficients
    # NOTE: coerces "intercept" to be 0 or 1 to grab relevant row of glm coefficient output
    betahat[j] = ifelse(is.na(log.fit.coef[1 + intercept, 1]), 0, log.fit.coef[1 + intercept, 1]) # beta-hat MLE (if na, just set to 0)
    shat2[j] = ifelse(is.na(log.fit.coef[1 + intercept, 2]), Inf, log.fit.coef[1 + intercept, 2]^2) # (std errof beta-hat MLE)^2 (if na, just set to Inf)
  }
  
  if (is.null(prior_weights)) {
    prior_weights = rep(1 / p, p)
  }
  
  if(optimize_V) {
    stop("Optimizing for prior variance not yet implemented for logistic case")
    #if(loglik.grad(0, betahat, shat2, prior_weights) < 0) {
    #  V = 0
    #} else {
    ##V.o = optim(par=log(V),fn=negloglik.logscale,gr = negloglik.grad.logscale,betahat=betahat,shat2=shat2,prior_weights=prior_weights,method="BFGS")
    ##if(V.o$convergence!=0){
    ##  warning("optimization over prior variance failed to converge")
    ##}
    #  V.u = uniroot(negloglik.grad.logscale, c(-10, 10), extendInt = "upX", betahat = betahat, shat2 = shat2, prior_weights = prior_weights)
    #  V = exp(V.u$root)
    #}
  }
  
  lbf = dnorm(betahat, 0, sqrt(V + shat2), log = TRUE) - dnorm(betahat, 0, sqrt(shat2), log = TRUE)
  #log(bf) on each SNP
  
  lbf[is.infinite(shat2)] = 0 # deal with special case of infinite shat2 (eg happens if X does not vary)
  
  maxlbf = max(lbf)
  w = exp(lbf - maxlbf) # w is proportional to BF, but subtract max for numerical stability
  # posterior prob on each SNP
  w_weighted = w * prior_weights
  weighted_sum_w = sum(w_weighted)
  alpha = w_weighted / weighted_sum_w
  post_var = 1 / ((1 / shat2) + (1 / V)) # posterior variance
  post_mean = (1 / shat2) * post_var * betahat # posterior mean
  post_mean2 = post_var + post_mean^2 # posterior second moment
  # BF for single effect model
  lbf_model = maxlbf + log(weighted_sum_w)
  # NOTE: Need to double check below
  loglik = lbf_model + log(1/2)*length(Y) # loglik of 0/1 response Y under p = .5
  return(list(alpha = alpha, mu = post_mean, mu2 = post_mean2, lbf = lbf, lbf_model = lbf_model, V = V, loglik = loglik))
}