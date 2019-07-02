### VB LOWER BOUND SUSIE BELOW

### NOTE: This is basic code, and I did not attempt to mirror the level of numerical sophistication in the susie functions
### If this idea is worth pursuing further, then this code can be improved

calc_Q = function(X, Sigma2, Mu, Alpha, Z, delta) {
  # X is data matrix
  # Sigma2[j, l] is the posterior variance for b_l when entry j is selected, p x L
  # Mu[j, l] is the posterior mean for b_l when entry j is selected, p x L
  # Alpha[j, l] is the posterior probability selecting entry j from b_l, p x L
  # Z is matrix of covariates (e.g. column for intercept, top 10 PCs, etc)
  # delta is current estimate for effects of Z variables
  
  ASU2 = Alpha * (Sigma2 + Mu^2) # [j, l] = alpha[j, l] * (Sigma2[j, l] + Mu[j, l]^2)
  
  Q = rowSums(X^2 %*% ASU2) # start w/ sum_l E_(ql)[(x_i' b_l)^2]
  
  # now add 2 sum_l sum_{k>l} (x_i' b_l_post)(x_i' b_k_post)
  # maybe try to think of a more efficient way of doing this?
  b_post_mat = Mu * Alpha # each column is posterior mean for b_l, p x L
  X_b_post = X %*% b_post_mat # [i, l] = x_i' b_l_post
  for (i in 1:nrow(X)) {
    xb_i = as.numeric(X_b_post[i, ])
    xb_outer = xb_i %o% xb_i # outer-product, [j, k] = (x_i'b_j_post)(x_i' b_k_post)
    diag(xb_outer) = 0 # remove diagonal, where j = k
    Q[i] = Q[i] + sum(xb_outer)
  }
  
  # now, add other terms with Z and delta
  Q = Q + as.numeric(2 * (X %*% rowSums(b_post_mat)) * (Z %*% delta)) + as.numeric((Z %*% delta)^2)
  
  return(Q)
  
}

g = function(x) { # ilogit function
  1 / (1 + exp(-x))
}


update_xi = function(X, Sigma2, Mu, Alpha, Z, delta) {
  # X is data matrix
  # Sigma2[j, l] is the posterior variance for b_l when entry j is selected, p x L
  # Mu[j, l] is the posterior mean for b_l when entry j is selected, p x L
  # Alpha[j, l] is the posterior probability selecting entry j from b_l, p x L
  # Z is matrix of covariates (e.g. column for intercept, top 10 PCs, etc)
  # delta is current estimate for effects of Z variables
  
  Q = calc_Q(X, Sigma2, Mu, Alpha, Z, delta)
  
  xi = sqrt(Q)
  
  return(xi)
  
}

# update_b_l_OLD = function(X, Y, xi, prior_weights, V, Sigma2, Mu, Alpha, l, Z, delta) {
#   # X is data matrix
#   # Y is binary response
#   # xi is lower-bound approximation parameters
#   # prior_weights is prior probabilities for selecting j (p-vector)
#   # V is prior variance
#   # Sigma2[j, l] is the posterior variance for b_l when entry j is selected, p x L
#   # Mu[j, l] is the posterior mean for b_l when entry j is selected, p x L
#   # Alpha[j, l] is the posterior probability selecting entry j from b_l, p x L
#   # l is index to update
#   # Z is matrix of covariates (e.g. column for intercept, top 10 PCs, etc)
#   # delta is current estimate for effects of Z variables
# 
#   b_post_mat = Mu * Alpha # each column is posterior mean for b_l, p x L
#   b_post_not_l = rowSums(as.matrix(b_post_mat[, -l], nrow = nrow(Mu))) # posterior, sum_(k != l) b_k_post
#   g_xi = g(xi) # vector of g(xi_i), pre-compute once
# 
# 
#   for (j in 1:ncol(X)) {
#     common_denom = (1 / V) + sum(as.numeric(X[, j]^2) * (g_xi - .5) / xi) # appears many times, compute once
# 
#     # update Alpha[j, l]
#     num = sum(X[, j] * (Y - .5 -  ((1 / xi) * (g_xi - .5) * ((X %*% b_post_not_l) + (Z %*% delta))))) # numerator in exp
#     denom = 2*common_denom # denominator in exp
#     Alpha[j, l] = log(prior_weights[j]) + (num^2 / denom) + (1/2)*log(1 / common_denom) # on log-scale, for stability
# 
#     # update Mu[j, l]
#     Mu[j, l] = num / common_denom
# 
#     # update Sigma[j, l]
#     Sigma2[j, l] = 1 / common_denom
#   }
#   Alpha[, l] = exp(Alpha[, l] - max(Alpha[, l])) # remove max for stability, everything still proportional
#   Alpha[, l] = Alpha[, l] / sum(Alpha[, l]) # normalize, sum to 1
# 
#   return(list(Sigma2 = Sigma2, Mu = Mu, Alpha = Alpha))
# 
# }


update_b_l = function(X, Y, xi, prior_weights, V, Sigma2, Mu, Alpha, l, Z, delta) {
  # X is data matrix
  # Y is binary response
  # xi is lower-bound approximation parameters
  # prior_weights is prior probabilities for selecting j (p-vector)
  # V is prior variance
  # Sigma2[j, l] is the posterior variance for b_l when entry j is selected, p x L
  # Mu[j, l] is the posterior mean for b_l when entry j is selected, p x L
  # Alpha[j, l] is the posterior probability selecting entry j from b_l, p x L
  # l is index to update
  # Z is matrix of covariates (e.g. column for intercept, top 10 PCs, etc)
  # delta is current estimate for effects of Z variables
  
  b_post_mat = Mu * Alpha # each column is posterior mean for b_l, p x L
  b_post_not_l = rowSums(as.matrix(b_post_mat[, -l], nrow = nrow(Mu))) # posterior, sum_(k != l) b_k_post
  g_xi = g(xi) # vector of g(xi_i), pre-compute once
  g_xi_5_xi = ((g_xi - .5) / xi) # vector of (g(xi_i) - .5) / xi_i, pre-compute once
  g_xi_5_xi[xi == 0] = .25 # case of 0/0 (e.g. x_i is all 0), use L'Hopital
  
  common_denoms = (1 / V) + (t(X^2) %*% g_xi_5_xi) # appears many times, compute once, posterior precisions
  
  # update Alpha[, l]
  nums = t(X) %*% (Y - .5 -  (g_xi_5_xi * ((X %*% b_post_not_l) + (Z %*% delta)))) # numerator in exp
  Alpha[, l] = as.numeric(log(prior_weights) + (nums^2 / (2*common_denoms)) - (1/2)*log(common_denoms))
  Alpha[, l] = exp(Alpha[, l] - max(Alpha[, l])) # remove max for stability, everything still proportional
  Alpha[, l] = Alpha[, l] / sum(Alpha[, l]) # normalize, sum to 1
  
  # update Mu[, l]
  Mu[, l] = as.numeric(nums / common_denoms)
  
  # update Sigma[, l]
  Sigma2[, l] = as.numeric(1 / common_denoms)
  
  return(list(Sigma2 = Sigma2, Mu = Mu, Alpha = Alpha))
  
}


update_delta = function(X, Y, xi, Mu, Alpha, Z) {
  # X is data matrix
  # Y is binary response
  # xi is lower-bound approximation parameters
  # Mu[j, l] is the posterior mean for b_l when entry j is selected, p x L
  # Alpha[j, l] is the posterior probability selecting entry j from b_l, p x L
  # Z is matrix of covariates (e.g. column for intercept, top 10 PCs, etc)
  
  g_xi = g(xi) # vector of g(xi_i), pre-compute once
  g_xi_5_xi = ((g_xi - .5) / xi) # vector of (g(xi_i) - .5) / xi_i, pre-compute once
  g_xi_5_xi[xi == 0] = .25
  
  #D = Matrix::Diagonal(x = g_xi_5_xi / 2)
  D = Matrix::Diagonal(x = g_xi_5_xi)
  ZtDZ = t(Z) %*% D %*% Z
  DXB = D %*% X %*% rowSums(Mu * Alpha)
  RHS = t(Z) %*% (Y - .5 - DXB) # RHS of system LL'delta = RHS
  
  # NOTE: the below system can/should be solved w/ Cholesky and forward/backward substitution (but if Z has small # of columns, shouldn't matter)
  #Lt = chol(ZtDZ)
  #u = forwardsolve(t(Lt), RHS)
  #delta = backsolve(Lt, u)
  ## COULD USE chol2inv!!!!
  #delta = solve(ZtDZ, RHS)
  delta = solve(as.matrix(ZtDZ), as.matrix(RHS))
  
  return(as.numeric(delta))
  
}

# update_int_VB = function(X, Y, xi, Mu, Alpha, prior_var) {
#   g_xi = g(xi) # vector of g(xi_i), pre-compute once
#   g_xi_5_xi = ((g_xi - .5) / xi) # vector of (g(xi_i) - .5) / xi_i, pre-compute once
#   g_xi_5_xi[xi == 0] = .25
#   d = g_xi_5_xi
#   
#   b_post = rowSums(Mu * Alpha)
#   
#   num = sum(Y - .5 - d*(X %*% b_post))
#   denom = (1 / prior_var) + sum(d)
#   
#   mu_int = num / denom
#   sigma2_int = 1 / denom
#   
#   return(list(mu_int = mu_int, sigma2_int = sigma2_int))
#   
# }

# update_b_l_VB_int = function(X, Y, xi, prior_weights, V, Sigma2, Mu, Alpha, l, mu_int) {
#   # X is data matrix
#   # Y is binary response
#   # xi is lower-bound approximation parameters
#   # prior_weights is prior probabilities for selecting j (p-vector)
#   # V is prior variance
#   # Sigma2[j, l] is the posterior variance for b_l when entry j is selected, p x L
#   # Mu[j, l] is the posterior mean for b_l when entry j is selected, p x L
#   # Alpha[j, l] is the posterior probability selecting entry j from b_l, p x L
#   # l is index to update
#   # Z is matrix of covariates (e.g. column for intercept, top 10 PCs, etc)
#   # mu_int is current mean of intercept distribution
#   
#   b_post_mat = Mu * Alpha # each column is posterior mean for b_l, p x L
#   b_post_not_l = rowSums(as.matrix(b_post_mat[, -l], nrow = nrow(Mu))) # posterior, sum_(k != l) b_k_post
#   g_xi = g(xi) # vector of g(xi_i), pre-compute once
#   g_xi_5_xi = ((g_xi - .5) / xi) # vector of (g(xi_i) - .5) / xi_i, pre-compute once
#   g_xi_5_xi[xi == 0] = .25 # case of 0/0 (e.g. x_i is all 0), use L'Hopital
#   
#   common_denoms = (1 / V) + (t(X^2) %*% g_xi_5_xi) # appears many times, compute once, posterior precisions
#   
#   # update Alpha[, l]
#   nums = t(X) %*% (Y - .5 -  (g_xi_5_xi * ((X %*% b_post_not_l) + mu_int))) # numerator in exp
#   Alpha[, l] = as.numeric(log(prior_weights) + (nums^2 / (2*common_denoms)) - (1/2)*log(common_denoms))
#   Alpha[, l] = exp(Alpha[, l] - max(Alpha[, l])) # remove max for stability, everything still proportional
#   Alpha[, l] = Alpha[, l] / sum(Alpha[, l]) # normalize, sum to 1
#   
#   # update Mu[, l]
#   Mu[, l] = as.numeric(nums / common_denoms)
#   
#   # update Sigma[, l]
#   Sigma2[, l] = as.numeric(1 / common_denoms)
#   
#   return(list(Sigma2 = Sigma2, Mu = Mu, Alpha = Alpha))
#   
# }

# calc_Q_VB_int = function(X, Sigma2, Mu, Alpha, mu_int, sigma2_int) {
#   # X is data matrix
#   # Sigma2[j, l] is the posterior variance for b_l when entry j is selected, p x L
#   # Mu[j, l] is the posterior mean for b_l when entry j is selected, p x L
#   # Alpha[j, l] is the posterior probability selecting entry j from b_l, p x L
#   # Z is matrix of covariates (e.g. column for intercept, top 10 PCs, etc)
#   # delta is current estimate for effects of Z variables
#   
#   ASU2 = Alpha * (Sigma2 + Mu^2) # [j, l] = alpha[j, l] * (Sigma2[j, l] + Mu[j, l]^2)
#   
#   Q = rowSums(X^2 %*% ASU2) # start w/ sum_l E_(ql)[(x_i' b_l)^2]
#   
#   # now add 2 sum_l sum_{k>l} (x_i' b_l_post)(x_i' b_k_post)
#   # maybe try to think of a more efficient way of doing this?
#   b_post_mat = Mu * Alpha # each column is posterior mean for b_l, p x L
#   X_b_post = X %*% b_post_mat # [i, l] = x_i' b_l_post
#   for (i in 1:nrow(X)) {
#     xb_i = as.numeric(X_b_post[i, ])
#     xb_outer = xb_i %o% xb_i # outer-product, [j, k] = (x_i'b_j_post)(x_i' b_k_post)
#     diag(xb_outer) = 0 # remove diagonal, where j = k
#     Q[i] = Q[i] + sum(xb_outer)
#   }
#   
#   # now, add other terms with Z and delta
#   Q = Q + as.numeric(2 * (X %*% rowSums(b_post_mat)) * mu_int) + (sigma2_int + mu_int^2)
#   
#   return(Q)
#   
# }


# update_xi_VB_int = function(X, Sigma2, Mu, Alpha, mu_int, sigma2_int) {
#   # X is data matrix
#   # Sigma2[j, l] is the posterior variance for b_l when entry j is selected, p x L
#   # Mu[j, l] is the posterior mean for b_l when entry j is selected, p x L
#   # Alpha[j, l] is the posterior probability selecting entry j from b_l, p x L
#   # Z is matrix of covariates (e.g. column for intercept, top 10 PCs, etc)
#   # delta is current estimate for effects of Z variables
#   
#   Q = calc_Q_VB_int(X, Sigma2, Mu, Alpha, mu_int, sigma2_int)
#   
#   xi = sqrt(Q)
#   
#   return(xi)
#   
# }


update_all = function(X, Y, xi, prior_weights, V, Sigma2, Mu, Alpha, Z, delta, estimate_prior_variance) {
  # X is data matrix
  # Y is binary response
  # xi is lower-bound approximation parameters
  # prior_weights is prior probabilities for selecting j (p-vector)
  # V is prior variance
  # Sigma2[j, l] is the posterior variance for b_l when entry j is selected, p x L
  # Mu[j, l] is the posterior mean for b_l when entry j is selected, p x L
  # Alpha[j, l] is the posterior probability selecting entry j from b_l, p x L
  # Z is matrix of covariates (e.g. column for intercept, top 10 PCs, etc)
  # delta is current estimate for effects of Z variables
  # estimate_prior_variance is logical for if prior variance, V, should be estimated
  
  # update delta
  if (any(Z != 0)) { # if covariates and/or intercept
    delta = update_delta(X, Y, xi, Mu, Alpha, Z)
  }
  
  # now, iterate over l = 1:L
  for (l in 1:ncol(Mu)) {
    res_l = update_b_l(X, Y, xi, prior_weights, V, Sigma2, Mu, Alpha, l, Z, delta)
    Sigma2 = res_l$Sigma2
    Mu = res_l$Mu
    Alpha = res_l$Alpha
  }
  
  # now, update xi
  xi = update_xi(X, Sigma2, Mu, Alpha, Z, delta)
  
  if (estimate_prior_variance == TRUE) {
    ASU2 = Alpha * (Sigma2 + Mu^2) # [j, l] = alpha[j, l] * (Sigma2[j, l] + Mu[j, l]^2)
    V = sum(ASU2) / ncol(Mu) # ncol(Mu) = L = sum_l sum_j alpha_jl
  }

  
  return(list(Sigma2 = Sigma2, Mu = Mu, Alpha = Alpha, delta = delta, xi = xi, V = V))
  
}

# This function updates delta and xi simultaneously
# update_all_alternate = function(X, Y, xi, prior_weights, V, Sigma2, Mu, Alpha, Z, delta, tol = 1e-6) {
#   # X is data matrix
#   # Y is binary response
#   # xi is lower-bound approximation parameters
#   # prior_weights is prior probabilities for selecting j (p-vector)
#   # V is prior variance
#   # Sigma2[j, l] is the posterior variance for b_l when entry j is selected, p x L
#   # Mu[j, l] is the posterior mean for b_l when entry j is selected, p x L
#   # Alpha[j, l] is the posterior probability selecting entry j from b_l, p x L
#   # Z is matrix of covariates (e.g. column for intercept, top 10 PCs, etc)
#   # delta is current estimate for effects of Z variables
#   # tol is tolerance for updating delta and xi
#   
#   # iterate over l = 1:L
#   for (l in 1:ncol(Mu)) {
#     res_l = update_b_l(X, Y, xi, prior_weights, V, Sigma2, Mu, Alpha, l, Z, delta)
#     Sigma2 = res_l$Sigma2
#     Mu = res_l$Mu
#     Alpha = res_l$Alpha
#   }
#   
#   # now, coordinate ascet of delta and xi (alternate until converge)
#   # first, pre-compute expensive part of Q (doesn't change)
#   ASU2 = Alpha * (Sigma2 + Mu^2) # [j, l] = alpha[j, l] * (Sigma2[j, l] + Mu[j, l]^2)
#   
#   Q_fixed = rowSums(X^2 %*% ASU2) # start w/ sum_l E_(ql)[(x_i' b_l)^2]
#   
#   # now add 2 sum_l sum_{k>l} (x_i' b_l_post)(x_i' b_k_post)
#   # maybe try to think of a more efficient way of doing this?
#   b_post_mat = Mu * Alpha # each column is posterior mean for b_l, p x L
#   X_b_post = X %*% b_post_mat # [i, l] = x_i' b_l_post
#   for (i in 1:nrow(X)) {
#     xb_i = as.numeric(X_b_post[i, ])
#     xb_outer = xb_i %o% xb_i # outer-product, [j, k] = (x_i'b_j_post)(x_i' b_k_post)
#     diag(xb_outer) = 0 # remove diagonal, where j = k
#     Q_fixed[i] = Q_fixed[i] + sum(xb_outer)
#   }
#   
#   
#   xi_old = rep(Inf, length(xi))
#   while (sqrt(sum((xi - xi_old)^2)) > tol) {
#     xi_old = xi
#     delta = update_delta(X, Y, xi, Mu, Alpha, Z)
#     Q = Q_fixed + as.numeric(2 * (X %*% rowSums(b_post_mat)) * (Z %*% delta)) + as.numeric((Z %*% delta)^2)
#     xi = sqrt(Q)
#   }
#   
#   return(list(Sigma2 = Sigma2, Mu = Mu, Alpha = Alpha, delta = delta, xi = xi))
#   
# }

# update_all_VB_int = function(X, Y, xi, prior_weights, V, Sigma2, Mu, Alpha, int_prior_var) {
#   # X is data matrix
#   # Y is binary response
#   # xi is lower-bound approximation parameters
#   # prior_weights is prior probabilities for selecting j (p-vector)
#   # V is prior variance
#   # Sigma2[j, l] is the posterior variance for b_l when entry j is selected, p x L
#   # Mu[j, l] is the posterior mean for b_l when entry j is selected, p x L
#   # Alpha[j, l] is the posterior probability selecting entry j from b_l, p x L
#   # Z is matrix of covariates (e.g. column for intercept, top 10 PCs, etc)
#   # mu_int is current mean estimate for intercept
#   # sigma2_int is current variance estimate for intercept
#   
#   # update delta
#   int = update_int_VB(X, Y, xi, Mu, Alpha, int_prior_var)
#   
#   # now, iterate over l = 1:L
#   for (l in 1:ncol(Mu)) {
#     res_l = update_b_l_VB_int(X, Y, xi, prior_weights, V, Sigma2, Mu, Alpha, l, int$mu_int)
#     Sigma2 = res_l$Sigma2
#     Mu = res_l$Mu
#     Alpha = res_l$Alpha
#   }
#   
#   # now, update xi
#   xi = update_xi_VB_int(X, Sigma2, Mu, Alpha, int$mu_int, int$sigma2_int)
#   
#   
#   return(list(Sigma2 = Sigma2, Mu = Mu, Alpha = Alpha, int = int, xi = xi))
#   
# }

# Y is vector of binary response (n x 1, can be sparse)
# X is matrix of variables (n x p, can be sparse)
# L is number of non-zero effects, positive integer
# V is prior variance on non-zero effect size, positive real number
# prior_weights is a vector of prior inclusion probabilities (p x 1)
# intercept is a logical of if the intercept should be fitted (default to TRUE)
# Z is a vector of covariates to be controlled for (non-penalized effect estimates, n x q, can be sparse).
# NOTE: Z should NOT include an intercept column
# estimate_prior_variance is a logical if prior variance, V, should be estimated. If "TRUE", supplied value of V is the initial value
# tol is the convergence criterior measuring the change in the ELBO
# maxit is the maximum number of iterations
susie_logistic_VB = function(Y, X, L = 10, V = 1, prior_weights = NULL, intercept = TRUE, Z = NULL, estimate_prior_variance = FALSE, tol = 1e-3, maxit = 1000) {
  p = ncol(X)
  n = nrow(X)
  
  if (is.null(Z)) {
    if (intercept == TRUE) {
      Z = matrix(1, nrow = n, ncol = 1)
    } else {
      Z = matrix(0, nrow = n, ncol = 1) # if no intercept and no control covariates, set to 0
    }
  } else {
    col_variances = apply(Z, MARGIN = 2, var)
    if (any(col_variances == 0)) { # is constant column in Z matrix
      stop("Matrix 'Z' cannot have a constant column")
    }
    Z = cbind(matrix(1, nrow = n, ncol = 1), Z) # add intercept column
  }
  
  if (is.null(prior_weights)) {
    prior_weights = rep(1 / p, p)
  }
  
  # place to store posterior info for each l = 1, ..., L
  # initialize: could think of something better
  #delta = glm(Y ~ Z - 1, family = "binomial")$coef # initialize to regular GLM solution
  if (all(Z == 0)) { # if no covariates and no intercept
    delta = 0
  } else {
    delta = glm(as.numeric(Y) ~ Z - 1, family = "binomial")$coef # initialize to regular GLM solution w/ just Z (no X)
  }
  Alpha = matrix(prior_weights, nrow = p, ncol = L)
  #Alpha = t(MCMCpack::rdirichlet(L, prior_weights)) # alternate initialization method
  Mu = matrix(0, nrow = p, ncol = L)
  Sigma2 = matrix(V, nrow = p, ncol = L)
  xi = update_xi(X, Sigma2, Mu, Alpha, Z, delta)
  post_info = list(Sigma2 = Sigma2, Mu = Mu, Alpha = Alpha, delta = delta, xi = xi, V = V)
  
  #beta_post_init = matrix(Inf, nrow = p, ncol = L) # initialize
  #beta_post_init2 = beta_post_init
  beta_post = post_info$Alpha * post_info$Mu
  
  ELBOs = numeric(maxit + 1)
  ELBOs[1] = -Inf
  ELBOs[2] = calc_ELBO(Y, X, post_info$Alpha, post_info$Mu, post_info$Sigma2, post_info$V, prior_weights, Z, post_info$delta, post_info$xi)
  
  iter = 1
  #while((norm(beta_post - beta_post_init, "1") > tol) & (norm(beta_post - beta_post_init2, "1") > tol)) { # repeat until posterior means
  while(abs(ELBOs[iter+1] - ELBOs[iter]) > tol) { # repeat until ELBO increase is negligible
    post_info = update_all(X, Y, post_info$xi, prior_weights, post_info$V, post_info$Sigma2, post_info$Mu, post_info$Alpha, Z, post_info$delta, estimate_prior_variance)
    
    beta_post = post_info$Alpha * post_info$Mu
    
    iter = iter + 1
    ELBOs[iter + 1] = calc_ELBO(Y, X, post_info$Alpha, post_info$Mu, post_info$Sigma2, post_info$V, prior_weights, Z, post_info$delta, post_info$xi)
    if (iter > maxit) {
      stop("Maximum number of iterations reached")
    }
  }
  
  # change output format to match the GLM version's output
  int = NULL
  delta = post_info$delta
  if (intercept == TRUE) {
    int = post_info$delta[1]
    delta = post_info$delta[-1]
  }
  if (length(delta) == 0) {
    delta = 0
  }
  res = list(post_alpha = post_info$Alpha, post_mu = post_info$Mu, post_sigma = post_info$Sigma2, intercept = int, delta = delta, xi = post_info$xi, ELBOs = ELBOs[2:(iter+1)], prior_variance = post_info$V)
  
  return(res)
  
}


# susie_logistic_VB_alternate = function(Y, X, L = 10, V = 1, prior_weights = NULL, Z = 1, tol = 1e-3, maxit = 1000) {
#   p = ncol(X)
#   n = nrow(X)
#   
#   if (Z == 1) { # run with intercept
#     Z = matrix(1, nrow = n, ncol = 1)
#   }
#   
#   if (is.null(prior_weights)) {
#     prior_weights = rep(1 / p, p)
#   }
#   
#   # place to store posterior info for each l = 1, ..., L
#   # initialize: could think of something better
#   #delta = glm(Y ~ Z - 1, family = "binomial")$coef # initialize to regular GLM solution
#   delta = glm(as.numeric(Y) ~ Z - 1, family = "binomial")$coef # initialize to regular GLM solution w/ just Z (no X)
#   Alpha = matrix(prior_weights, nrow = p, ncol = L)
#   Mu = matrix(0, nrow = p, ncol = L)
#   Sigma2 = matrix(V, nrow = p, ncol = L)
#   xi = update_xi(X, Sigma2, Mu, Alpha, Z, delta)
#   post_info = list(Sigma2 = Sigma2, Mu = Mu, Alpha = Alpha, delta = delta, xi = xi)
#   
#   #beta_post_init = matrix(Inf, nrow = p, ncol = L) # initialize
#   #beta_post_init2 = beta_post_init
#   beta_post = post_info$Alpha * post_info$Mu
#   
#   ELBOs = numeric(maxit + 1)
#   ELBOs[1] = -Inf
#   ELBOs[2] = calc_ELBO(Y, X, post_info$Alpha, post_info$Mu, post_info$Sigma2, V, prior_weights, Z, post_info$delta, post_info$xi)
#   
#   iter = 1
#   #while((norm(beta_post - beta_post_init, "1") > tol) & (norm(beta_post - beta_post_init2, "1") > tol)) { # repeat until posterior means
#   while(abs(ELBOs[iter+1] - ELBOs[iter]) > tol) { # repeat until ELBO increase is negligible
#     post_info = update_all_alternate(X, Y, post_info$xi, prior_weights, V, post_info$Sigma2, post_info$Mu, post_info$Alpha, Z, post_info$delta, 1e-6)
#     
#     beta_post = post_info$Alpha * post_info$Mu
#     
#     iter = iter + 1
#     ELBOs[iter + 1] = calc_ELBO(Y, X, post_info$Alpha, post_info$Mu, post_info$Sigma2, V, prior_weights, Z, post_info$delta, post_info$xi)
#     if (iter > maxit) {
#       stop("Maximum number of iterations reached")
#     }
#   }
#   
#   # change output format to match the GLM version's output
#   res = list(post_alpha = post_info$Alpha, post_mu = post_info$Mu, post_sigma = post_info$Sigma2, delta = post_info$delta, xi = post_info$xi, ELBOs = ELBOs[2:(iter+1)])
#   
#   return(res)
#   
# }
# 
# susie_logistic_VB_int = function(Y, X, L = 10, V = 1, prior_weights = NULL, int_prior_var = 100, tol = 1e-3, maxit = 1000) {
#   p = ncol(X)
#   n = nrow(X)
#   
#   if (is.null(prior_weights)) {
#     prior_weights = rep(1 / p, p)
#   }
#   
#   # place to store posterior info for each l = 1, ..., L
#   # initialize: could think of something better
#   glm.init = glm(as.numeric(Y) ~ 1, family = "binomial")
#   int = list(mu_int = summary(glm.init)$coef[1], sigma2_int = 1 / ((1 / int_prior_var) + (1 / summary(glm.init)$coef[2])))
#   
#   Alpha = matrix(prior_weights, nrow = p, ncol = L)
#   Mu = matrix(0, nrow = p, ncol = L)
#   Sigma2 = matrix(V, nrow = p, ncol = L)
#   xi = update_xi_VB_int(X, Sigma2, Mu, Alpha, int$mu_int, int$sigma2_int)
#   post_info = list(Sigma2 = Sigma2, Mu = Mu, Alpha = Alpha, int = int, xi = xi)
#   
#   #beta_post_init = matrix(Inf, nrow = p, ncol = L) # initialize
#   #beta_post_init2 = beta_post_init
#   beta_post = post_info$Alpha * post_info$Mu
#   
#   ELBOs = numeric(maxit + 1)
#   ELBOs[1] = -Inf
#   ELBOs[2] = calc_ELBO_VB_int(Y, X, post_info$Alpha, post_info$Mu, post_info$Sigma2, V, prior_weights, post_info$int, post_info$xi, int_prior_var)
#   
#   iter = 1
#   #while((norm(beta_post - beta_post_init, "1") > tol) & (norm(beta_post - beta_post_init2, "1") > tol)) { # repeat until posterior means
#   while(abs(ELBOs[iter+1] - ELBOs[iter]) > tol) { # repeat until ELBO increase is negligible
#     post_info = update_all_VB_int(X, Y, post_info$xi, prior_weights, V, post_info$Sigma2, post_info$Mu, post_info$Alpha, int_prior_var)
#     
#     beta_post = post_info$Alpha * post_info$Mu
#     
#     iter = iter + 1
#     ELBOs[iter + 1] = calc_ELBO_VB_int(Y, X, post_info$Alpha, post_info$Mu, post_info$Sigma2, V, prior_weights, post_info$int, post_info$xi, int_prior_var)
#     if (iter > maxit) {
#       stop("Maximum number of iterations reached")
#     }
#   }
#   
#   # change output format to match the GLM version's output
#   res = list(post_alpha = post_info$Alpha, post_mu = post_info$Mu, post_sigma = post_info$Sigma2, int = post_info$int, xi = post_info$xi, ELBOs = ELBOs[2:(iter+1)])
#   
#   return(res)
#   
# }


# calculate the variational lower bound
# CAREFUL: Not sure what to do when Alpha[j, l] = 0 (we get 0*ln(0)). I will set this to 0
# Note: This expression is only valid when xi has been updated to be sqrt(Q), where Q is the 
# expectation of the square of the linear predictor under the variational posterior (what we update xi to nomrally)
calc_ELBO = function(Y, X, Alpha, Mu, Sigma2, V, prior_weights, Z, delta, xi) {
  p = nrow(Mu)
  L = ncol(Mu)
  P = matrix(prior_weights, nrow = p, ncol = L)
  b_post = rowSums(Alpha * Mu)
  
  expected_log_lik = sum(log(g(xi)) + (Y - .5) * as.numeric(((X %*% b_post) + (Z %*% delta))) - (xi / 2))
  KL_div_vb_prior = Alpha * (log(Alpha) - log(P) + (log(V) / 2) - (log(Sigma2) / 2) - .5 + ((Sigma2 + Mu^2) / (2 * V)))
  KL_div_vb_prior[Alpha == 0] = 0
  KL_div_vb_prior = sum(KL_div_vb_prior)
  
  if (KL_div_vb_prior < 0) { # to diagnose any issues
    warning("KL Divergence calculated to be < 0")
  }
  
  # TERMS BELOW WRONG B/C FORGOT MULTIPLICIATIVE FACTOR OF alpha_jl
  #KL_div_vb_prior = sum(log(Alpha)) - L*sum(log(prior_weights)) + (L*p*log(V) / 2) - sum(log(Sigma2) / 2) - (L*p / 2) + (sum(Sigma2 + Mu^2) / (2*V))
  #KL_div_vb_prior = sum(log(Alpha)) - sum(log(Sigma2) / 2) + (sum(Sigma2 + Mu^2) / (2*V)) # up to a constant
  
  ELBO = expected_log_lik - KL_div_vb_prior
  return(ELBO)
}



# calc_ELBO_VB_int = function(Y, X, Alpha, Mu, Sigma2, V, prior_weights, int, xi, int_prior_var) {
#   p = nrow(Mu)
#   L = ncol(Mu)
#   P = matrix(prior_weights, nrow = p, ncol = L)
#   b_post = rowSums(Alpha * Mu)
#   
#   expected_log_lik = sum(log(g(xi)) + (Y - .5) * as.numeric(((X %*% b_post) + int$mu_int)) - (xi / 2))
#   KL_div_vb_prior = Alpha * (log(Alpha) - log(P) + (log(V) / 2) - (log(Sigma2) / 2) - .5 + ((Sigma2 + Mu^2) / (2 * V)))
#   KL_div_vb_prior[Alpha == 0] = 0
#   KL_div_vb_prior = sum(KL_div_vb_prior) - (log(int$sigma2_int) / 2) + ((int$sigma2_int + int$mu_int^2) / (2*int_prior_var))
#   
#   if (KL_div_vb_prior < 0) { # to diagnose any issues
#     warning("KL Divergence calculated to be < 0")
#   }
#   
#   # TERMS BELOW WRONG B/C FORGOT MULTIPLICIATIVE FACTOR OF alpha_jl
#   #KL_div_vb_prior = sum(log(Alpha)) - L*sum(log(prior_weights)) + (L*p*log(V) / 2) - sum(log(Sigma2) / 2) - (L*p / 2) + (sum(Sigma2 + Mu^2) / (2*V))
#   #KL_div_vb_prior = sum(log(Alpha)) - sum(log(Sigma2) / 2) + (sum(Sigma2 + Mu^2) / (2*V)) # up to a constant
#   
#   ELBO = expected_log_lik - KL_div_vb_prior
#   return(ELBO)
# }
