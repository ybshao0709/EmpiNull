#' Estimate Empirical Null Parameters with provider-level Confounding Factors 
#'
#' @param Obs Numeric vector of observed counts
#' @param Exp Numeric vector of expected counts
#' @param W TBD
#' @param p.grid Numeric vector of candidate \eqn{p_0} values for grid search.
#'   Default: \code{seq(0.4, 0.999, 0.001)}.
#' @param cutoff Numeric scalar in \code{(0,1)}; quantile used to form initial
#'   truncation bounds. Default: \code{0.975}.
#' @param family Character scalar; model family. Supported:
#'   \code{"poisson"} (default) or \code{"approx poisson"}.
#'
#' @return A list with elements:
#' \item{nu_hat}{Estimated coefficient vector \eqn{\hat{\nu}} (length \code{P}).}
#' \item{sig2_hat}{Estimated variance component \eqn{\hat{\sigma}^2}.}
#' \item{p0}{Selected null proportion \eqn{p_0} (grid minimum).}
#' \item{nu_var}{Estimated covariance matrix of \eqn{\hat{\nu}}.}
#' \item{W}{The centered design matrix used in fitting.}
#' \item{(unnamed)}{Integer indices of observations outside the truncation window
#'   (i.e., those with \code{Z_FE < aorig} or \code{Z_FE > borig}).}
#' @export
empirical_null_CF_prov <- function(Obs, Exp, W,
                                   p.grid = seq(0.4, 0.999, 0.001),
                                   cutoff = 0.975,
                                   family = "poisson") {
  W <- as.matrix(W)
  W <- W - matrix(apply(W, 2, mean), nrow = nrow(W), ncol = ncol(W), byrow = TRUE)
  
  P <- ncol(W)
  Z_FE <- (Obs - Exp) / sqrt(Exp)
  cval <- qnorm(cutoff)
  niter <- length(p.grid)
  N <- length(Z_FE)
  eval.p.grid <- rep(0, niter)
  W_pre <- W * sqrt(Exp)
  
  get_pdr <- function(Wp, arg) {
    return((as.matrix(Wp) %*% arg[1:P])[, 1])
  }
  
  get_mean_var <- function(arg, ntilde, W, family) {
    if (family == "poisson") {
      temp <- exp(get_pdr(W, arg) + arg[P + 1] / 2)
      m <- sqrt(ntilde) * (temp - 1)
      v <- temp * (1 + temp * (exp(arg[P + 1]) - 1) * ntilde)
    } else if (family == "approx poisson") {
      pdrW <- get_pdr(W, arg)
      pdrWpre <- get_pdr(W_pre, arg)
      m <- pdrWpre
      v <- 1 + pdrW + arg[P + 1] * ntilde
    }
    return(list(m, v))
  }
  
  sigma_est <- rlm(Z_FE ~ 0 + W_pre, method = "M", scale.est = "MAD", psi = psi.huber)
  nu_init <- sigma_est$coefficients
  pdrW <- get_pdr(W, nu_init)
  varphi_init <- (sigma_est$s^2 - 1 - mean(pdrW)) / mean(Exp)
  initial <- c(nu_init, varphi_init)
  
  temp <- exp(get_pdr(W, initial) + initial[P + 1] / 2)
  m_init <- sqrt(Exp) * (temp - 1)
  v_init <- temp * (1 + temp * (exp(initial[P + 1]) - 1) * Exp)
  aorig <- m_init - cval * sqrt(v_init)
  borig <- m_init + cval * sqrt(v_init)
  
  negloglik <- function(arg) {
    idx_in <- which(Z_FE >= aorig & Z_FE <= borig)
    Z_FE0 <- Z_FE[idx_in]
    Exp0 <- Exp[idx_in]
    W0 <- W[idx_in, , drop = FALSE]
    
    null_mean_var <- get_mean_var(arg, Exp0, W0, family)
    m0 <- null_mean_var[[1]]
    v0 <- null_mean_var[[2]]
    N0 <- length(Exp0)
    
    idx_out <- setdiff(seq_along(Z_FE), idx_in)
    Exp1 <- Exp[idx_out]
    W1 <- W[idx_out, , drop = FALSE]
    aorig1 <- aorig[idx_out]
    borig1 <- borig[idx_out]
    
    out_mean_var <- get_mean_var(arg, Exp1, W1, family)
    m1 <- out_mean_var[[1]]
    v1 <- out_mean_var[[2]]
    
    v0 <- pmax(v0, 1e-4)
    v1 <- pmax(v1, 1e-4)
    
    Q <- pnorm(borig1, mean = m1, sd = sqrt(v1)) - pnorm(aorig1, mean = m1, sd = sqrt(v1))
    loglik <- N0 * log(p0) +
      sum(log(dnorm(Z_FE0, mean = m0, sd = sqrt(v0)))) +
      sum(log(1 - p0 * Q))
    
    return(-loglik)
  }
  
  for (i in seq_len(niter)) {
    p0 <- p.grid[i]
    result.optim <- optim(par = initial, fn = negloglik)
    eval.p.grid[i] <- result.optim$value
  }
  
  p0 <- p.grid[which.min(eval.p.grid)]
  result.optim <- optim(par = initial, fn = negloglik)
  
  idx_in_final <- which(Z_FE >= aorig & Z_FE <= borig)
  Exp0 <- Exp[idx_in_final]
  W0 <- W[idx_in_final, , drop = FALSE]
  W_pre0 <- W_pre[idx_in_final, , drop = FALSE]
  pdrW0 <- get_pdr(W0, result.optim$par)
  Omega <- diag(1 + pdrW0 + result.optim$par[P + 1] * Exp0)
  
  XtX_inv <- solve(t(W_pre0) %*% W_pre0)
  nu_var_est <- XtX_inv %*% t(W_pre0) %*% Omega %*% W_pre0 %*% XtX_inv
  
  return(list(
    nu_hat = result.optim$par[1:P],
    sig2_hat = result.optim$par[P + 1],
    p0 = p0,
    nu_var = nu_var_est,
    W = W,
    out_idx = setdiff(seq_along(Z_FE), idx_in_final)
  ))
}