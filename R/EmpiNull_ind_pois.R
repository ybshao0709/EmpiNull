#' Estimate Empirical Null Parameters for the Poisson Exact Models
#'
#' @param Z A numeric vector of observed Z-scores.
#' @param Xbar A numeric matrix (n_clusters × p) of cluster means.
#' @param E A numeric vector, where \eqn{E_i = \sum_j \exp(X_{ij}^T \beta)}.
#' @param E2 A numeric vector, where \eqn{E_i^{(2)} = \sum_j \exp(X_{ij}^T \beta)^2}.
#' @param p_grid A numeric vector of candidate p0 values to search over 
#'   (default: \code{seq(0.8, 0.99, 0.005)}).
#' @param init Optional initial parameter vector 
#'   \code{[ zeta (length p), sigmaAlpha2, sigmaEps2 ]}. 
#'   If \code{NULL}, the function calls \code{init_trunc_pois()} to construct one.
#' @param cutoff Quantile cutoff for initializing truncation bounds 
#'   in \code{init_trunc_pois} (default = 0.95).
#'
#' @return A named list with elements:
#'   \item{p0}{Estimated null proportion.}
#'   \item{negloglik}{Minimum negative log-likelihood.}
#'   \item{zeta01..zeta_p}{Estimated regression coefficients.}
#'   \item{sigmaAlpha2}{Estimated variance component \eqn{\sigma_\alpha^2}.}
#'   \item{sigmaEps2}{Estimated variance component \eqn{\sigma_\epsilon^2}.}
#' @export
empirical_null_ind_pois <- function(Z, Xbar, E, E2,
                                    p_grid = seq(0.8,0.99,0.005),
                                    cutoff = 0.95) {
  p <- ifelse(is.null(ncol(Xbar)), 1, ncol(Xbar))

  init_info <- init_trunc_pois(Z, Xbar, E, E2, cutoff)
  trunc_l   <- init_info$trunc_lower
  trunc_u   <- init_info$trunc_upper
  
  if (is.null(init)) {
    init <- c(
      init_info$zeta_init,
      sigmaAlpha2 = init_info$sigma2_init/2,
      sigmaEps2   = init_info$sigma2_init/2
    )
  }
  
  results <- lapply(p_grid, function(p0_val) {
    fit <- optim(
      par     = init,
      fn      = negloglik_pois2,
      Z       = Z,
      Xbar    = Xbar,
      E       = E,
      E2      = E2,
      trunc_l = trunc_l,
      trunc_u = trunc_u,
      p0      = p0_val,
      method  = "L-BFGS-B",
      lower   = c(rep(-Inf,p), 1e-6, 1e-6),
      upper   = c(rep( Inf,p), Inf,   Inf)
    )
    c(p0 = p0_val, negloglik = fit$value, fit$par)
  })
  
  # 4) assemble & pick best
  df <- do.call(rbind, results)
  colnames(df) <- c("p0", "negloglik",
                    sprintf("zeta%02d", 1:p),
                    "sigmaAlpha2", "sigmaEps2")
  best <- df[ which.min(df[,"negloglik"]), ]
  as.list(best)
}




init_trunc_pois <- function(Z, Xbar, E, E2, cutoff = 0.975) {
  n      <- length(Z)
  c_val  <- qnorm(cutoff)
  
  x_pre     <- Xbar * sqrt(E)
  sigma_est <- rlm(Z ~ 0 + x_pre, method = "M", scale.est = "MAD", psi = psi.huber)
  zeta_init <- sigma_est$coefficients
  lp        <- (as.matrix(Xbar) %*% as.matrix(zeta_init))[, 1]
  varA_init <- (sigma_est$s^2 - 1 - mean(lp)) / mean(E)
  
  m_init <- sqrt(E) * lp
  v_init <- 1 + lp * E2 / E + varA_init * E
  lower  <- m_init - c_val * sqrt(v_init)
  upper  <- m_init + c_val * sqrt(v_init)
  
  list(
    zeta_init   = zeta_init,
    sigma2_init = varA_init,
    trunc_lower = lower,
    trunc_upper = upper
  )
}