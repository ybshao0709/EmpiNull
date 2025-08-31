#' Estimate Individualized Empirical Null Parameters (CRE Model)
#'
#' Estimates a null proportion (p0) and a variance component (phi) for an
#' individualized empirical null model.
#'
#' The truncation interval for defining the null set is determined with the
#' following priority:
#' 1. If `range` is provided, a sample-size-dependent interval is used.
#' 2. If `range` is missing but `sigma_est` is provided, a uniform interval is used.
#' 3. If both are missing, `range` is estimated using robust regression.
#'
#' @param z A numeric vector of Z-scores.
#' @param n A numeric vector of effective sample sizes corresponding to each Z-score.
#' @param p_grid A numeric vector specifying the grid of null proportions (p0)
#'   to search. Defaults to a sequence from 0.8 to 0.999.
#' @param phi_grid A numeric vector specifying the grid of the variance
#'   component (phi) to search. Defaults to a sequence from 0.01 to 5.
#' @param cutoff A numeric value specifying the quantile for the truncation
#'   interval width. Defaults to 0.95.
#' @param sigma_est Optional. A pre-estimated standard deviation. If provided,
#'   it defines a uniform truncation interval for all Z-scores.
#' @param range Optional. A parameter to define a sample-size-dependent
#'   truncation interval. Takes precedence over `sigma_est`.
#'
#' @return A list containing the estimated parameters:
#'   \item{phi}{The square root of the estimated variance component `phi`.}
#'   \item{p}{The estimated null proportion `p0`.}
#'
#' @importFrom MASS rlm
#' @importFrom stats qnorm median IQR
#' @export
empirical_null_ind_CF_CRE <- function(z, n, p_grid = seq(0.8, 0.999, 0.005),
                                      phi_grid = seq(0.01, 5, 0.01),
                                      cutoff = 0.95, sigma_est, range) {
  N <- length(z)
  c_val <- qnorm(cutoff)
  
  if (missing(sigma_est) & missing(range)) {
    rlm_est <- MASS::rlm(z ~ 1)
    range <- (rlm_est$s^2 - 1) / median(n)
    xlim <- c_val * sqrt(n * range)
  } else if (!missing(sigma_est) & missing(range)) {
    xlim <- c_val * sigma_est
  } else {
    xlim <- c_val * sqrt(n * range)
  }
  
  trunc_lower <- -xlim
  trunc_upper <- xlim
  
  eval_res <- sapply(phi_grid, function(g) {
    eval_p_grid <- sapply(p_grid, function(p) {
      negloglik_confounder_CRE(arg = g, z, n, trunc_lower, trunc_upper, p0 = p)
    })
    return(c(p = p_grid[which.min(eval_p_grid)], neglik = min(eval_p_grid)))
  })
  
  phi_selected <- phi_grid[which.min(eval_res["neglik", ])]
  p_selected <- eval_res["p", which.min(eval_res["neglik", ])]
  
  return(list(phi = sqrt(phi_selected), p = p_selected))
}