#' Estimate Individualized Empirical Null Parameters with Confounding Factors
#'
#' Estimates parameters for an individualized empirical null model that explicitly
#' accounts for confounding factors (`x_bar`). 
#'
#' @param z A numeric vector of Z-scores.
#' @param n A numeric vector of effective sample sizes corresponding to each Z-score.
#' @param x_bar A numeric matrix or vector of confounders. Each row corresponds
#'   to an observation in `z`.
#' @param p_grid A numeric vector specifying the grid of null proportions (p0)
#'   to search. Defaults to a sequence from 0.8 to 0.999.
#' @param phi_grid A numeric vector specifying the grid of the variance
#'   component (phi) to search. Defaults to a sequence from 0.01 to 5.
#' @param cutoff A numeric value specifying the quantile for the truncation
#'   interval width. Defaults to 0.95.
#' @param sigma_est Optional. A pre-estimated standard deviation. If provided,
#'   it defines a uniform width for the truncation interval.
#' @param range Optional. A parameter to define a sample-size-dependent width
#'   for the truncation interval. Takes precedence over `sigma_est`.
#'
#' @return A list containing the final estimated parameters:
#'   \item{zeta}{A numeric vector containing the estimated coefficients for the
#'   confounders in `x_bar`.}
#'   \item{phi}{The estimated variance component `phi`. Note: this is the raw
#'   `phi`, not its square root.}
#'   \item{p0}{The estimated null proportion `p0`.}
#'
#' @importFrom MASS rlm
#' @importFrom stats qnorm median IQR optim optimize
#' @export
empirical_null_ind_CF <- function(z, n, x_bar, p_grid = seq(0.8, 0.999, 0.005),
                                  phi_grid = seq(0.01, 5, 0.01),
                                  cutoff = 0.95, sigma_est, range) {
  N <- length(z)
  c_val <- qnorm(cutoff)
  rlm_est <- MASS::rlm(z ~ x_bar:sqrt(n) - 1)
  
  if (missing(sigma_est) & missing(range)) {
    range <- (rlm_est$s^2 - 1) / median(n)
    xlim <- c_val * sqrt(1 + n * range)
  } else if (!missing(sigma_est) & missing(range)) {
    xlim <- c_val * sigma_est
  } else {
    xlim <- c_val * sqrt(1 + n * range)
  }
  
  trunc_lower <- -xlim + as.matrix(sqrt(n) * x_bar) %*% as.matrix(rlm_est$coefficients)
  trunc_upper <- xlim + as.matrix(sqrt(n) * x_bar) %*% as.matrix(rlm_est$coefficients)
  
  eval_res <- sapply(phi_grid, function(g) {
    eval_p_grid <- sapply(p_grid, function(p) {
      ifelse(is.null(dim(x_bar)),
             optimize(f = negloglik_confounder_gs, interval = c(-10, 10),
                      z = z, x_bar = x_bar, n = n, trunc_lower = trunc_lower,
                      trunc_upper = trunc_upper, phi = g, p0 = p)$objective,
             optim(par = c(rlm_est$coefficients), fn = negloglik_confounder_gs,
                   z = z, x_bar = x_bar, n = n, trunc_lower = trunc_lower,
                   trunc_upper = trunc_upper, phi = g, p0 = p,
                   hessian = FALSE, method = "Nelder-Mead")$value)
    })
    return(c(p = p_grid[which.min(eval_p_grid)], neglik = min(eval_p_grid)))
  })
  
  phi_selected <- phi_grid[which.min(eval_res["neglik", ])]
  p_selected <- eval_res["p", which.min(eval_res["neglik", ])]
  zeta <- ifelse(is.null(dim(x_bar)),
                 optimize(f = negloglik_confounder_gs, interval = c(-10, 10),
                          z = z, x_bar = x_bar, n = n, trunc_lower = trunc_lower,
                          trunc_upper = trunc_upper, phi = phi_selected, p0 = p_selected)$minimum,
                 optim(par = c(rlm_est$coefficients), fn = negloglik_confounder_gs,
                       z = z, x_bar = x_bar, n = n, trunc_lower = trunc_lower,
                       trunc_upper = trunc_upper, phi = phi_selected, p0 = p_selected,
                       hessian = FALSE, method = "Nelder-Mead")$par)
  return(list(zeta = zeta,
              phi = phi_selected,
              p0 = p_selected))
}