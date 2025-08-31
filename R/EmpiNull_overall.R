#' Estimate Parameters of an Overall Empirical Null Distribution
#'
#' This function estimates the parameters of an overall empirical null distribution
#' (mean, standard deviation, and null proportion p0) from a vector of Z-scores.
#' It uses a truncated normal distribution and optimizes the negative log-likelihood
#' over a specified grid of null proportions.
#'
#' @param z A numeric vector of Z-scores.
#' @param p_grid A numeric vector specifying the grid of null proportions (p0)
#'   to search. Defaults to a sequence from 0.8 to 0.999.
#' @param xlim An optional numeric vector of length two specifying the center
#'   and half-width (`c(center, half_width)`) of the truncation interval used
#'   to define the null distribution. If `NULL` (the default), the interval is
#'   calculated automatically based on the median and IQR of `z`.
#'
#' @return A list containing one element:
#'   \item{est}{A numeric vector of length three containing the estimated mean,
#'   standard deviation, and the selected null proportion (p0).}
#'
#' @export
empirical_null_overall <- function(z, p_grid = seq(0.8, 0.999, 0.005), xlim = NULL) {
  
  N <- length(z)
  
  if (is.null(xlim)) {
    b <- ifelse(N > 500000, 1, 4.3 * exp(-0.26 * log(N, 10)))
    xlim <- c(median(z), b * IQR(z) / (2 * qnorm(0.75)))
  }
  trunc_lower <- xlim[1] - xlim[2]
  trunc_upper <- xlim[1] + xlim[2]
  
  z0 <- z[which(z >= trunc_lower & z <= trunc_upper)]
  
  eval_p_grid <- sapply(p_grid, function(p) {
    optim(par = c(mean(z0), sd(z0)), fn = negloglik_overall,
          z = z, trunc_lower = trunc_lower,
          trunc_upper = trunc_upper, p0 = p,
          hessian = FALSE, method = "Nelder-Mead")$value
  })
  
  p0_selected <- p_grid[which.min(eval_p_grid)]
  res_optim <- optim(par = c(mean(z0), sd(z0)), fn = negloglik_overall,
                     z = z, trunc_lower = trunc_lower,
                     trunc_upper = trunc_upper, p0 = p0_selected,
                     hessian = FALSE, method = "Nelder-Mead")
  
  # Corrected to return the selected p0, not the last one from the grid.
  return(list(est = c(res_optim$par, p0_selected)))
}