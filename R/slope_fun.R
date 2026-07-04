
# FUNCTION TO COMPUTE THE SLOPE OF A NUMERIC VECTOR #########################
################################################################################
################################################################################

#' Compute the linear slope of a numeric sequence
#'
#' Computes the slope of a simple linear regression of a numeric vector
#' against its index (`seq_along(x)`). Non-finite (`NA`, `NaN`, `Inf`) values
#' are removed prior to computation, but the remaining values keep their
#' original positions as the regressor, so removing values does not distort
#' the trend. If fewer than two finite values remain, the function returns `0`.
#'
#' @param x Numeric vector.
#'
#' @details
#' The slope is estimated from the model
#' \eqn{x_i = \beta_0 + \beta_1 i + \varepsilon_i},
#' where \eqn{i = 1, \dots, n}. The function returns the estimated slope
#' \eqn{\beta_1}.
#'
#' This summary is useful for characterizing monotonic trends in ordered
#' risk values along a path.
#'
#' @return
#' A numeric scalar giving the slope of the fitted linear trend.
#'
#' @examples
#' slope_fun(c(1, 2, 3, 4))
#' slope_fun(c(4, 3, 2, 1))
#' slope_fun(c(NA, 1, 2, Inf, 3))
#'
#' @export
slope_fun <- function(x) {

  keep <- is.finite(x)
  if (sum(keep) <= 1) return(0)

  # keep the original positions as regressor so dropped values do not
  # compress the index and bias the slope
  xi <- which(keep)
  x <- x[keep]
  n <- length(x)
  as.numeric((n * sum(xi * x) - sum(xi) * sum(x)) / (n * sum(xi^2) - sum(xi)^2))
}
