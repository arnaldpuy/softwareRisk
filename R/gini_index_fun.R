
# FUNCTION TO COMPUTE THE GINI INDEX OF A NUMERIC VECTOR #######################
################################################################################
################################################################################

#' Compute the Gini index of a numeric vector
#'
#' Computes the Gini index (a measure of inequality) for a numeric vector.
#' Non-finite (`NA`, `NaN`, `Inf`) values are removed prior to computation.
#' If fewer than two finite values remain, the function returns `0`.
#'
#' @param x Numeric vector.
#'
#' @details
#' The Gini index ranges from 0 (perfect equality) to 1 (maximal inequality).
#'
#' @return
#' A numeric scalar giving the Gini index of `x`.
#'
#' @examples
#' gini_index_fun(c(1, 1, 1, 1))
#' gini_index_fun(c(1, 2, 3, 4))
#' gini_index_fun(c(NA, 1, 2, Inf, 3))
#'
#' @export
#' @importFrom ineq Gini
gini_index_fun <- function(x) {

  x <- x[is.finite(x) & !is.na(x)]
  if (length(x) <= 1) return(0)

  ineq::Gini(x)
}
