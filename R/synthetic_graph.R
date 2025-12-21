#' Synthetic citation graph for software risk examples
#'
#' A synthetic directed graph with cyclomatic complexity values to
#' illustrate the functions of the package.
#'
#' The graph is stored as a \code{tbl_graph} object with:
#' \itemize{
#'   \item Node attributes: \code{name}, \code{cyclo}
#'   \item Directed edges defined by \code{from} → \code{to}
#' }
#'
#' @format
#' A \code{tbl_graph} with:
#' \describe{
#'   \item{nodes}{55 nodes}
#'   \item{edges}{122 directed edges}
#' }
#'
#' @usage data(synthetic_graph)
#'
#' @keywords datasets
#'
"synthetic_graph"
