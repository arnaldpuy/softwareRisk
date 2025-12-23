# Internal helper --------------------------------------------------------------

#' @keywords internal
risk_ua_sa_fun <- function(cyclo_sc, indeg_sc, btw_sc, sample_matrix, N, params, order) {

  # Extract only alpha, beta, gamma columns from the sample matrix
  mat <- sample_matrix[, c("alpha", "beta", "gamma"), drop = FALSE]

  # Risk score for each row
  y <- mat[, 1] * cyclo_sc + mat[, 2] * indeg_sc + mat[, 3] * btw_sc

  # Sobol indices
  ind <- sensobol::sobol_indices(Y = y, N = N, params = params, order = order)

  # Wrap up and output
  out <- list(ua = y[1:(2 * N)], sa = ind)
  out
}


#' Uncertainty and sensitivity analysis of node and path risk
#'
#' Runs a full variance-based uncertainty and sensitivity analysis for
#' node risk scores using the results returned by
#' [all_paths_fun()] and the functions provided by the \pkg{sensobol} package (Puy et al. 2022).
#'
#' To assess how sensitive the risk scores are to the choice of weighting
#' parameters, this function explores many alternative combinations of
#' weights in \eqn{(\alpha, \beta, \gamma)} and examines how resulting risk
#' scores vary. The weights are randomly sampled under the constraint that they
#' always sum to one, so each represents a different plausible balance between
#' complexity (cyclomatic complexity), connectivity (in-degree) and centrality (betweenness).
#'
#' For each node, risk scores are repeatedly recalculated using these
#' alternative weight combinations, producing a distribution of possible
#' outcomes. This distribution is then used to quantify both overall
#' uncertainty in the risk scores and the relative influence of each
#' weight on the results.
#'
#' @param all_paths_out A list produced by [all_paths_fun()] with elements
#'   \code{nodes} and \code{paths}. \code{nodes} must contain columns
#'   \code{name}, \code{cyclomatic_complexity}, \code{indeg}, \code{btw};
#'   \code{paths} must contain \code{path_id} and \code{path_nodes}.
#' @param N Integer. Base sample size used for Sobol' matrices.
#' @param order Passed to [sensobol::sobol_matrices()] and [sensobol::sobol_indices()] to
#'   control which Sobol indices are computed (e.g., first/total/second order),
#'   depending on your implementation.
#'
#' @details
#'
#' For more information about the uncertainty and sensitivity analysis and the output
#' of [uncertainty_fun()], see the \pkg{sensobol} package (Puy et al. 2022).
#'
#' The returned node table includes the following columns:
#' \itemize{
#'   \item \code{name}: name of the node.
#'   \item \code{uncertainty_analysis}: numeric vector showing the uncertainty in the risk scores for each node.
#'   \item \code{sensitivity_indices}: a \code{data.table} showing the results of the sensitivity analysis.
#' }
#'
#' The returned paths table includes:
#' \itemize{
#'   \item \code{path_str}: sequence of function calls for each path.
#'   \item \code{hops}: Number of edges.
#'   \item \code{uncertainty analysis}: numeric vector showing the uncertainty in the risk score for each path.
#' }
#'
#' @references
#' Puy, A., Lo Piano, S., Saltelli, A., and Levin, S. A. (2022).
#' *sensobol: An R Package to Compute Variance-Based Sensitivity Indices*.
#' Journal of Statistical Software, 102(5), 1--37.
#' doi:10.18637/jss.v102.i05
#'
#' @return A named list with:
#' \describe{
#'   \item{nodes}{A data.table of node results.}
#'   \item{paths}{A data.table of path results.}
#' }
#'
#' @examples
# synthetic_graph is a tidygraph::tbl_graph with node attribute "cyclomatic_complexity"
#' data(synthetic_graph)
#' out <- all_paths_fun(graph = synthetic_graph, alpha = 0.6, beta = 0.3,
#' gamma = 0.1, complexity_col = "cyclo")
#' results <- uncertainty_fun(all_paths_out = out, N = 2^10, order = "first")
#' results$nodes
#' results$paths
#'
#' @export
#' @importFrom scales rescale
#' @importFrom data.table as.data.table
#' @importFrom rlang ":="
#' @importFrom sensobol sobol_matrices sobol_indices
uncertainty_fun <- function(all_paths_out, N, order) {

  # ---- validate input --------------------------------------------------------
  if (!is.list(all_paths_out) || !all(c("nodes", "paths") %in% names(all_paths_out))) {
    stop("`all_paths_out` must be the output of all_paths_fun() (a list with $nodes and $paths).",
         call. = FALSE)
  }

  nodes_tbl <- all_paths_out$nodes
  paths_tbl <- all_paths_out$paths

  if (!is.data.frame(nodes_tbl) || !is.data.frame(paths_tbl)) {
    stop("`all_paths_out$nodes` and `all_paths_out$paths` must be data.frames/tibbles.",
         call. = FALSE)
  }

  required_nodes <- c("name", "cyclomatic_complexity", "indeg", "btw")
  missing_nodes <- setdiff(required_nodes, names(nodes_tbl))
  if (length(missing_nodes) > 0) {
    stop("`all_paths_out$nodes` is missing required columns: ",
         paste(missing_nodes, collapse = ", "),
         call. = FALSE)
  }

  required_paths <- c("path_id", "path_nodes", "path_str", "hops")
  missing_paths <- setdiff(required_paths, names(paths_tbl))
  if (length(missing_paths) > 0) {
    stop("`all_paths_out$paths` is missing required columns: ",
         paste(missing_paths, collapse = ", "),
         call. = FALSE)
  }

  # ---- rescale node metrics to [0,1] (NO NSE) --------------------------------
  node_dt <- data.table::as.data.table(nodes_tbl)

  data.table::set(
    node_dt, j = "cyclo_sc",
    value = scales::rescale(node_dt[["cyclomatic_complexity"]])
  )
  data.table::set(
    node_dt, j = "indeg_sc",
    value = scales::rescale(node_dt[["indeg"]])
  )
  data.table::set(
    node_dt, j = "btw_sc",
    value = scales::rescale(node_dt[["btw"]])
  )

  # ---- Sobol sampling for alpha/beta/gamma ----------------------------------
  params <- c("a_raw", "b_raw", "c_raw")
  mat <- sensobol::sobol_matrices(N = N, params = params, order = order)
  s <- rowSums(mat)

  alpha <- mat[, "a_raw"] / s
  beta  <- mat[, "b_raw"] / s
  gamma <- mat[, "c_raw"] / s

  mat <- cbind(mat, alpha = alpha, beta = beta, gamma = gamma)

  # ---- run UA/SA per node ----------------------------------------------------
  ua_list <- vector("list", nrow(node_dt))
  sa_list <- vector("list", nrow(node_dt))

  for (i in seq_len(nrow(node_dt))) {
    tmp <- risk_ua_sa_fun(
      cyclo_sc = node_dt[["cyclo_sc"]][i],
      indeg_sc = node_dt[["indeg_sc"]][i],
      btw_sc   = node_dt[["btw_sc"]][i],
      sample_matrix = mat,
      N = N,
      params = params,
      order = order
    )
    ua_list[[i]] <- tmp[["ua"]]
    sa_list[[i]] <- tmp[["sa"]]$results
  }

  data.table::set(node_dt, j = "uncertainty_analysis", value = ua_list)
  data.table::set(node_dt, j = "sensitivity_analysis", value = sa_list)

  # ---- propagate UA draws to paths ------------------------------------------
  path_prob_from_nodes <- function(risks) 1 - prod(1 - risks)

  paths_dt <- data.table::as.data.table(paths_tbl)

  P_k <- vector("list", nrow(paths_dt))

  node_names <- node_dt[["name"]]

  for (i in seq_len(nrow(paths_dt))) {

    path_nodes_i <- paths_dt[["path_nodes"]][[i]]
    idx <- match(path_nodes_i, node_names)
    idx <- idx[!is.na(idx)]

    if (length(idx) == 0) {
      P_k[[i]] <- rep(NA_real_, 2 * N)
      next
    }

    ua_mat <- do.call(rbind, node_dt[["uncertainty_analysis"]][idx])
    P_k[[i]] <- apply(ua_mat, 2, path_prob_from_nodes)
  }

  # store path UA as list-column------------------------------------------------
  data.table::set(paths_dt, j = "uncertainty_analysis", value = P_k)

  # ---- final outputs ---------------------------------------------------------
  node_dt_out <- node_dt[, c("name", "uncertainty_analysis", "sensitivity_analysis"), with = FALSE]
  paths_out   <- paths_dt[, c("path_id", "path_str", "hops", "uncertainty_analysis"), with = FALSE]

  list(nodes = node_dt_out, paths = paths_out)
}
