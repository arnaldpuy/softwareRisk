
# FUNCTION TO COMPUTE PATH-AWARE NODE CRITICALITY ##############################
################################################################################
################################################################################

#' Path-aware node criticality (risk load)
#'
#' Summarizes, for each node, its participation in the entry-to-sink paths
#' enumerated by [all_paths_fun()]. Whereas the node risk score captures how
#' error-prone a function is in isolation, the *risk load* captures how many
#' risky execution chains depend on it: a function of moderate complexity can
#' still be critical if most high-risk paths route through it.
#'
#' For node \eqn{v} the function computes:
#' \itemize{
#'   \item `n_paths`: the number of paths that contain \eqn{v},
#'   \item `path_share`: `n_paths` divided by the total number of paths,
#'   \item `risk_load`: \eqn{\sum_{k : v \in k} P_k}, the sum of
#'     `path_risk_score` over the paths containing \eqn{v},
#'   \item `risk_load_share`: `risk_load` divided by \eqn{\sum_k P_k},
#'   \item `mean_path_risk`: the mean `path_risk_score` of the paths
#'     containing \eqn{v}.
#' }
#'
#' The output also reports the rank of each node by `risk_load` and by its
#' standalone `risk_score`, so that structurally exposed functions can be
#' told apart from merely complex ones: a node ranking much higher on
#' `risk_load` than on `risk_score` is a chokepoint rather than a hotspot.
#'
#' @param all_paths_out A list produced by [all_paths_fun()] with elements
#'   `nodes` and `paths`. `nodes` must contain `name` and `risk_score`;
#'   `paths` must contain `path_nodes` and `path_risk_score`.
#'
#' @return A tibble with one row per node appearing on at least one path,
#'   sorted by decreasing `risk_load`, with columns `name`, `risk_score`,
#'   `n_paths`, `path_share`, `risk_load`, `risk_load_share`,
#'   `mean_path_risk`, `rank_risk_load` and `rank_risk_score`. Nodes on no
#'   path (e.g., disconnected helpers) are omitted. If `all_paths_out$paths`
#'   is empty, an empty tibble with these columns is returned.
#'
#' @examples
#' data(synthetic_graph)
#' out <- all_paths_fun(graph = synthetic_graph, alpha = 0.6, beta = 0.3,
#'                      gamma = 0.1, complexity_col = "cyclo")
#' node_exposure_fun(out)
#'
#' @seealso [path_fix_heatmap()] for the per-path effect of fixing single
#'   nodes, and [fix_portfolio_fun()] for a budgeted selection of fixes.
#'
#' @export
#' @importFrom tibble tibble
node_exposure_fun <- function(all_paths_out) {

  # ---- validate input --------------------------------------------------------

  validate_all_paths_out(all_paths_out,
                         required_cols = c("name", "risk_score"))

  nodes_tbl <- all_paths_out$nodes
  paths_tbl <- all_paths_out$paths

  empty_out <- tibble::tibble(
    name = character(0),
    risk_score = numeric(0),
    n_paths = integer(0),
    path_share = numeric(0),
    risk_load = numeric(0),
    risk_load_share = numeric(0),
    mean_path_risk = numeric(0),
    rank_risk_load = integer(0),
    rank_risk_score = integer(0)
  )

  if (nrow(paths_tbl) == 0) {
    return(empty_out)
  }

  required_paths <- c("path_nodes", "path_risk_score")
  missing_paths <- setdiff(required_paths, names(paths_tbl))
  if (length(missing_paths) > 0) {
    stop("`all_paths_out$paths` is missing required columns: ",
         paste(missing_paths, collapse = ", "),
         call. = FALSE)
  }

  # ---- accumulate per-node participation --------------------------------------

  path_nodes <- paths_tbl[["path_nodes"]]
  path_risk  <- paths_tbl[["path_risk_score"]]
  n_total    <- length(path_nodes)
  total_risk <- sum(path_risk)

  # one row per (path, node) occurrence; a node cannot repeat within a
  # simple path
  node_occ <- rep(seq_len(n_total), lengths(path_nodes))
  node_nm  <- unlist(path_nodes, use.names = FALSE)

  n_paths   <- tapply(node_occ, node_nm, length)
  risk_load <- tapply(path_risk[node_occ], node_nm, sum)

  nm <- names(n_paths)

  out <- tibble::tibble(
    name = nm,
    risk_score = nodes_tbl$risk_score[match(nm, nodes_tbl$name)],
    n_paths = as.integer(n_paths),
    path_share = as.numeric(n_paths) / n_total,
    risk_load = as.numeric(risk_load),
    risk_load_share = if (total_risk > 0) as.numeric(risk_load) / total_risk
                      else rep(0, length(nm)),
    mean_path_risk = as.numeric(risk_load) / as.numeric(n_paths)
  )

  out <- out[order(-out$risk_load, out$name), , drop = FALSE]
  out$rank_risk_load  <- seq_len(nrow(out))
  out$rank_risk_score <- as.integer(rank(-out$risk_score, ties.method = "min"))

  out
}
