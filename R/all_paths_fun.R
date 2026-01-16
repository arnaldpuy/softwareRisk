
# FUNCTION TO COMPUTE ALL PATHS AND PATH-LEVEL METRICS #########################
################################################################################
################################################################################

#' Enumerate entry-to-sink call paths and compute risk metrics at the node and path level
#'
#' Given a directed call graph (`tidygraph::tbl_graph`) with a node attribute for
#' cyclomatic complexity, this function:
#' \itemize{
#'   \item computes node-level metrics (in-degree, out-degree, betweenness),
#'   \item calculates a node risk score as a weighted combination of rescaled metrics,
#'   \item enumerates all simple paths from entry nodes (in-degree = 0) to sink
#'   nodes (out-degree = 0),
#'   \item computes path-level summaries and a path-level risk score.
#'   \item calculates a gini index and the slope of risk at the path-level.
#' }
#'
#' The normalized node metrics are computed using `scales::rescale()` and denoted
#' by a tilde \eqn{\tilde{\cdot}}.
#'
#' If `risk_form = "additive"`, the risk score for node \eqn{v_i} is computed as
#'
#' \deqn{r_{(v_i)} = \alpha\,\tilde{C}_{(v_i)} + \beta\, \tilde{d}_{(v_i)}^{\mathrm{in}} + \gamma\,\tilde{b}_{(v_i)}\,,}
#'
#' where the weights \eqn{\alpha}, \eqn{\beta} and \eqn{\gamma} reflect the relative
#' importance of complexity, coupling and network position, with the constraint
#' \eqn{\alpha + \beta + \gamma = 1}.
#'
#' If `risk_form = "power_mean"`, the risk score for node \eqn{v_i} is computed as
#' the weighted power mean of normalized metrics:
#'
#' \deqn{r_{(v_i)} =
#' \left(\alpha\,\tilde{C}_{(v_i)}^{p} + \beta\,\tilde{d}_{(v_i)}^{\mathrm{in}\,p} + \gamma\,\tilde{b}_{(v_i)}^{p}\right)^{1/p}\,,}
#'
#' where \eqn{p} is the power-mean parameter. In the limit \eqn{p \to 0}, this reduces
#' to a weighted geometric mean, implemented with a small constant \eqn{\epsilon} to
#' ensure numerical stability:
#'
#' \deqn{r_{(v_i)} =
#' \exp\left(\alpha\log(\max(\tilde{C}_{(v_i)}, \epsilon)) +
#' \beta\log(\max(\tilde{d}_{(v_i)}^{\mathrm{in}}, \epsilon)) +
#' \gamma\log(\max(\tilde{b}_{(v_i)}, \epsilon))\right)\,.}
#'
#' The path-level risk score is calculated as
#'
#' \deqn{P_k = 1 - \prod_{i=1}^{n_k} (1 - r_{k(v_i)})\,,}
#'
#' where \eqn{r_{k(v_i)}} is the risk of the \eqn{i}-th function in path \eqn{k} and
#' \eqn{n_k} is the number of functions in that path. The equation behaves like a
#' saturating OR-operator: \eqn{P_k} is at least as large as the maximum individual
#' function risk and monotonically increases as more functions on the path become risky,
#' approaching 1 when several functions have high risk scores.
#'
#' The Gini index of path \eqn{k} is computed as
#'
#' \deqn{G_k = \frac{\sum_i \sum_j |r_{k(v_i)} - r_{k(v_j)}|}{2 n_k^2 \bar{r}_k}\,,}
#'
#' where \eqn{\bar{r}_k} is the mean node risk in path \eqn{k}.
#'
#' Finally, the trend of risk is defined by the slope of the regression
#'
#' \deqn{r_{k(v_i)} = \theta_{0k} + \theta_{1k}\, i + \epsilon_i \,,}
#'
#' where \eqn{r_{k(v_i)}} is the risk score of the function at position \eqn{i} along
#' path \eqn{k} (ordered from upstream to downstream execution) and \eqn{\epsilon_i} is a residual term.
#'
#' @param graph A directed `tidygraph::tbl_graph`. Graph nodes must have a `name`
#'   attribute (i.e., `igraph::V(as.igraph(graph))$name`) and a numeric node
#'   attribute specified by `complexity_col`.
#' @param alpha,beta,gamma Numeric non-negative weights for the risk score,
#'   constrained such that `alpha + beta + gamma == 1` (within `weight_tol`).
#' @param risk_form Character. Form used to combine normalized node metrics into a node risk score.
#'   One of `"additive"` or `"power_mean"`. Default `"additive"`.
#' @param p Numeric scalar. Power parameter used when `risk_form = "power_mean"`.
#'   Must be finite and lie in the interval \eqn{[-1, 2]}. Default `1`.
#' @param eps Numeric. Small positive constant \eqn{\epsilon} used for numerical stability in the
#'   \eqn{p \to 0} (geometric mean) case. Default `1e-12`.
#' @param complexity_col Character scalar. Name of the node attribute containing
#'   cyclomatic complexity. Default `"cyclo"`.
#' @param weight_tol Numeric tolerance for enforcing the weight-sum constraint.
#'   Default `1e-8`.
#'
#' @details
#' The returned `paths` tibble includes `path_cc`, a list-column giving the numeric
#' vector of cyclomatic complexity values in the same order as `path_nodes` /
#' `path_str`.
#'
#' @return A named list with two tibbles:
#' \describe{
#'   \item{nodes}{Node-level metrics with columns `name`, `cyclomatic_complexity`,
#'     `indeg` (in-degree), `outdeg` (out-degree), `btw` (betweenness), `risk_score`.}
#'   \item{paths}{Path-level metrics with columns `path_id`, `path_nodes`,
#'     `path_str`, `hops`, `path_risk_score`, `path_cc`, `gini_node_risk`,
#'     `risk_slope`, `risk_mean`, `risk_sum`}
#' }
#'
#' @examples
#' # synthetic_graph is a tidygraph::tbl_graph with node attribute "cyclo"
#' data(synthetic_graph)
#'
#' # additive risk (default)
#' out1 <- all_paths_fun(
#'   graph = synthetic_graph,
#'   alpha = 0.6, beta = 0.3, gamma = 0.1,
#'   risk_form = "additive",
#'   complexity_col = "cyclo"
#' )
#'
#' # power-mean risk (p = 0 ~ weighted geometric mean)
#' out2 <- all_paths_fun(
#'   graph = synthetic_graph,
#'   alpha = 0.6, beta = 0.3, gamma = 0.1,
#'   risk_form = "power_mean",
#'   p = 0,
#'   eps = 1e-12,
#'   complexity_col = "cyclo"
#' )
#'
#' out1$nodes
#' out1$paths
#'
#'
#' @export
#' @importFrom tidygraph as.igraph
#' @importFrom igraph V degree betweenness distances all_simple_paths vertex_attr vertex_attr_names
#' @importFrom tibble tibble
#' @importFrom purrr flatten map_dfr
#' @importFrom scales rescale
#'

all_paths_fun <- function(graph,
                          alpha = 0.6,
                          beta  = 0.3,
                          gamma = 0.1,
                          risk_form = c("additive", "power_mean"),
                          p = 1,
                          eps = 1e-12,
                          complexity_col = "cyclo",
                          weight_tol = 1e-8) {

  risk_form <- match.arg(risk_form)

  # ---- validate power-mean parameter ----------------------------------------

  if (risk_form == "power_mean") {
    if (!is.numeric(p) || length(p) != 1L || !is.finite(p)) {
      stop("`p` must be a single finite numeric value.", call. = FALSE)
    }

    if (p < -1 || p > 2) {
      stop("For `risk_form = 'power_mean'`, `p` must lie in the interval [-1, 2].",
           call. = FALSE)
    }
  }

  # ---- validate graph --------------------------------------------------------

  if (!inherits(graph, "tbl_graph")) {
    stop("`graph` must be a tidygraph::tbl_graph.", call. = FALSE)
  }

  ig <- tidygraph::as.igraph(graph)

  v_names_all <- igraph::V(ig)$name
  if (is.null(v_names_all)) {
    stop("Graph nodes must have a `name` attribute.", call. = FALSE)
  }

  # ---- validate weights ------------------------------------------------------

  wsum <- alpha + beta + gamma
  if (!is.finite(wsum) || abs(wsum - 1) > weight_tol) {
    stop("Weights must satisfy `alpha + beta + gamma = 1`.", call. = FALSE)
  }
  if (any(c(alpha, beta, gamma) < 0)) {
    stop("`alpha`, `beta`, and `gamma` must be non-negative.", call. = FALSE)
  }

  # ---- cyclomatic complexity (igraph-safe) -----------------------------------

  v_attr_names <- igraph::vertex_attr_names(ig)

  if (!(complexity_col %in% v_attr_names)) {
    stop(
      "Node attribute `", complexity_col, "` not found in graph. ",
      "Available attributes are: ",
      paste(v_attr_names, collapse = ", "),
      call. = FALSE
    )
  }

  cc <- igraph::vertex_attr(ig, complexity_col)

  if (!is.numeric(cc)) {
    stop("`", complexity_col, "` must be numeric.", call. = FALSE)
  }

  # ---- node-level metrics ----------------------------------------------------

  indeg  <- igraph::degree(ig, mode = "in")
  outdeg <- igraph::degree(ig, mode = "out")

  # RAW (unnormalized) betweenness
  btw <- igraph::betweenness(
    ig,
    directed = TRUE,
    normalized = FALSE
  )

  # normalize to [0,1] once
  cc_n    <- scales::rescale(cc)
  indeg_n <- scales::rescale(indeg)
  btw_n   <- scales::rescale(btw)

  # ---- node risk score: additive vs power mean -------------------------------

  if (risk_form == "additive") {

    risk_score <- alpha * cc_n + beta * indeg_n + gamma * btw_n

  } else if (risk_form == "power_mean") {

    if (!is.finite(p)) stop("`p` must be finite for `risk_form = 'power_mean'`.", call. = FALSE)

    if (abs(p) < 1e-12) {
      # weighted geometric mean (p -> 0), with eps to avoid log(0)
      risk_score <- exp(
        alpha * log(pmax(cc_n, eps)) +
          beta  * log(pmax(indeg_n, eps)) +
          gamma * log(pmax(btw_n, eps))
      )
    } else {
      risk_score <- (alpha * (cc_n^p) + beta * (indeg_n^p) + gamma * (btw_n^p))^(1 / p)
    }

    # numerical safety: keep in [0,1]
    risk_score <- pmin(1, pmax(0, risk_score))

  }

  nodes_tbl <- tibble::tibble(
    name = v_names_all,
    cyclomatic_complexity = cc,
    indeg = as.numeric(indeg),
    outdeg = as.numeric(outdeg),
    btw = as.numeric(btw),
    risk_score = as.numeric(risk_score)
  )

  # ---- identify entry & sink nodes -------------------------------------------

  entry_ids <- which(indeg == 0)
  sink_ids  <- which(outdeg == 0)

  if (length(entry_ids) == 0 || length(sink_ids) == 0) {
    return(list(nodes = nodes_tbl, paths = tibble::tibble()))
  }

  # ---- distance matrix -------------------------------------------------------

  dist_mat <- igraph::distances(ig, v = entry_ids, to = sink_ids, mode = "out")

  pairs_st <- expand.grid(
    s = seq_along(entry_ids),
    t = seq_along(sink_ids),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  pairs_st$s <- as.integer(pairs_st$s)
  pairs_st$t <- as.integer(pairs_st$t)

  keep1 <- entry_ids[pairs_st$s] != sink_ids[pairs_st$t]
  pairs_st <- pairs_st[keep1, , drop = FALSE]

  keep2 <- dist_mat[cbind(pairs_st$s, pairs_st$t)] < Inf
  pairs_st <- pairs_st[keep2, , drop = FALSE]

  if (nrow(pairs_st) == 0) {
    return(list(nodes = nodes_tbl, paths = tibble::tibble()))
  }

  # ---- enumerate all simple paths --------------------------------------------

  all_paths_nested <- lapply(seq_len(nrow(pairs_st)), function(i) {
    igraph::all_simple_paths(
      ig,
      from = entry_ids[pairs_st$s[i]],
      to   = sink_ids[pairs_st$t[i]],
      mode = "out"
    )
  })

  all_paths <- purrr::flatten(all_paths_nested)
  if (length(all_paths) == 0) {
    return(list(nodes = nodes_tbl, paths = tibble::tibble()))
  }

  # ---- path-level metrics ----------------------------------------------------

  paths_tbl <- purrr::map_dfr(seq_along(all_paths), function(k) {

    v <- all_paths[[k]]
    v_names <- names(v)

    # risk along the path
    ps <- nodes_tbl$risk_score[match(v_names, nodes_tbl$name)]

    # cyclomatic complexity along the path (numeric vector, ordered)
    cc_vec <- nodes_tbl$cyclomatic_complexity[match(v_names, nodes_tbl$name)]

    path_risk_score <- if (length(ps) == 0) 0 else
      1 - prod(1 - ps, na.rm = TRUE)

    tibble::tibble(
      path_id = k,
      path_nodes = list(v_names),
      path_str = paste(v_names, collapse = " \u2192 "),
      hops = length(v_names) - 1,
      path_risk_score = path_risk_score,
      path_cc = list(as.numeric(cc_vec)),
      gini_node_risk = gini_index_fun(ps),
      risk_slope = slope_fun(ps),
      risk_mean = mean(ps, na.rm = TRUE),
      risk_sum = sum(ps, na.rm = TRUE)
    )
  })

  list(
    nodes = nodes_tbl,
    paths = paths_tbl
  )
}









