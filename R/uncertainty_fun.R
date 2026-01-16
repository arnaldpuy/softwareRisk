# Internal helper --------------------------------------------------------------

#' Internal helper for uncertainty and sensitivity analysis of node risk
#'
#' Computes a vector of node-level risk scores under a Sobol' sampling design and
#' returns both the uncertainty-analysis draws and Sobol sensitivity indices.
#'
#' For `risk_form = "additive"`, the node risk score is
#' \deqn{r = \alpha\,\tilde{C} + \beta\,\tilde{d}^{\mathrm{in}} + \gamma\,\tilde{b}\,,}
#' where \eqn{\tilde{C}}, \eqn{\tilde{d}^{\mathrm{in}}}, and \eqn{\tilde{b}} are scaled
#' inputs and \eqn{\alpha + \beta + \gamma = 1}.
#'
#' For `risk_form = "power_mean"`, the node risk score is computed as a (weighted) power mean
#' with exponent \eqn{p}:
#' \deqn{r =
#' \left(\alpha\,\tilde{C}^{p} + \beta\,(\tilde{d}^{\mathrm{in}})^{p} + \gamma\,\tilde{b}^{p}\right)^{1/p}\,.}
#'
#' In the limit \eqn{p \to 0}, this reduces to a weighted geometric mean, implemented
#' with a small constant \eqn{\epsilon} to avoid \eqn{\log(0)}:
#' \deqn{r = \exp\left(\alpha\log(\max(\tilde{C},\epsilon)) +
#' \beta\log(\max(\tilde{d}^{\mathrm{in}},\epsilon)) +
#' \gamma\log(\max(\tilde{b},\epsilon))\right)\,.}
#'
#' @keywords internal
risk_ua_sa_fun <- function(cyclo_sc, indeg_sc, btw_sc,
                           sample_matrix, N, params, order,
                           risk_form = c("additive", "power_mean"),
                           eps = 1e-12) {

  risk_form <- match.arg(risk_form)

  if (!is.numeric(eps) || length(eps) != 1L || !is.finite(eps) || eps <= 0) {
    stop("`eps` must be a single positive finite numeric value.", call. = FALSE)
  }

  if (risk_form == "additive") {

    # Extract only alpha, beta, gamma columns from the sample matrix
    mat <- sample_matrix[, c("alpha", "beta", "gamma"), drop = FALSE]

    # Risk score for each row
    y <- mat[, 1] * cyclo_sc + mat[, 2] * indeg_sc + mat[, 3] * btw_sc

  } else {

    # Extract only alpha, beta, gamma, p columns from the sample matrix
    mat <- sample_matrix[, c("alpha", "beta", "gamma", "p"), drop = FALSE]

    a <- mat[, 1]
    b <- mat[, 2]
    g <- mat[, 3]
    p <- mat[, 4]

    # Handle p ~ 0 using weighted geometric mean (vectorized; p can vary by row)
    tol <- 1e-12
    y <- numeric(length(p))

    idx0 <- abs(p) < tol
    if (any(idx0)) {
      y[idx0] <- exp(
        a[idx0] * log(pmax(cyclo_sc, eps)) +
          b[idx0] * log(pmax(indeg_sc, eps)) +
          g[idx0] * log(pmax(btw_sc, eps))
      )
    }

    # Standard power mean for p != 0
    idx1 <- !idx0
    if (any(idx1)) {
      y[idx1] <- (a[idx1] * (cyclo_sc^p[idx1]) +
                    b[idx1] * (indeg_sc^p[idx1]) +
                    g[idx1] * (btw_sc^p[idx1]))^(1 / p[idx1])
    }

    # numerical safety: keep in [0,1] (inputs are scaled, so this is appropriate)
    y <- pmin(1, pmax(0, y))
  }

  ind <- sensobol::sobol_indices(Y = y, N = N, params = params, order = order)

  list(
    ua = y[1:N],   # instead of y[1:(2 * N)]
    sa = ind,
    mu = mean(y[1:N]),
    V  = stats::var(y[1:N])
  )
}

#' Uncertainty and sensitivity analysis of node and path risk
#'
#' Runs a full variance-based uncertainty and sensitivity analysis for node risk
#' scores using the results returned by [all_paths_fun()] and the functions provided
#' by the \pkg{sensobol} package (Puy et al. 2022).
#'
#' To assess how sensitive the risk scores are to the choice of weighting parameters,
#' this function explores many alternative combinations of weights and (optionally)
#' a power-mean parameter, and examines how resulting risk scores vary.
#'
#' When `risk_form = "additive"`, uncertainty is induced by sampling weight triplets
#' \eqn{(\alpha, \beta, \gamma)} under the constraint \eqn{\alpha + \beta + \gamma = 1},
#' representing different plausible balances between complexity, connectivity and centrality.
#'
#' When `risk_form = "power_mean"`, uncertainty is induced by sampling both the weights
#' \eqn{(\alpha, \beta, \gamma)} (renormalized to sum to 1) and a power parameter \eqn{p}
#' used in the node-risk definition:
#' \deqn{r =
#' \left(\alpha\,\tilde{C}^{p} + \beta\,(\tilde{d}^{\mathrm{in}})^{p} + \gamma\,\tilde{b}^{p}\right)^{1/p}\,.}
#'
#' For each node, risk scores are repeatedly recalculated using the sampled parameter
#' combinations, producing a distribution of possible outcomes. This distribution is
#' then used to quantify uncertainty in the risk scores and compute Sobol' sensitivity
#' indices for each sampled parameter.
#'
#' Path-level uncertainty is obtained by propagating node-level uncertainty draws through
#' the path aggregation function:
#' \deqn{P_k = 1 - \prod_{i=1}^{n_k} (1 - r_{k(v_i)})\,,}
#' where \eqn{r_{k(v_i)}} are node risks along path \eqn{k}.
#'
#' All uncertainty metrics are computed from the first N
#' Sobol draws (matrix A), while sensitivity indices use the full Sobol' design.
#'
#' @param all_paths_out A list produced by [all_paths_fun()] with elements `nodes` and `paths`.
#'   `nodes` must contain columns `name`, `cyclomatic_complexity`, `indeg`, `btw`;
#'   `paths` must contain `path_id`, `path_nodes`, `path_str`, and `hops`.
#' @param N Integer. Base sample size used for Sobol' matrices.
#' @param order Passed to `sensobol::sobol_matrices()` and `sensobol::sobol_indices()` to control
#'   which Sobol indices are computed (e.g., first/total/second order), depending on your implementation.
#' @param risk_form Character. Risk definition used in the uncertainty/sensitivity analysis.
#'   One of `"additive"` or `"power_mean"`. Default `"additive"`.
#'
#' @details
#' For more information about the uncertainty and sensitivity analysis and the output of
#' this function, see the \pkg{sensobol} package (Puy et al. 2022).
#'
#' The returned node table includes the following columns:
#' \itemize{
#'   \item `name`: name of the node.
#'   \item `uncertainty_analysis`: numeric vector giving the uncertainty draws in the node risk score.
#'   \item `sensitivity_analysis`: object returned by `sensobol::sobol_indices()` (per-node).
#' }
#'
#' The returned paths table includes:
#' \itemize{
#'   \item `path_id`: path identifier.
#'   \item `path_str`: sequence of function calls for each path.
#'   \item `hops`: number of edges.
#'   \item `uncertainty_analysis`: numeric vector giving the uncertainty draws in the path risk score.
#'   \item `gini_index`: numeric vector giving the uncertainty draws in the gini index.
#'   \item `risk_trend`: numeric vector giving the uncertainty draws in the risk trend.
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
#' \donttest{
#' data(synthetic_graph)
#' out <- all_paths_fun(graph = synthetic_graph, alpha = 0.6, beta = 0.3,
#'                      gamma = 0.1, complexity_col = "cyclo")
#'
#' # Additive risk (increase N to at least 2^10 for a proper UA/SA)
#' results1 <- uncertainty_fun(all_paths_out = out, N = 2^2, order = "first",
#'                             risk_form = "additive")
#'
#' # Power-mean risk (increase N to at least 2^10 for a proper UA/SA)
#' results2 <- uncertainty_fun(all_paths_out = out, N = 2^2, order = "first",
#'                             risk_form = "power_mean")
#'
#' results1$nodes
#' results1$paths
#' }
#'
#' @export
#' @importFrom scales rescale
#' @importFrom data.table as.data.table
#' @importFrom rlang ":="
#' @importFrom sensobol sobol_matrices sobol_indices
uncertainty_fun <- function(all_paths_out, N, order,
                             risk_form = c("additive", "power_mean")) {

  risk_form <- match.arg(risk_form)

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

  # ---- Sobol sampling --------------------------------------------------------

  if (risk_form == "additive") {

    params <- c("a_raw", "b_raw", "c_raw")
    mat <- sensobol::sobol_matrices(N = N, params = params, order = order)
    s <- rowSums(mat)

    alpha <- mat[, "a_raw"] / s
    beta  <- mat[, "b_raw"] / s
    gamma <- mat[, "c_raw"] / s

    mat <- cbind(mat, alpha = alpha, beta = beta, gamma = gamma)

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
        order = order,
        risk_form = "additive"
      )
      ua_list[[i]] <- tmp[["ua"]]
      sa_list[[i]] <- tmp[["sa"]]
    }

  } else {

    params <- c("a_raw", "b_raw", "c_raw", "p")
    mat <- sensobol::sobol_matrices(N = N, params = params, order = order)

    s <- rowSums(mat[, c(1, 2, 3)])
    alpha <- mat[, "a_raw"] / s
    beta  <- mat[, "b_raw"] / s
    gamma <- mat[, "c_raw"] / s

    # Map Sobol p in [0,1] to p in [-1,2]
    p <- -1 + 3 * mat[, "p"]

    mat <- cbind(mat, alpha = alpha, beta = beta, gamma = gamma, p = p)

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
        order = order,
        risk_form = "power_mean"
      )
      ua_list[[i]] <- tmp[["ua"]]
      sa_list[[i]] <- tmp[["sa"]]
    }

  }

  data.table::set(node_dt, j = "uncertainty_analysis", value = ua_list)
  data.table::set(node_dt, j = "sensitivity_analysis", value = sa_list)

  # ---- propagate UA draws to paths ------------------------------------------
  path_prob_from_nodes <- function(risks) 1 - prod(1 - risks)

  # helper: slope per draw (vectorised over draws)
  slope_per_draw_from_nodes <- function(ua_mat) {
    n_nodes <- nrow(ua_mat)
    if (n_nodes < 2L) {
      # slope is undefined or zero with < 2 points; return zeros
      return(rep(0, ncol(ua_mat)))
    }

    x <- seq_len(n_nodes)
    x_mean <- mean(x)
    x_centered <- x - x_mean
    denom <- sum(x_centered^2)

    y_means <- colMeans(ua_mat)
    y_centered <- sweep(ua_mat, 2, y_means, FUN = "-")

    numer <- as.numeric(crossprod(x_centered, y_centered))
    numer / denom
  }

  # helper: Gini per draw (column-wise)
  gini_per_draw_from_nodes <- function(ua_mat) {
    apply(ua_mat, 2, gini_index_fun)
  }

  paths_dt <- data.table::as.data.table(paths_tbl)

  P_k        <- vector("list", nrow(paths_dt))
  gini_list  <- vector("list", nrow(paths_dt))
  trend_list <- vector("list", nrow(paths_dt))

  node_names <- node_dt[["name"]]

  for (i in seq_len(nrow(paths_dt))) {

    path_nodes_i <- paths_dt[["path_nodes"]][[i]]
    idx <- match(path_nodes_i, node_names)
    idx <- idx[!is.na(idx)]

    if (length(idx) == 0) {
      # no valid nodes; fall back to NA vectors using length of first UA
      len_draws <- length(node_dt[["uncertainty_analysis"]][[1]])
      P_k[[i]]        <- rep(NA_real_, len_draws)
      gini_list[[i]]  <- rep(NA_real_, len_draws)
      trend_list[[i]] <- rep(NA_real_, len_draws)
      next
    }

    # rows = nodes on path, cols = draws
    ua_mat <- do.call(rbind, node_dt[["uncertainty_analysis"]][idx])
    n_draws <- ncol(ua_mat)

    # path risk per draw
    P_k[[i]] <- apply(ua_mat, 2, path_prob_from_nodes)

    # Gini and slope per draw (vectors of length n_draws)
    gini_list[[i]]  <- gini_per_draw_from_nodes(ua_mat)
    trend_list[[i]] <- slope_per_draw_from_nodes(ua_mat)
  }

  data.table::set(paths_dt, j = "uncertainty_analysis", value = P_k)
  data.table::set(paths_dt, j = "gini_index",           value = gini_list)
  data.table::set(paths_dt, j = "risk_trend",           value = trend_list)

  # ---- final outputs ---------------------------------------------------------
  node_dt_out <- node_dt[, c("name",
                             "uncertainty_analysis",
                             "sensitivity_analysis"),
                         with = FALSE]

  paths_out   <- paths_dt[, c("path_id",
                              "path_str",
                              "hops",
                              "uncertainty_analysis",
                              "gini_index",   # list-column: vector per draw
                              "risk_trend"),  # list-column: vector per draw
                          with = FALSE]

  list(nodes = node_dt_out, paths = paths_out)
}
