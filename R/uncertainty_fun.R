# Internal helper --------------------------------------------------------------

#' Internal helper for uncertainty and sensitivity analysis of node risk
#'
#' Computes a vector of node-level risk scores under a Sobol' sampling design and
#' returns both the uncertainty-analysis draws and Sobol' sensitivity indices.
#'
#' The node risk score is computed as a (weighted) power mean
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
#' Internally the Sobol' design samples four independent \eqn{U(0,1)} values
#' (`a_raw`, `b_raw`, `c_raw`, `p_raw`). These are transformed before the model is
#' evaluated: the three weight draws are normalised to sum to one, yielding
#' \eqn{\alpha}, \eqn{\beta}, \eqn{\gamma}; and `p_raw` is mapped to
#' \eqn{p \in [-1, 2]}. The sensitivity indices are attributed to the
#' transformed parameters and reported under the labels `alpha`, `beta`, `gamma`,
#' and `p`.  The raw names are used only for constructing the Sobol' design matrix.
#'
#' @param params_labels Character vector of length 4 giving the labels to use in the
#'   `sensobol::sobol_indices()` output (typically `c("alpha","beta","gamma","p")`).
#'
#' @keywords internal
risk_ua_sa_fun <- function(cyclo_sc, indeg_sc, btw_sc,
                           sample_matrix, N, params, params_labels, order,
                           eps = 1e-12) {

  if (!is.numeric(eps) || length(eps) != 1L || !is.finite(eps) || eps <= 0) {
    stop("`eps` must be a single positive finite numeric value.", call. = FALSE)
  }

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

  ind <- sensobol::sobol_indices(Y = y, N = N, params = params_labels, order = order)

  list(
    ua = y[1:N],
    sa = ind,
    mu = mean(y[1:N]),
    V  = stats::var(y[1:N])
  )
}

#' Uncertainty and sensitivity analysis of node and path risk
#'
#' Runs a full variance-based uncertainty and sensitivity analysis (UA/SA) for node
#' risk scores using the results returned by [all_paths_fun()] and the functions
#' provided by the \pkg{sensobol} package (Puy et al. 2022).
#'
#' Uncertainty is induced by jointly sampling the weights
#' \eqn{(\alpha, \beta, \gamma)} (renormalized to sum to 1) and the power parameter
#' \eqn{p \in [-1, 2]} used in the node-risk definition:
#' \deqn{r =
#' \left(\alpha\,\tilde{C}^{p} + \beta\,(\tilde{d}^{\mathrm{in}})^{p} + \gamma\,\tilde{b}^{p}\right)^{1/p}\,.}
#'
#' For each node, risk scores are repeatedly recalculated using the sampled parameter
#' combinations, producing a distribution of possible outcomes. Sobol' first-order
#' and/or total-order sensitivity indices are then computed for all four parameters
#' (\eqn{\alpha}, \eqn{\beta}, \eqn{\gamma}, and \eqn{p}), quantifying how much of
#' the variance in the node risk score is attributable to each parameter.
#'
#' **Parameter labels and the Sobol' design.**
#' Internally the design samples four independent \eqn{U(0,1)} values
#' (`a_raw`, `b_raw`, `c_raw`, `p_raw`) because the Sobol' quasi-random sequence
#' requires independent uniform inputs.  Before evaluating the risk model, the raw
#' draws are transformed: the three weight draws are normalised to sum to one,
#' yielding \eqn{\alpha}, \eqn{\beta}, \eqn{\gamma}; and `p_raw` is mapped linearly
#' to \eqn{p \in [-1, 2]}.  The sensitivity indices are then attributed to the
#' *transformed* parameters and the output labels them as `alpha`, `beta`, `gamma`,
#' and `p` rather than the internal raw names, so the results are directly
#' interpretable in terms of the model parameters.
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
#' @param eps Numeric. Small positive constant \eqn{\epsilon} used for numerical stability
#'   in the \eqn{p \to 0} evaluation. Default `1e-12`.
#'
#' @details
#' For more information about the uncertainty and sensitivity analysis and the output of
#' this function, see the \pkg{sensobol} package (Puy et al. 2022).
#'
#' The returned node table includes the following columns:
#' \itemize{
#'   \item `name`: name of the node.
#'   \item `uncertainty_analysis`: numeric vector of length \eqn{N} giving the
#'     uncertainty draws in the node risk score (from Sobol matrix A).
#'   \item `sensitivity_analysis`: object returned by `sensobol::sobol_indices()`
#'     for that node, containing Sobol' indices labelled `alpha`, `beta`, `gamma`,
#'     and `p`.  These correspond to the normalised weights and the power-mean
#'     exponent respectively.  The indices are computed on the transformed
#'     parameters (see Details).
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
#'   \item{nodes}{A tibble of node results.}
#'   \item{paths}{A tibble of path results.}
#' }
#'
#' @examples
#' \donttest{
#' data(synthetic_graph)
#' out <- all_paths_fun(graph = synthetic_graph, alpha = 0.6, beta = 0.3,
#'                      gamma = 0.1, complexity_col = "cyclo")
#'
#' # Power-mean risk (increase N to at least 2^10 for a proper UA/SA)
#' results <- uncertainty_fun(all_paths_out = out, N = 2^2, order = "first")
#'
#' results$nodes
#' results$paths
#' }
#'
#' @export
#' @importFrom scales rescale
#' @importFrom sensobol sobol_matrices sobol_indices
#' @importFrom tibble tibble
uncertainty_fun <- function(all_paths_out, N, order, eps = 1e-12) {

  # ---- validate input --------------------------------------------------------
  validate_all_paths_out(all_paths_out,
                         required_cols = c("name", "cyclomatic_complexity", "indeg", "btw"))

  nodes_tbl <- all_paths_out$nodes
  paths_tbl <- all_paths_out$paths

  required_paths <- c("path_id", "path_nodes", "path_str", "hops")
  missing_paths <- setdiff(required_paths, names(paths_tbl))
  if (length(missing_paths) > 0) {
    stop("`all_paths_out$paths` is missing required columns: ",
         paste(missing_paths, collapse = ", "),
         call. = FALSE)
  }

  # ---- rescale node metrics to [0,1] -----------------------------------------
  cyclo_sc <- scales::rescale(nodes_tbl[["cyclomatic_complexity"]])
  indeg_sc <- scales::rescale(nodes_tbl[["indeg"]])
  btw_sc   <- scales::rescale(nodes_tbl[["btw"]])

  # ---- Sobol sampling --------------------------------------------------------

  # params_raw: column names used internally in the Sobol' design matrix.
  # The three weight draws (a_raw, b_raw, c_raw) are independent U(0,1) samples
  # that are normalised to sum to one, yielding alpha, beta, gamma.  p_raw is a
  # U(0,1) draw mapped to p in [-1, 2].  Sensitivity indices are reported under
  # the transformed-parameter labels (params_labels) so the output is
  # interpretable directly in terms of alpha, beta, gamma and p.
  params_raw    <- c("a_raw", "b_raw", "c_raw", "p_raw")
  params_labels <- c("alpha", "beta", "gamma", "p")

  mat <- sensobol::sobol_matrices(N = N, params = params_raw, order = order)

  s <- rowSums(mat[, c("a_raw", "b_raw", "c_raw")])
  alpha <- mat[, "a_raw"] / s
  beta  <- mat[, "b_raw"] / s
  gamma <- mat[, "c_raw"] / s

  # Map Sobol p_raw in [0,1] to p in [-1,2]
  p <- -1 + 3 * mat[, "p_raw"]

  mat <- cbind(mat, alpha = alpha, beta = beta, gamma = gamma, p = p)

  n_nodes <- nrow(nodes_tbl)
  ua_list <- vector("list", n_nodes)
  sa_list <- vector("list", n_nodes)

  for (i in seq_len(n_nodes)) {
    tmp <- risk_ua_sa_fun(
      cyclo_sc = cyclo_sc[i],
      indeg_sc = indeg_sc[i],
      btw_sc   = btw_sc[i],
      sample_matrix = mat,
      N = N,
      params = params_raw,
      params_labels = params_labels,
      order = order,
      eps = eps
    )
    ua_list[[i]] <- tmp[["ua"]]
    sa_list[[i]] <- tmp[["sa"]]
  }

  node_names <- nodes_tbl[["name"]]

  # ---- propagate UA draws to paths ------------------------------------------
  path_prob_from_nodes <- function(risks) 1 - prod(1 - risks)

  # helper: Gini per draw (column-wise)
  gini_per_draw_from_nodes <- function(ua_mat) {
    apply(ua_mat, 2, gini_index_fun)
  }

  n_paths <- nrow(paths_tbl)
  P_k        <- vector("list", n_paths)
  gini_list  <- vector("list", n_paths)
  trend_list <- vector("list", n_paths)

  for (i in seq_len(n_paths)) {

    path_nodes_i <- paths_tbl[["path_nodes"]][[i]]
    idx <- match(path_nodes_i, node_names)
    idx <- idx[!is.na(idx)]

    if (length(idx) == 0) {
      len_draws <- length(ua_list[[1]])
      P_k[[i]]        <- rep(NA_real_, len_draws)
      gini_list[[i]]  <- rep(NA_real_, len_draws)
      trend_list[[i]] <- rep(NA_real_, len_draws)
      next
    }

    # rows = nodes on path, cols = draws
    ua_mat <- do.call(rbind, ua_list[idx])
    n_draws <- ncol(ua_mat)

    # path risk per draw
    P_k[[i]] <- apply(ua_mat, 2, path_prob_from_nodes)

    # Gini and slope per draw (vectors of length n_draws)
    gini_list[[i]] <- gini_per_draw_from_nodes(ua_mat)

    # slope per draw using slope_fun
    if (nrow(ua_mat) < 2L) {
      trend_list[[i]] <- rep(0, n_draws)
    } else {
      trend_list[[i]] <- apply(ua_mat, 2, slope_fun)
    }
  }

  # ---- final outputs as tibbles ----------------------------------------------
  node_out <- tibble::tibble(
    name = node_names,
    uncertainty_analysis = ua_list,
    sensitivity_analysis = sa_list
  )

  paths_out <- tibble::tibble(
    path_id = paths_tbl[["path_id"]],
    path_str = paths_tbl[["path_str"]],
    hops = paths_tbl[["hops"]],
    uncertainty_analysis = P_k,
    gini_index = gini_list,
    risk_trend = trend_list
  )

  list(nodes = node_out, paths = paths_out)
}
