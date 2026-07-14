
# FUNCTION TO SELECT A BUDGETED PORTFOLIO OF NODE FIXES ########################
################################################################################
################################################################################

#' Greedy selection of the most risk-reducing node fixes
#'
#' Selects, under a budget of `budget` refactoring interventions, the set of
#' nodes whose fixing (risk set to 0) most reduces the path-level risk of the
#' system. Whereas [path_fix_heatmap()] shows the effect of fixing single
#' nodes on single paths, this function answers the budgeted question
#' directly: *given resources to refactor* \eqn{k} *functions, which*
#' \eqn{k}*?*
#'
#' At each step the algorithm evaluates, for every remaining candidate node,
#' the objective obtained by adding that node to the already-selected set,
#' fixes the best node, and repeats. Path risk is recomputed under the
#' independence assumption used throughout the package:
#' \deqn{P_k = 1 - \prod_{i=1}^{n_k} (1 - r_{k(v_i)})\,,}
#' with \eqn{r = 0} for fixed nodes. Two objectives are available:
#' \itemize{
#'   \item `"total"`: minimize \eqn{\sum_k P_k}, the total path risk. This
#'     objective is monotone submodular in the set of fixed nodes, so the
#'     greedy solution is guaranteed to achieve at least
#'     \eqn{1 - 1/e \approx 63\%} of the reduction of the optimal set of the
#'     same size (Nemhauser et al. 1978).
#'   \item `"max"`: minimize \eqn{\max_k P_k}, the risk of the riskiest path.
#' }
#'
#' The greedy search stops early if no remaining node yields an improvement
#' larger than `tol`.
#'
#' @param all_paths_out A list produced by [all_paths_fun()] with elements
#'   `nodes` and `paths`. `nodes` must contain `name` and `risk_score`;
#'   `paths` must contain `path_nodes` and `path_risk_score`.
#' @param budget Integer. Maximum number of nodes to fix. Default `5`.
#' @param objective Character scalar, `"total"` (default) or `"max"`. See
#'   Details.
#' @param tol Numeric. Minimum improvement in the objective required to
#'   continue selecting nodes. Default `1e-12`.
#'
#' @return A list with two elements:
#' \itemize{
#'   \item `portfolio`: a tibble with one row per selected node, in selection
#'     order, with columns `step`, `node`, `objective_after` (value of the
#'     objective once the node is fixed), `delta` (improvement contributed by
#'     the node) and `cum_reduction` (cumulative fraction of the initial
#'     objective removed).
#'   \item `plot`: a \pkg{ggplot2} object showing the objective as a function
#'     of the number of fixed nodes, with the selected nodes labelled.
#' }
#'
#' @references
#' Nemhauser, G. L., Wolsey, L. A., and Fisher, M. L. (1978).
#' *An analysis of approximations for maximizing submodular set functions--I*.
#' Mathematical Programming, 14(1), 265--294. doi:10.1007/BF01588971
#'
#' @examples
#' data(synthetic_graph)
#' out <- all_paths_fun(graph = synthetic_graph, alpha = 0.6, beta = 0.3,
#'                      gamma = 0.1, complexity_col = "cyclo")
#' res <- fix_portfolio_fun(out, budget = 5)
#' res$portfolio
#' res$plot
#'
#' @seealso [path_fix_heatmap()] for the node-by-path improvement matrix and
#'   [node_exposure_fun()] for path-aware node criticality.
#'
#' @export
#' @importFrom tibble tibble
#' @importFrom rlang .data
#' @importFrom ggplot2 ggplot aes geom_line geom_point geom_text labs
fix_portfolio_fun <- function(all_paths_out, budget = 5,
                              objective = c("total", "max"),
                              tol = 1e-12) {

  # ---- validate input --------------------------------------------------------

  validate_all_paths_out(all_paths_out,
                         required_cols = c("name", "risk_score"))
  objective <- match.arg(objective)

  if (!is.numeric(budget) || length(budget) != 1L || budget < 1) {
    stop("`budget` must be a positive scalar integer.", call. = FALSE)
  }
  budget <- as.integer(budget)

  nodes_tbl <- all_paths_out$nodes
  paths_tbl <- all_paths_out$paths

  if (nrow(paths_tbl) == 0) {
    stop("`all_paths_out$paths` is empty; there are no paths to improve.",
         call. = FALSE)
  }

  required_paths <- c("path_nodes", "path_risk_score")
  missing_paths <- setdiff(required_paths, names(paths_tbl))
  if (length(missing_paths) > 0) {
    stop("`all_paths_out$paths` is missing required columns: ",
         paste(missing_paths, collapse = ", "),
         call. = FALSE)
  }

  # ---- set up the greedy search ------------------------------------------------

  r_map <- stats::setNames(nodes_tbl[["risk_score"]], nodes_tbl[["name"]])

  # per-path vector of node risks (updated in place as nodes get fixed)
  path_risks <- lapply(paths_tbl[["path_nodes"]], function(v) r_map[v])

  path_prob <- function(r) 1 - prod(1 - r)
  obj_fun <- if (objective == "total") {
    function(P) sum(P)
  } else {
    function(P) max(P)
  }

  P <- vapply(path_risks, path_prob, numeric(1))
  obj0 <- obj_fun(P)
  obj_current <- obj0

  # candidates: nodes appearing on at least one path
  path_of_node <- split(
    rep(seq_along(path_risks), lengths(paths_tbl[["path_nodes"]])),
    unlist(paths_tbl[["path_nodes"]], use.names = FALSE)
  )
  candidates <- names(path_of_node)

  steps <- vector("list", budget)

  for (s in seq_len(budget)) {

    if (length(candidates) == 0) break

    best_node  <- NA_character_
    best_obj   <- obj_current
    best_P_new <- NULL

    for (v in candidates) {
      idx <- path_of_node[[v]]
      P_new_v <- vapply(idx, function(k) {
        r <- path_risks[[k]]
        r[names(r) == v] <- 0
        path_prob(r)
      }, numeric(1))

      P_try <- P
      P_try[idx] <- P_new_v
      obj_try <- obj_fun(P_try)

      if (obj_try < best_obj - tol) {
        best_node  <- v
        best_obj   <- obj_try
        best_P_new <- P_try
      }
    }

    if (is.na(best_node)) {
      message("Stopping early at step ", s,
              ": no remaining node improves the objective.")
      break
    }

    # commit the fix
    for (k in path_of_node[[best_node]]) {
      r <- path_risks[[k]]
      r[names(r) == best_node] <- 0
      path_risks[[k]] <- r
    }
    P <- best_P_new
    delta <- obj_current - best_obj
    obj_current <- best_obj

    steps[[s]] <- tibble::tibble(
      step = s,
      node = best_node,
      objective_after = obj_current,
      delta = delta,
      cum_reduction = (obj0 - obj_current) / obj0
    )

    candidates <- setdiff(candidates, best_node)
  }

  steps <- steps[!vapply(steps, is.null, logical(1))]
  portfolio <- do.call(rbind, steps)

  if (is.null(portfolio)) {
    stop("No node improves the objective; nothing to select.", call. = FALSE)
  }

  # ---- plot: objective vs number of fixed nodes ---------------------------------

  curve_tbl <- tibble::tibble(
    step = c(0L, portfolio$step),
    objective = c(obj0, portfolio$objective_after),
    node = c("", portfolio$node)
  )

  y_lab <- if (objective == "total") {
    expression(sum(P[k], k))
  } else {
    expression(max(P[k]))
  }

  fig <- ggplot2::ggplot(curve_tbl,
                         ggplot2::aes(x = .data$step, y = .data$objective)) +
    ggplot2::geom_line(colour = "grey60") +
    ggplot2::geom_point(size = 1.5) +
    ggplot2::geom_text(ggplot2::aes(label = .data$node),
                       size = 2.3, vjust = -0.9, hjust = 0.2) +
    ggplot2::labs(x = "Number of fixed nodes", y = y_lab) +
    theme_AP()

  list(portfolio = portfolio, plot = fig)
}
