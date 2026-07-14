
# FUNCTIONS TO ASSESS THE ROBUSTNESS OF RISK RANKINGS ##########################
################################################################################
################################################################################

#' Robustness of the risk ranking under uncertainty
#'
#' Quantifies how stable the identification of the top-\eqn{k} riskiest paths
#' (or nodes) is across the uncertainty draws produced by [uncertainty_fun()].
#' Every draw corresponds to one plausible definition of risk (one sampled
#' combination of the weights \eqn{\alpha, \beta, \gamma} and the exponent
#' \eqn{p}); ranking the paths within each draw therefore shows which paths
#' are flagged as risky *regardless* of how risk is weighted and which ones
#' enter the top-\eqn{k} only under specific assumptions.
#'
#' For each item (path or node) the function reports the probability of
#' belonging to the top-\eqn{k} across draws and the median and the 5th and
#' 95th percentiles of its rank (rank 1 = riskiest; ties share the minimum
#' rank). As a global summary it also reports the mean Spearman correlation
#' between the risk values of each draw and the consensus (mean across
#' draws) risk values: values close to 1 indicate that the ranking barely
#' responds to the risk definition, values well below 1 that conclusions
#' depend on it.
#'
#' Draws that are `NA` for an item are treated as never placing that item in
#' the top-\eqn{k}.
#'
#' @param ua_sa_out A list returned by [uncertainty_fun()] with elements
#'   `nodes` and `paths`, each containing an `uncertainty_analysis`
#'   list-column of numeric draws.
#' @param top_k Integer. Size of the top set whose membership is tracked.
#'   Clamped to the number of items with a message if larger. Default `10`.
#' @param what Character scalar, `"paths"` (default) or `"nodes"`.
#'
#' @return A list with:
#' \itemize{
#'   \item `summary`: a tibble with one row per item, sorted by decreasing
#'     `top_k_prob`, with columns `id` (path `path_id` or node `name`),
#'     `path_str` (paths only), `mean_risk`, `top_k_prob`, `rank_median`,
#'     `rank_q05` and `rank_q95`.
#'   \item `consensus_correlation`: mean Spearman correlation between each
#'     draw and the consensus risk values.
#'   \item `top_k`, `what`: the settings used.
#' }
#'
#' @examples
#' \donttest{
#' data(synthetic_graph)
#' out <- all_paths_fun(graph = synthetic_graph, alpha = 0.6, beta = 0.3,
#'                      gamma = 0.1, complexity_col = "cyclo")
#' results <- uncertainty_fun(all_paths_out = out, N = 2^7, order = "first")
#'
#' rr <- rank_robustness_fun(results, top_k = 10, what = "paths")
#' rr$summary
#' rr$consensus_correlation
#' }
#'
#' @seealso [rank_robustness_plot()] to visualize the result and
#'   [uncertainty_fun()] to generate the draws.
#'
#' @export
#' @importFrom tibble tibble
#' @importFrom stats quantile cor median sd
rank_robustness_fun <- function(ua_sa_out, top_k = 10,
                                what = c("paths", "nodes")) {

  # ---- validate input --------------------------------------------------------

  what <- match.arg(what)

  if (!is.list(ua_sa_out) || !all(c("nodes", "paths") %in% names(ua_sa_out))) {
    stop("`ua_sa_out` must be the output of uncertainty_fun() ",
         "(a list with $nodes and $paths).", call. = FALSE)
  }

  tbl <- ua_sa_out[[what]]

  if (!is.data.frame(tbl) || !"uncertainty_analysis" %in% names(tbl)) {
    stop("`ua_sa_out$", what, "` must be a data.frame with an ",
         "`uncertainty_analysis` list-column.", call. = FALSE)
  }

  if (nrow(tbl) == 0) {
    stop("`ua_sa_out$", what, "` has no rows; nothing to rank.",
         call. = FALSE)
  }

  if (!is.numeric(top_k) || length(top_k) != 1L || top_k < 1) {
    stop("`top_k` must be a positive scalar integer.", call. = FALSE)
  }
  top_k <- as.integer(top_k)

  draws <- tbl[["uncertainty_analysis"]]
  lens <- lengths(draws)
  if (length(unique(lens)) != 1L) {
    stop("All `uncertainty_analysis` vectors must have the same length.",
         call. = FALSE)
  }

  n_items <- nrow(tbl)
  if (top_k > n_items) {
    message("`top_k` (", top_k, ") exceeds the number of ", what, " (",
            n_items, "); using top_k = ", n_items, ".")
    top_k <- n_items
  }

  # rows = items, cols = draws
  M <- do.call(rbind, draws)
  n_draws <- ncol(M)

  # ---- per-draw rankings -------------------------------------------------------

  # rank 1 = riskiest; NA draws are ranked last and never make the top-k
  rank_mat <- apply(M, 2, function(x) rank(-x, ties.method = "min",
                                           na.last = "keep"))
  in_top <- rank_mat <= top_k
  in_top[is.na(in_top)] <- FALSE

  top_k_prob  <- rowMeans(in_top)
  rank_median <- apply(rank_mat, 1, stats::median, na.rm = TRUE)
  rank_q05    <- apply(rank_mat, 1, stats::quantile, probs = 0.05,
                       na.rm = TRUE, names = FALSE)
  rank_q95    <- apply(rank_mat, 1, stats::quantile, probs = 0.95,
                       na.rm = TRUE, names = FALSE)

  # ---- consensus correlation ----------------------------------------------------

  consensus <- rowMeans(M, na.rm = TRUE)
  draw_cor <- apply(M, 2, function(x) {
    ok <- is.finite(x) & is.finite(consensus)
    if (sum(ok) < 3 || stats::sd(x[ok]) == 0 || stats::sd(consensus[ok]) == 0) {
      return(NA_real_)
    }
    stats::cor(x[ok], consensus[ok], method = "spearman")
  })
  consensus_correlation <- if (all(is.na(draw_cor))) NA_real_ else
    mean(draw_cor, na.rm = TRUE)

  # ---- assemble output ------------------------------------------------------------

  id <- if (what == "paths") tbl[["path_id"]] else tbl[["name"]]

  summary_tbl <- tibble::tibble(
    id = id,
    mean_risk = consensus,
    top_k_prob = top_k_prob,
    rank_median = as.numeric(rank_median),
    rank_q05 = as.numeric(rank_q05),
    rank_q95 = as.numeric(rank_q95)
  )

  if (what == "paths" && "path_str" %in% names(tbl)) {
    summary_tbl$path_str <- tbl[["path_str"]]
    summary_tbl <- summary_tbl[, c("id", "path_str", "mean_risk",
                                   "top_k_prob", "rank_median",
                                   "rank_q05", "rank_q95")]
  }

  ord <- order(-summary_tbl$top_k_prob, -summary_tbl$mean_risk)
  summary_tbl <- summary_tbl[ord, , drop = FALSE]

  list(
    summary = summary_tbl,
    consensus_correlation = consensus_correlation,
    top_k = top_k,
    what = what,
    n_draws = n_draws
  )
}

#' Plot the robustness of the risk ranking
#'
#' Displays, for the items most often ranked in the top-\eqn{k}, the
#' probability of top-\eqn{k} membership across the uncertainty draws
#' computed by [rank_robustness_fun()]. Items with probability close to 1
#' are flagged as risky under essentially any risk definition; intermediate
#' probabilities identify items whose criticality depends on how the risk
#' score is weighted.
#'
#' @param rank_out A list returned by [rank_robustness_fun()].
#' @param top_n Integer. Number of items (by decreasing `top_k_prob`) to
#'   display. Default `20`.
#'
#' @return A \pkg{ggplot2} object.
#'
#' @examples
#' \donttest{
#' data(synthetic_graph)
#' out <- all_paths_fun(graph = synthetic_graph, alpha = 0.6, beta = 0.3,
#'                      gamma = 0.1, complexity_col = "cyclo")
#' results <- uncertainty_fun(all_paths_out = out, N = 2^7, order = "first")
#' rr <- rank_robustness_fun(results, top_k = 10, what = "paths")
#' rank_robustness_plot(rr, top_n = 20)
#' }
#'
#' @export
#' @importFrom rlang .data
#' @importFrom ggplot2 ggplot aes geom_segment geom_point labs scale_x_continuous
rank_robustness_plot <- function(rank_out, top_n = 20) {

  # ---- validate input --------------------------------------------------------

  required <- c("summary", "top_k", "what")
  if (!is.list(rank_out) || !all(required %in% names(rank_out))) {
    stop("`rank_out` must be the output of rank_robustness_fun().",
         call. = FALSE)
  }

  if (!is.numeric(top_n) || length(top_n) != 1L || top_n < 1) {
    stop("`top_n` must be a positive scalar integer.", call. = FALSE)
  }
  top_n <- as.integer(top_n)

  s <- rank_out$summary
  s <- utils::head(s, top_n)
  s$id <- factor(as.character(s$id), levels = rev(as.character(s$id)))

  x_lab <- paste0("Probability of top-", rank_out$top_k, " membership")
  y_lab <- if (rank_out$what == "paths") "Path ID" else "Node"

  ggplot2::ggplot(s, ggplot2::aes(x = .data$top_k_prob, y = .data$id)) +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = .data$top_k_prob,
                                       yend = .data$id),
                          colour = "grey70", linewidth = 0.4) +
    ggplot2::geom_point(size = 1.5) +
    ggplot2::scale_x_continuous(limits = c(0, 1)) +
    ggplot2::labs(x = x_lab, y = y_lab) +
    theme_AP()
}
