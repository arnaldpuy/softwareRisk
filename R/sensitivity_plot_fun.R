
# FUNCTION TO PLOT SOBOL' SENSITIVITY INDICES ##################################
################################################################################
################################################################################

#' Plot Sobol' sensitivity indices of the node risk scores
#'
#' Visualizes the Sobol' indices stored in the `sensitivity_analysis`
#' list-column returned by [uncertainty_fun()], showing how much of the
#' variance in the node risk scores is driven by each parameter of the risk
#' definition (the weights \eqn{\alpha}, \eqn{\beta}, \eqn{\gamma} and the
#' power-mean exponent \eqn{p}).
#'
#' Two displays are available:
#' \itemize{
#'   \item With `nodes = NULL` (default), the distribution of the indices
#'     across *all* nodes is shown as one boxplot per parameter, answering
#'     the system-level question of which assumption dominates the risk
#'     ranking overall.
#'   \item With `nodes` set to a character vector of node names, the indices
#'     of those nodes are shown as dodged bars with their confidence
#'     intervals, one panel per node.
#' }
#'
#' Nodes whose risk score does not vary across draws (e.g., all metrics
#' zero) yield non-finite indices; these are dropped with a message.
#'
#' @param ua_sa_out A list returned by [uncertainty_fun()] with element
#'   `nodes` containing a `sensitivity_analysis` list-column.
#' @param nodes Optional character vector of node names to display
#'   individually. Default `NULL` (aggregate view across all nodes).
#' @param index Character scalar: `"both"` (default), `"Si"` (first-order)
#'   or `"Ti"` (total-order). Which indices to display.
#'
#' @return A \pkg{ggplot2} object.
#'
#' @examples
#' \donttest{
#' data(synthetic_graph)
#' out <- all_paths_fun(graph = synthetic_graph, alpha = 0.6, beta = 0.3,
#'                      gamma = 0.1, complexity_col = "cyclo")
#' results <- uncertainty_fun(all_paths_out = out, N = 2^7, order = "first")
#'
#' # system-level view: which assumption drives the risk scores?
#' sensitivity_plot_fun(results)
#'
#' # node-level view
#' sensitivity_plot_fun(results, nodes = results$nodes$name[1:4])
#' }
#'
#' @seealso [uncertainty_fun()] to generate the indices and the
#'   \pkg{sensobol} package (Puy et al. 2022) for their estimation.
#'
#' @references
#' Puy, A., Lo Piano, S., Saltelli, A., and Levin, S. A. (2022).
#' *sensobol: An R Package to Compute Variance-Based Sensitivity Indices*.
#' Journal of Statistical Software, 102(5), 1--37.
#' doi:10.18637/jss.v102.i05
#'
#' @export
#' @importFrom rlang .data
#' @importFrom ggplot2 ggplot aes geom_boxplot geom_col geom_errorbar
#'   facet_wrap position_dodge labs scale_x_discrete scale_fill_manual
sensitivity_plot_fun <- function(ua_sa_out, nodes = NULL,
                                 index = c("both", "Si", "Ti")) {

  # ---- validate input --------------------------------------------------------

  index <- match.arg(index)

  if (!is.list(ua_sa_out) || !"nodes" %in% names(ua_sa_out)) {
    stop("`ua_sa_out` must be the output of uncertainty_fun() and contain ",
         "a `$nodes` element.", call. = FALSE)
  }

  nodes_tbl <- ua_sa_out$nodes

  if (!is.data.frame(nodes_tbl) ||
      !all(c("name", "sensitivity_analysis") %in% names(nodes_tbl))) {
    stop("`ua_sa_out$nodes` must contain columns `name` and ",
         "`sensitivity_analysis`.", call. = FALSE)
  }

  if (!is.null(nodes)) {
    if (!is.character(nodes) || length(nodes) == 0) {
      stop("`nodes` must be a character vector of node names.", call. = FALSE)
    }
    missing_nodes <- setdiff(nodes, nodes_tbl$name)
    if (length(missing_nodes) > 0) {
      stop("Nodes not found in `ua_sa_out$nodes`: ",
           paste(missing_nodes, collapse = ", "), call. = FALSE)
    }
  }

  # ---- extract the indices into one table --------------------------------------

  keep <- if (is.null(nodes)) seq_len(nrow(nodes_tbl)) else
    match(nodes, nodes_tbl$name)

  ind_list <- lapply(keep, function(i) {
    si <- nodes_tbl$sensitivity_analysis[[i]]
    if (!is.list(si) || !"results" %in% names(si)) return(NULL)
    res <- as.data.frame(si$results)
    res$name <- nodes_tbl$name[i]
    res
  })
  ind_tbl <- do.call(rbind, ind_list)

  if (is.null(ind_tbl) || nrow(ind_tbl) == 0) {
    stop("No sensitivity indices found in `ua_sa_out$nodes`.", call. = FALSE)
  }

  if (index != "both") {
    ind_tbl <- ind_tbl[ind_tbl$sensitivity == index, , drop = FALSE]
  }

  n_bad <- sum(!is.finite(ind_tbl$original))
  if (n_bad > 0) {
    message(n_bad, " non-finite index value(s) dropped (nodes with ",
            "constant risk score).")
    ind_tbl <- ind_tbl[is.finite(ind_tbl$original), , drop = FALSE]
  }

  if (nrow(ind_tbl) == 0) {
    stop("All sensitivity indices are non-finite; nothing to plot.",
         call. = FALSE)
  }

  param_labels <- c(alpha = expression(alpha), beta = expression(beta),
                    gamma = expression(gamma), p = expression(p))
  fill_values <- c(Si = "grey35", Ti = "grey75")

  # ---- plot --------------------------------------------------------------------

  if (is.null(nodes)) {

    fig <- ggplot2::ggplot(
      ind_tbl,
      ggplot2::aes(x = .data$parameters, y = .data$original,
                   fill = .data$sensitivity)
    ) +
      ggplot2::geom_boxplot(linewidth = 0.3, outlier.size = 0.4)

  } else {

    fig <- ggplot2::ggplot(
      ind_tbl,
      ggplot2::aes(x = .data$parameters, y = .data$original,
                   fill = .data$sensitivity)
    ) +
      ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.8),
                        width = 0.7)

    if (all(c("low.ci", "high.ci") %in% names(ind_tbl))) {
      fig <- fig +
        ggplot2::geom_errorbar(
          ggplot2::aes(ymin = .data$low.ci, ymax = .data$high.ci),
          position = ggplot2::position_dodge(width = 0.8),
          width = 0.25, linewidth = 0.3
        )
    }

    fig <- fig + ggplot2::facet_wrap(~ name)
  }

  fig +
    ggplot2::scale_x_discrete(labels = param_labels) +
    ggplot2::scale_fill_manual(values = fill_values, name = "") +
    ggplot2::labs(x = "Parameter", y = "Sobol' index") +
    theme_AP()
}
