#' Plot path-level uncertainty for the top-risk paths
#'
#' Plot the top \code{n_paths} paths ranked by their mean risk score,
#' with horizontal error bars representing the uncertainty range
#' (minimum and maximum risk) computed from the Monte Carlo samples
#' stored in \code{uncertainty_analysis}.
#'
#' This function is designed to work with the \code{paths} component of
#' the output of [uncertainty_fun()]. For each path, it summarises the
#' vector of path risk values by computing the mean, minimum and maximum values, and then displays
#' these summaries for the \code{n_paths} most risky paths.
#'
#' @param ua_sa_out A list returned by [uncertainty_fun()] containing at
#'   least an element \code{$paths}, which must be a data frame with
#'   columns \code{path_id} and \code{uncertainty_analysis}. The column
#'   \code{uncertainty_analysis} is expected to be a list-column where
#'   each element is a numeric vector of path risk values obtained from
#'   Monte Carlo sampling.
#' @param n_paths Integer, number of top paths (by mean risk) to include
#'   in the plot. Defaults to 20.
#'
#' @return A \pkg{ggplot2} object.
#'
#' @examples
#' data(synthetic_graph)
#' out <- all_paths_fun(graph = synthetic_graph, alpha = 0.6, beta = 0.3,
#' gamma = 0.1, complexity_col = "cyclo")
#' results <- uncertainty_fun(all_paths_out = out, N = 2^10, order = "first")
#' path_uncertainty_plot(ua_sa_out = results, n_paths = 20)
#' @importFrom rlang .data
#' @export
path_uncertainty_plot <- function(ua_sa_out, n_paths = 20) {

  # ---- validate input --------------------------------------------------------
  if (!is.list(ua_sa_out) || !"paths" %in% names(ua_sa_out)) {
    stop("`ua_sa_out` must be the output of uncertainty_fun() and contain a `$paths` element.",
         call. = FALSE)
  }

  paths_tbl <- ua_sa_out$paths

  if (!is.data.frame(paths_tbl)) {
    stop("`ua_sa_out$paths` must be a data.frame / tibble.", call. = FALSE)
  }

  required_cols <- c("path_id", "uncertainty_analysis")
  missing_cols  <- setdiff(required_cols, names(paths_tbl))
  if (length(missing_cols) > 0) {
    stop("`ua_sa_out$paths` is missing required columns: ",
         paste(missing_cols, collapse = ", "),
         call. = FALSE)
  }

  if (!is.numeric(n_paths) || length(n_paths) != 1L || n_paths <= 0) {
    stop("`n_paths` must be a positive scalar integer.", call. = FALSE)
  }
  n_paths <- as.integer(n_paths)

  ua_col <- paths_tbl[["uncertainty_analysis"]]
  if (!is.list(ua_col)) {
    stop("`uncertainty_analysis` must be a list-column of numeric vectors.",
         call. = FALSE)
  }

  # ---- compute mean, min, max from uncertainty_analysis ----------------------

  P_k_mean <- vapply(
    ua_col,
    function(x) {
      x_num <- as.numeric(x)
      mean(x_num, na.rm = TRUE)
    },
    numeric(1L)
  )

  P_k_min <- vapply(
    ua_col,
    function(x) {
      x_num <- as.numeric(x)
      min(x_num, na.rm = TRUE)
    },
    numeric(1L)
  )

  P_k_max <- vapply(
    ua_col,
    function(x) {
      x_num <- as.numeric(x)
      max(x_num, na.rm = TRUE)
    },
    numeric(1L)
  )

  paths_tbl[["P_k_mean"]] <- P_k_mean
  paths_tbl[["P_k_min"]]  <- P_k_min
  paths_tbl[["P_k_max"]]  <- P_k_max

  # ---- select top paths by mean risk -----------------------------------------

  ord <- order(paths_tbl[["P_k_mean"]], decreasing = TRUE)
  paths_ordered <- paths_tbl[ord, , drop = FALSE]
  top_paths <- utils::head(paths_ordered, n_paths)

  if (nrow(top_paths) == 0L) {
    stop("No paths available after computing summaries. Check `uncertainty_analysis`.",
         call. = FALSE)
  }

  # ---- build plot: sort y inside aes() ---------------------------------------

  p <- ggplot2::ggplot(
    top_paths,
    ggplot2::aes(
      x = .data$P_k_mean,
      y = stats::reorder(.data$path_id, .data$P_k_mean)
    )
  ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        xmin = .data$P_k_min,
        xmax = .data$P_k_max
      ),
      height = 0.2
    ) +
    ggplot2::geom_point(size = 1) +
    ggplot2::labs(
      y = "Path ID",
      x = expression(P[k])
    ) +
    theme_AP() +
    ggplot2::scale_x_continuous(
      breaks = scales::breaks_pretty(n = 3)
    ) +
    ggplot2::theme(
      axis.text.y     = ggplot2::element_text(size = 4),
      legend.position = "none"
    )

  p
}
