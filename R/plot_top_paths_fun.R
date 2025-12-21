
#' Plot the top risky call paths on a Sugiyama layout
#'
#' Visualizes the most risky entry-to-sink paths (by decreasing `path_risk_score`)
#' computed by [all_paths_fun()]. Edges that occur on the top paths are
#' highlighted, with edge colour mapped to the mean path risk and edge width
#' mapped to the number of top paths using that edge. Nodes on the top paths
#' are emphasized, with node size mapped to in-degree and node fill mapped to
#' binned cyclomatic complexity.
#'
#' @param graph A directed `tidygraph::tbl_graph` representing the call graph to
#'   plot (typically the same graph used as input to [all_paths_fun()]).
#' @param all_paths_out Output from [all_paths_fun()], i.e. a list
#'   with elements `nodes` and `paths`. For backward compatibility, a `paths`
#'   tibble can also be supplied directly; in that case node metrics are derived
#'   from `graph` where possible.
#' @param model.name Character scalar used in the plot title (e.g., model name).
#' @param language Character scalar used in the plot title (e.g., language name).
#' @param top_n Integer. Number of highest-risk paths to display (default 10).
#' @param alpha_non_top Numeric in \deqn{[0, 1]}. Alpha (transparency) for edges that are
#'   not on the top-risk paths. Smaller values fade background edges more.
#' @details
#' The function selects the `top_n` paths by sorting `paths_tbl` on
#' `path_risk_score` (descending). For those paths, it:
#' \itemize{
#'   \item builds an edge list from `path_nodes`,
#'   \item marks graph edges that appear on at least one top path,
#'   \item computes `path_freq` (how many top paths include each edge),
#'   \item computes `risk_mean_path` (mean of `risk_sum` across top paths that
#'   include each edge),
#'   \item highlights nodes that appear on any top path.
#' }
#'
#' Node fills are based on `cyclomatic_complexity` using breaks
#' `(-Inf, 10]`, `(10, 20]`, `(20, 50]`, `(50, Inf]` as per Watson & McCabe (1996).
#'
#' This function relies on external theming/label objects `theme_AP()` and
#' `lab_expr` being available in the calling environment or package namespace.
#'
#' @references
#' Watson, A. H. and McCabe, T. J. (1996).
#' *Structured Testing: A Testing Methodology Using the Cyclomatic Complexity
#' Metric*.
#' NIST Special Publication 500-235, National Institute of Standards and
#' Technology, Gaithersburg, MD.
#' doi:10.6028/NIST.SP.500-235
#'
#' @return A `ggplot` object (invisibly). The plot is also printed as a side effect.
#'
#' @examples
#' data(synthetic_graph)
#' out <- all_paths_fun(graph = synthetic_graph, alpha = 0.6, beta = 0.3,
#' gamma = 0.1, complexity_col = "cyclo")
#' p <- plot_top_paths_fun(synthetic_graph, out, model.name = "MyModel", language = "R", top_n = 10)
#' p
#'
#' @export
#' @importFrom tidygraph activate as.igraph
#' @importFrom igraph V degree vertex_attr as_data_frame
#' @importFrom dplyr arrange desc slice_head mutate group_by summarise left_join if_else n
#' @importFrom purrr imap_dfr
#' @importFrom tibble tibble as_tibble
#' @importFrom ggraph ggraph geom_edge_link0 geom_node_point scale_edge_colour_gradient scale_edge_width
#' @importFrom ggplot2 aes scale_size_continuous scale_fill_manual guides labs theme element_blank ggtitle margin
#' @importFrom grid arrow unit
#' @importFrom rlang .data
#
plot_top_paths_fun <- function(graph,
                               all_paths_out,
                               model.name = "",
                               language = "",
                               top_n = 10,
                               alpha_non_top = 0.05) {

  # ---- Accept outputs from all_paths_fun() -----------------------------------
  # all_paths_out can be:
  #   (1) list(nodes = <tibble>, paths = <tibble>)  [recommended]
  #   (2) a paths tibble (back-compat)
  #   (3) explicit list(paths_tbl=..., nodes_tbl=...) (tolerated)

  if (is.list(all_paths_out) && !is.data.frame(all_paths_out)) {

    nodes_tbl <- all_paths_out$nodes %||% all_paths_out$nodes_tbl
    paths_tbl <- all_paths_out$paths %||% all_paths_out$paths_tbl

  } else {

    nodes_tbl <- NULL
    paths_tbl <- all_paths_out
  }

  if (is.null(paths_tbl) || !is.data.frame(paths_tbl) || nrow(paths_tbl) == 0) {
    message("No paths in `paths_tbl`; skipping plot for: ", model.name)
    return(invisible(NULL))
  }

  # If nodes_tbl is missing, derive a minimal one from graph (fallback) -------

  if (is.null(nodes_tbl) || !is.data.frame(nodes_tbl) || nrow(nodes_tbl) == 0) {
    ig_tmp <- tidygraph::as.igraph(graph)
    nodes_tbl <- tibble::tibble(
      name = igraph::V(ig_tmp)$name,
      indeg = as.numeric(igraph::degree(ig_tmp, mode = "in")),
      cyclomatic_complexity = igraph::vertex_attr(ig_tmp, "cyclomatic_complexity")
    )
  }

  # ---- Identify top risky paths ----------------------------------------------

  k <- min(top_n, nrow(paths_tbl))

  top_k_paths <- paths_tbl |>
    dplyr::arrange(dplyr::desc(.data$path_risk_score)) |>
    dplyr::slice_head(n = k)

  # ---- Build edge list from top-k paths --------------------------------------

  path_edges_all <- purrr::imap_dfr(top_k_paths$path_nodes, function(nodes_vec, pid) {
    tibble::tibble(
      from = utils::head(nodes_vec, -1),
      to   = utils::tail(nodes_vec, -1),
      path_id = pid,
      risk_sum = top_k_paths$risk_sum[pid]
    )
  })

  # ---- Edge marks: which edges appear on top paths ---------------------------

  ig2 <- tidygraph::as.igraph(graph)

  edge_df_names <- igraph::as_data_frame(ig2, what = "edges") |>
    dplyr::mutate(.edge_idx = dplyr::row_number())

  path_edges_collapsed <- path_edges_all |>
    dplyr::group_by(.data$from, .data$to) |>
    dplyr::summarise(
      path_freq = dplyr::n(),
      risk_mean_path = mean(.data$risk_sum, na.rm = TRUE),
      .groups = "drop"
    )

  edge_marks <- edge_df_names |>
    dplyr::left_join(path_edges_collapsed, by = c("from", "to")) |>
    dplyr::mutate(
      on_top_path = !is.na(.data$path_freq),
      path_freq = dplyr::if_else(is.na(.data$path_freq), 0L, as.integer(.data$path_freq)),
      risk_mean_path = dplyr::if_else(is.na(.data$risk_mean_path), 0, .data$risk_mean_path)
    )

  # ---- Annotate graph with edge + node attributes ----------------------------

  graph_sugi <- graph |>
    tidygraph::activate("edges") |>
    dplyr::mutate(
      on_top_path = edge_marks$on_top_path,
      path_freq = edge_marks$path_freq,
      risk_mean_path = edge_marks$risk_mean_path
    ) |>
    tidygraph::activate("nodes") |>
    dplyr::mutate(
      indeg = nodes_tbl$indeg[match(.data$name, nodes_tbl$name)],
      cyclomatic_complexity = nodes_tbl$cyclomatic_complexity[match(.data$name, nodes_tbl$name)],
      risk_score = nodes_tbl$risk_score[match(.data$name, nodes_tbl$name)],
      btw = nodes_tbl$btw[match(.data$name, nodes_tbl$name)],
      cyclo_class = dplyr::case_when(
        .data$cyclomatic_complexity <= 10 ~ "green",
        .data$cyclomatic_complexity <= 20 ~ "orange",
        .data$cyclomatic_complexity <= 50 ~ "red",
        .data$cyclomatic_complexity >  50 ~ "purple",
        TRUE ~ "grey"
      ),
      complexity_category = cut(
        .data$cyclomatic_complexity,
        breaks = c(-Inf, 10, 20, 50, Inf),
        labels = c("b1", "b2", "b3", "b4")
      )
    )

  risky_nodes <- unique(unlist(top_k_paths$path_nodes))

  graph_sugi <- graph_sugi |>
    tidygraph::activate("nodes") |>
    dplyr::mutate(on_top_node = .data$name %in% risky_nodes)

  # ---- Legend sizes: min/mid/max indegree among risky nodes ------------------
  node_df <- graph_sugi |>
    tidygraph::activate("nodes") |>
    tibble::as_tibble()

  risky_indeg <- node_df$indeg[node_df$on_top_node & is.finite(node_df$indeg)]
  risky_indeg <- sort(unique(risky_indeg))

  if (length(risky_indeg) == 0) {
    risky_indeg <- sort(unique(node_df$indeg[is.finite(node_df$indeg)]))
  }

  if (length(risky_indeg) >= 3) {
    min_indeg <- min(risky_indeg)
    max_indeg <- max(risky_indeg)
    mid_indeg <- risky_indeg[ceiling(length(risky_indeg) / 2)]
    legend_breaks <- c(min_indeg, mid_indeg, max_indeg)
  } else {
    legend_breaks <- risky_indeg
  }

  legend_labels <- round(legend_breaks, 1)

  # ---- Plot ------------------------------------------------------------------
  p_sugi <- ggraph::ggraph(graph_sugi, layout = "sugiyama") +
    ggraph::geom_edge_link0(
      ggplot2::aes(filter = !.data$on_top_path),
      colour = "grey80",
      alpha = alpha_non_top,
      linewidth = 0.3
    ) +
    ggraph::geom_edge_link0(
      ggplot2::aes(
        filter = .data$on_top_path,
        colour = .data$risk_mean_path,
        linewidth = pmin(pmax(.data$path_freq, 0.5), 3)
      ),
      alpha = 0.9,
      arrow = grid::arrow(length = grid::unit(1, "mm"))
    ) +
    ggraph::scale_edge_colour_gradient(
      low = "orange",
      high = "red",
      guide = "none"
    ) +
    ggraph::scale_edge_width(range = c(0.3, 2.2), guide = "none") +
    ggraph::geom_node_point(size = 1.2, colour = "#BDBDBD", alpha = 0.35, show.legend = FALSE) +
    ggraph::geom_node_point(
      ggplot2::aes(filter = .data$on_top_node, size = .data$indeg, fill = .data$complexity_category),
      shape = 21, alpha = 0.95, show.legend = TRUE
    ) +
    ggplot2::scale_size_continuous(breaks = legend_breaks, labels = legend_labels, name = "indegree") +
    ggplot2::scale_fill_manual(
      values = c("yellowgreen", "orange", "red", "purple"),
      labels =  c(b1 = expression(C %in% "(" * 0 * ", 10" * "]"),
                  b2 = expression(C %in% "(" * 10 * ", 20" * "]"),
                  b3 = expression(C %in% "(" * 20 * ", 50" * "]"),
                  b4 = expression(C %in% "(" * 50 * ", " * infinity * ")")),
      name = ""
    ) +
    theme_AP() +
    ggplot2::labs(x = "", y = "") +
    ggplot2::theme(
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      legend.position = "right",
      plot.margin = ggplot2::margin(1, 1, 1, 1),
      panel.spacing = grid::unit(0.1, "lines")
    ) +
    ggplot2::ggtitle(paste(model.name, ": ", language, sep = ""))

  print(p_sugi)
  invisible(p_sugi)
}
