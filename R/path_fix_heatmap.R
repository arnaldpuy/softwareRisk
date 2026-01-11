#' Path-level improvement from fixing high-risk nodes
#'
#' Compute how much the risk score of the riskiest paths would
#' decrease if selected high-risk nodes were made perfectly reliable
#' (risk fixed to 0), and visualise the result as a heatmap.
#'
#' For each of the top \code{n_nodes} nodes ranked by \code{risk_score} and
#' the top \code{k_paths} paths ranked by \code{path_risk_score}, the function
#' sets the risk of that node to 0 along the path (for all its occurrences)
#' and recomputes the path risk score under the independence assumption,
#' using
#'
#' \deqn{P_k = 1 - \prod_{i=1}^{n_k} (1 - r_{k(v_i)})}
#'
#' The improvement
#'
#' \deqn{\Delta P_k = R_{\mathrm{orig}} - R_{\mathrm{fix}}}
#'
#' is used as the fill value in the heatmap cells.
#'
#' Bright cells indicate nodes that act as chokepoints for a given path.
#' Rows with many bright cells correspond to nodes whose refactoring would
#' improve many risky paths (global chokepoints), while columns with a few
#' very bright cells correspond to paths dominated by a single risky node.
#'
#' @param all_paths_out A list returned by [all_paths_fun()], with elements
#'   \code{nodes} and \code{paths}.
#'   \code{nodes} must contain at least the columns \code{name} and
#'   \code{risk_score}. \code{paths} must contain at least the columns
#'   \code{path_id}, \code{path_nodes} (list-column of node names) and
#'   \code{path_risk_score}.
#' @param n_nodes Integer, number of top-risk nodes (by \code{risk_score})
#'   to include as rows in the heatmap. Defaults to 20.
#' @param k_paths Integer, number of top-risk paths (by \code{path_risk_score})
#'   to include as columns in the heatmap. Defaults to 20.
#'
#' @return A list with two elements:
#' \itemize{
#'   \item \code{delta_tbl}: a tibble with columns \code{node}, \code{path_id}
#'     and \code{deltaR}, containing the reduction in path risk score
#'     when fixing the node in that path.
#'   \item \code{plot}: a \pkg{ggplot2} object containing the heatmap.
#' }
#'
#' @examples
#' data(synthetic_graph)
#' out <- all_paths_fun(graph = synthetic_graph, alpha = 0.6, beta = 0.3,
#' gamma = 0.1, complexity_col = "cyclo")
#' res <- path_fix_heatmap(all_paths_out = out, n_nodes = 20, k_paths = 20)
#' res
#'
#' @importFrom rlang .data
#' @export
path_fix_heatmap <- function(all_paths_out, n_nodes = 20, k_paths = 20) {

  # ---- validate input --------------------------------------------------------
  if (!is.list(all_paths_out) || !all(c("nodes", "paths") %in% names(all_paths_out))) {
    stop("`all_paths_out` must be the output of all_paths_fun() (a list with $nodes and $paths).",
         call. = FALSE)
  }

  node_df   <- all_paths_out$nodes
  paths_tbl <- all_paths_out$paths

  if (!is.data.frame(node_df) || !is.data.frame(paths_tbl)) {
    stop("`all_paths_out$nodes` and `all_paths_out$paths` must be data.frames/tibbles.",
         call. = FALSE)
  }

  required_nodes <- c("name", "risk_score")
  missing_nodes  <- setdiff(required_nodes, names(node_df))
  if (length(missing_nodes) > 0) {
    stop("`all_paths_out$nodes` is missing required columns: ",
         paste(missing_nodes, collapse = ", "),
         call. = FALSE)
  }

  required_paths <- c("path_id", "path_nodes", "path_risk_score")
  missing_paths  <- setdiff(required_paths, names(paths_tbl))
  if (length(missing_paths) > 0) {
    stop("`all_paths_out$paths` is missing required columns: ",
         paste(missing_paths, collapse = ", "),
         call. = FALSE)
  }

  if (!is.numeric(n_nodes) || length(n_nodes) != 1L || n_nodes <= 0) {
    stop("`n_nodes` must be a positive scalar integer.", call. = FALSE)
  }
  if (!is.numeric(k_paths) || length(k_paths) != 1L || k_paths <= 0) {
    stop("`k_paths` must be a positive scalar integer.", call. = FALSE)
  }

  n_nodes <- as.integer(n_nodes)
  k_paths <- as.integer(k_paths)

  # ---- select top nodes and paths --------------------------------------------

  ord_nodes <- order(node_df[["risk_score"]], decreasing = TRUE)
  node_df_ordered <- node_df[ord_nodes, , drop = FALSE]
  top_nodes <- utils::head(node_df_ordered, n_nodes)

  ord_paths <- order(paths_tbl[["path_risk_score"]], decreasing = TRUE)
  paths_tbl_ordered <- paths_tbl[ord_paths, , drop = FALSE]
  top_paths <- utils::head(paths_tbl_ordered, k_paths)

  if (nrow(top_nodes) == 0L) {
    stop("No nodes available after filtering. Check `risk_score` in `all_paths_out$nodes`.",
         call. = FALSE)
  }
  if (nrow(top_paths) == 0L) {
    stop("No paths available after filtering. Check `path_risk_score` in `all_paths_out$paths`.",
         call. = FALSE)
  }

  # ---- named lookup for node risk scores -------------------------------------

  r_map <- stats::setNames(node_df[["risk_score"]], node_df[["name"]])

  # ---- compute deltaR for each (node, path) pair -----------------------------

  delta_tbl <- purrr::map_dfr(seq_len(nrow(top_nodes)), function(i) {
    node_name <- top_nodes[["name"]][i]

    purrr::map_dfr(seq_len(nrow(top_paths)), function(j) {
      pid        <- top_paths[["path_id"]][j]
      nodes_vec  <- top_paths[["path_nodes"]][[j]]
      R_orig     <- top_paths[["path_risk_score"]][j]

      # node not in this path: no improvement
      if (!node_name %in% nodes_vec) {
        tibble::tibble(
          node    = node_name,
          path_id = pid,
          deltaR  = 0
        )
      } else {
        r_vec <- r_map[nodes_vec]

        idxs <- which(nodes_vec == node_name)

        r_vec_fix <- r_vec
        r_vec_fix[idxs] <- 0

        R_fix <- if (length(r_vec_fix) == 0L) 0 else 1 - prod(1 - r_vec_fix)

        tibble::tibble(
          node    = node_name,
          path_id = pid,
          deltaR  = R_orig - R_fix
        )
      }
    })
  })

  # ---- order factors for plotting --------------------------------------------

  # keep explicit factor for nodes (so highest-risk at top)
  delta_tbl[["node"]] <- factor(delta_tbl[["node"]], levels = rev(top_nodes[["name"]]))
  # keep explicit factor for paths (left-to-right ordering)
  delta_tbl[["path_id"]] <- factor(delta_tbl[["path_id"]], levels = top_paths[["path_id"]])

  # ---- build heatmap (using .data pronoun) -----------------------------------

  fig_heat <- ggplot2::ggplot(
    delta_tbl,
    ggplot2::aes(
      x    = .data$path_id,
      y    = .data$node,
      fill = .data$deltaR
    )
  ) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_viridis_c(
      option = "C",
      name   = expression(Delta * R[k]),
      trans  = "sqrt"
    ) +
    theme_AP() +
    ggplot2::theme(
      axis.text.x  = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1, size = 5),
      axis.text.y  = ggplot2::element_text(size = 5),
      axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 5)),
      axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 5))
    ) +
    ggplot2::labs(x = "Path ID", y = "Node ID")

  list(
    delta_tbl = delta_tbl,
    plot      = fig_heat
  )
}
