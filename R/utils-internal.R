#' Validate the output of all_paths_fun
#'
#' @param x Object to validate.
#' @param required_cols Optional character vector of columns that must be present
#'   in \code{x$nodes}.
#'
#' @return Invisible \code{TRUE} on success; throws an error otherwise.
#' @keywords internal
validate_all_paths_out <- function(x, required_cols = NULL) {
  if (!is.list(x) || !all(c("nodes", "paths") %in% names(x))) {
    stop("`all_paths_out` must be the output of all_paths_fun() (a list with $nodes and $paths).",
         call. = FALSE)
  }

  if (!is.data.frame(x$nodes) || !is.data.frame(x$paths)) {
    stop("`all_paths_out$nodes` and `all_paths_out$paths` must be data.frames/tibbles.",
         call. = FALSE)
  }

  if (!is.null(required_cols)) {
    missing <- setdiff(required_cols, names(x$nodes))
    if (length(missing) > 0) {
      stop("`all_paths_out$nodes` is missing required columns: ",
           paste(missing, collapse = ", "),
           call. = FALSE)
    }
  }

  invisible(TRUE)
}
