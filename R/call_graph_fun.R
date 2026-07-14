
# FUNCTIONS TO BUILD THE CALL GRAPH ############################################
################################################################################
################################################################################

# ---- internal: cyclomatic complexity of a function ---------------------------

#' Cyclomatic complexity of an R function
#'
#' Computes McCabe's cyclomatic complexity by walking the abstract syntax tree
#' of the function body (default arguments included). The complexity is
#' \eqn{1} plus one unit for each decision point: `if`, `for`, `while`,
#' `repeat`, `&&`, `||`, and each `switch()` alternative beyond the first.
#'
#' @param f A function.
#' @return An integer scalar (>= 1).
#' @keywords internal
#' @noRd
cyclocomp_fun <- function(f) {

  count_expr <- function(e) {
    if (is.call(e)) {
      head <- e[[1]]
      inc <- 0L
      if (is.symbol(head)) {
        op <- as.character(head)
        if (op %in% c("if", "for", "while", "repeat", "&&", "||")) {
          inc <- 1L
        } else if (op == "switch") {
          # each alternative beyond the first is an extra branch
          inc <- max(0L, length(e) - 3L)
        }
      }
      inc + sum(vapply(as.list(e)[-1], count_expr, integer(1)),
                count_expr(head))
    } else if (is.pairlist(e)) {
      sum(vapply(as.list(e), count_expr, integer(1)))
    } else {
      0L
    }
  }

  b <- body(f)
  frm <- formals(f)
  n <- 1L
  if (!is.null(b)) n <- n + count_expr(b)
  if (!is.null(frm)) {
    frm_list <- as.list(frm)
    for (i in seq_along(frm_list)) {
      # skip empty defaults (the missing-argument sentinel is an empty symbol);
      # indexed access avoids binding the sentinel to a variable, whose bare
      # lookup would error
      if (is.symbol(frm_list[[i]]) && !nzchar(as.character(frm_list[[i]]))) next
      n <- n + count_expr(frm_list[[i]])
    }
  }
  n
}

# ---- internal: extract top-level function definitions from parsed code ------

#' Extract `name <- function(...)` definitions from parsed expressions
#'
#' Handles `<-`, `=` and `<<-` assignments whose right-hand side is a
#' function literal. Later definitions of the same name overwrite earlier
#' ones. The function literals are evaluated in `baseenv()`, which creates
#' the closures without executing any package code.
#'
#' @param exprs A list of parsed expressions.
#' @return A named list of functions.
#' @keywords internal
#' @noRd
extract_fun_defs <- function(exprs) {

  defs <- list()

  for (e in exprs) {
    if (is.call(e) && length(e) == 3L && is.symbol(e[[1]]) &&
        as.character(e[[1]]) %in% c("<-", "=", "<<-")) {

      lhs <- e[[2]]
      rhs <- e[[3]]

      if (is.symbol(lhs) && is.call(rhs) && is.symbol(rhs[[1]]) &&
          as.character(rhs[[1]]) == "function") {
        # evaluating a function literal only builds the closure; the body
        # is not executed
        defs[[as.character(lhs)]] <- eval(rhs, envir = baseenv())
      }
    }
  }

  defs
}

# ---- internal: edges among a named set of functions --------------------------

#' Detect calls among a named list of functions
#'
#' Uses `codetools::findGlobals()` to detect, for each function, which
#' other functions of the set it calls. Symbols used as values (e.g., a
#' function passed to `lapply()`) are treated as calls when they match a
#' function in the set.
#'
#' @param funs Named list of functions.
#' @return A data.frame with columns `from` (caller) and `to` (callee).
#' @keywords internal
#' @noRd
detect_call_edges <- function(funs) {

  def_names <- names(funs)

  edge_list <- lapply(def_names, function(nm) {
    g <- codetools::findGlobals(funs[[nm]], merge = FALSE)
    callees <- union(
      intersect(g$functions, def_names),
      intersect(g$variables, def_names)
    )
    if (length(callees) == 0) return(NULL)
    data.frame(from = nm, to = callees, stringsAsFactors = FALSE)
  })

  edges <- do.call(rbind, edge_list)
  if (is.null(edges)) {
    edges <- data.frame(from = character(0), to = character(0),
                        stringsAsFactors = FALSE)
  }
  edges
}

#' Build a call graph from an R package or a directory of R scripts
#'
#' Constructs the directed call graph required by [all_paths_fun()]
#' automatically, so that no manual preparation of edge lists or complexity
#' spreadsheets is needed for `R` code. Each node is a function, each edge
#' `from -> to` is a call from function `from` (caller) to function `to`
#' (callee), and each node carries a `cyclo` attribute with its cyclomatic
#' complexity.
#'
#' @param pkg Character scalar. Name of an installed package whose namespace
#'   functions (exported and internal) are analyzed. Exactly one of `pkg` or
#'   `dir` must be supplied.
#' @param dir Character scalar. Path to a directory containing `.R` files.
#'   Files are parsed (searched recursively) and all top-level
#'   `name <- function(...)` definitions (also with `=` or `<<-`) are
#'   collected. The code is never executed: function literals are only
#'   parsed and turned into closures.
#' @param exclude Optional character vector of function names to drop from
#'   the graph (e.g., generated or vendored code).
#' @param keep_self_loops Logical. Keep self-recursion edges (a function
#'   calling itself)? These edges do not affect the simple-path enumeration
#'   of [all_paths_fun()] but change the in-degree. Default `FALSE`.
#'
#' @details
#' Call detection relies on `codetools::findGlobals()`, i.e., on static
#' analysis of each function body. Two limitations follow: (i) calls built
#' at run time (e.g., `do.call(paste0("f", i), ...)`, `get()`, method
#' dispatch) cannot be detected; (ii) a function of the set that is merely
#' *referenced* (e.g., passed to `lapply()`) is treated as called, since
#' in a dependency sense the caller relies on it.
#'
#' Cyclomatic complexity is computed internally as \eqn{1} plus one unit per
#' decision point (`if`, `for`, `while`, `repeat`, `&&`, `||`, and each
#' `switch()` alternative beyond the first), following McCabe (1976).
#'
#' @references
#' McCabe, T. J. (1976). *A Complexity Measure*.
#' IEEE Transactions on Software Engineering, SE-2(4), 308--320.
#' doi:10.1109/TSE.1976.233837
#'
#' @return A directed `tidygraph::tbl_graph` with node attributes `name` and
#'   `cyclo`, ready to be passed to [all_paths_fun()].
#'
#' @seealso [read_call_graph()] to import a call graph prepared outside `R`
#'   (e.g., for Fortran, C or Python models).
#'
#' @examples
#' # build the call graph of a directory of R scripts
#' td <- file.path(tempdir(), "cg_example")
#' dir.create(td, showWarnings = FALSE)
#' writeLines(c(
#'   "load_data <- function(x) x",
#'   "clean_data <- function(x) load_data(x)",
#'   "calc_scores <- function(x) if (length(x) > 0) mean(x) else 0",
#'   "compute_risk <- function(x) calc_scores(clean_data(x))"
#' ), file.path(td, "model.R"))
#'
#' g <- call_graph_fun(dir = td)
#' g
#'
#' # the result feeds directly into the pipeline
#' out <- all_paths_fun(g, alpha = 0.6, beta = 0.3, gamma = 0.1)
#' out$paths
#'
#' # build the call graph of an installed package
#' \donttest{
#' g_pkg <- call_graph_fun(pkg = "softwareRisk")
#' g_pkg
#' }
#'
#' @export
#' @importFrom tibble tibble
#' @importFrom tidygraph tbl_graph
call_graph_fun <- function(pkg = NULL, dir = NULL, exclude = NULL,
                           keep_self_loops = FALSE) {

  # ---- validate input ---------------------------------------------------------

  if (is.null(pkg) == is.null(dir)) {
    stop("Supply exactly one of `pkg` or `dir`.", call. = FALSE)
  }

  if (!is.null(exclude) && !is.character(exclude)) {
    stop("`exclude` must be a character vector of function names.",
         call. = FALSE)
  }

  # ---- collect the functions --------------------------------------------------

  if (!is.null(pkg)) {

    if (!is.character(pkg) || length(pkg) != 1L) {
      stop("`pkg` must be a single package name.", call. = FALSE)
    }
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Package `", pkg, "` is not installed.", call. = FALSE)
    }

    ns <- asNamespace(pkg)
    obj_names <- ls(ns, all.names = TRUE)
    funs <- mget(obj_names, envir = ns, mode = "function",
                 ifnotfound = list(NULL))
    funs <- Filter(function(f) is.function(f) && !is.primitive(f), funs)

  } else {

    if (!is.character(dir) || length(dir) != 1L || !dir.exists(dir)) {
      stop("`dir` must be the path to an existing directory.", call. = FALSE)
    }

    files <- list.files(dir, pattern = "\\.[Rr]$", full.names = TRUE,
                        recursive = TRUE)
    if (length(files) == 0) {
      stop("No `.R` files found in `", dir, "`.", call. = FALSE)
    }

    exprs <- unlist(lapply(files, function(f) as.list(parse(f))),
                    recursive = FALSE)
    funs <- extract_fun_defs(exprs)
  }

  if (!is.null(exclude)) {
    funs <- funs[setdiff(names(funs), exclude)]
  }

  if (length(funs) == 0) {
    stop("No function definitions found.", call. = FALSE)
  }

  # ---- edges and complexity ---------------------------------------------------

  edges <- detect_call_edges(funs)

  if (!keep_self_loops) {
    edges <- edges[edges$from != edges$to, , drop = FALSE]
  }

  cyclo <- vapply(funs, cyclocomp_fun, integer(1))

  nodes_tbl <- tibble::tibble(
    name = names(funs),
    cyclo = as.numeric(cyclo)
  )

  message("Call graph: ", nrow(nodes_tbl), " functions, ",
          nrow(edges), " calls.")

  tidygraph::tbl_graph(nodes = nodes_tbl, edges = edges, directed = TRUE)
}

#' Import a call graph from edge-list and complexity tables
#'
#' Builds the directed call graph required by [all_paths_fun()] from the two
#' datasets described in the package vignette: an edge list of function calls
#' and a table of per-function cyclomatic complexity. The input is validated
#' so that common preparation errors (misspelled columns, functions missing
#' from the complexity table, duplicated or non-numeric entries) fail early
#' with an informative message. This is the recommended entry point for
#' models written in languages other than `R` (e.g., Fortran, C, Python),
#' whose edge lists and complexity values are extracted with external tools.
#'
#' @param edges A data.frame (or path to a `.csv` file) with one row per
#'   function call. Must contain the columns given by `from_col` (caller)
#'   and `to_col` (callee).
#' @param metrics A data.frame (or path to a `.csv` file) with one row per
#'   function. Must contain the columns given by `name_col` (function name)
#'   and `complexity_col` (cyclomatic complexity).
#' @param from_col,to_col Character scalars. Names of the caller and callee
#'   columns in `edges`. Defaults `"from"` and `"to"`.
#' @param name_col Character scalar. Name of the function-name column in
#'   `metrics`. Default `"name"`.
#' @param complexity_col Character scalar. Name of the complexity column in
#'   `metrics`. Default `"cyclo"`. The column keeps this name in the
#'   returned graph, matching the default of [all_paths_fun()].
#'
#' @return A directed `tidygraph::tbl_graph` with node attributes `name` and
#'   `complexity_col`, ready to be passed to [all_paths_fun()].
#'
#' @seealso [call_graph_fun()] to build the graph automatically from `R`
#'   source code.
#'
#' @examples
#' calls_df <- data.frame(
#'   from = c("clean_data", "compute_risk", "compute_risk", "calc_scores"),
#'   to   = c("load_data", "clean_data", "calc_scores", "trim_mean")
#' )
#' cyclo_df <- data.frame(
#'   name  = c("clean_data", "load_data", "compute_risk", "calc_scores",
#'             "trim_mean"),
#'   cyclo = c(6, 3, 12, 5, 2)
#' )
#'
#' g <- read_call_graph(edges = calls_df, metrics = cyclo_df)
#' g
#'
#' @export
#' @importFrom tibble tibble
#' @importFrom tidygraph tbl_graph
#' @importFrom utils read.csv
read_call_graph <- function(edges, metrics,
                            from_col = "from", to_col = "to",
                            name_col = "name", complexity_col = "cyclo") {

  # ---- read files if paths are given -------------------------------------------

  read_input <- function(x, what) {
    if (is.character(x) && length(x) == 1L) {
      if (!file.exists(x)) {
        stop("`", what, "` file not found: ", x, call. = FALSE)
      }
      x <- utils::read.csv(x, stringsAsFactors = FALSE)
    }
    if (!is.data.frame(x)) {
      stop("`", what, "` must be a data.frame or the path to a `.csv` file.",
           call. = FALSE)
    }
    x
  }

  edges   <- read_input(edges, "edges")
  metrics <- read_input(metrics, "metrics")

  # ---- validate columns ---------------------------------------------------------

  missing_edge_cols <- setdiff(c(from_col, to_col), names(edges))
  if (length(missing_edge_cols) > 0) {
    stop("`edges` is missing required columns: ",
         paste(missing_edge_cols, collapse = ", "), call. = FALSE)
  }

  missing_metric_cols <- setdiff(c(name_col, complexity_col), names(metrics))
  if (length(missing_metric_cols) > 0) {
    stop("`metrics` is missing required columns: ",
         paste(missing_metric_cols, collapse = ", "), call. = FALSE)
  }

  from <- as.character(edges[[from_col]])
  to   <- as.character(edges[[to_col]])
  nm   <- as.character(metrics[[name_col]])
  cc   <- metrics[[complexity_col]]

  # ---- validate content ----------------------------------------------------------

  if (anyNA(from) || anyNA(to)) {
    stop("`edges` contains missing caller or callee names.", call. = FALSE)
  }

  dup <- nm[duplicated(nm)]
  if (length(dup) > 0) {
    stop("`metrics` contains duplicated function names: ",
         paste(unique(dup), collapse = ", "), call. = FALSE)
  }

  if (!is.numeric(cc) || anyNA(cc)) {
    stop("`", complexity_col, "` must be numeric with no missing values.",
         call. = FALSE)
  }
  if (any(cc < 1)) {
    stop("`", complexity_col, "` must be >= 1 (McCabe complexity of a ",
         "function with no decision points is 1).", call. = FALSE)
  }

  missing_nodes <- setdiff(union(from, to), nm)
  if (length(missing_nodes) > 0) {
    stop("Functions present in `edges` but missing from `metrics`: ",
         paste(missing_nodes, collapse = ", "), call. = FALSE)
  }

  n_self <- sum(from == to)
  if (n_self > 0) {
    message(n_self, " self-loop(s) (recursive calls) kept in the graph.")
  }

  nodes_tbl <- tibble::tibble(name = nm)
  nodes_tbl[[complexity_col]] <- as.numeric(cc)

  edges_tbl <- tibble::tibble(from = from, to = to)

  tidygraph::tbl_graph(nodes = nodes_tbl, edges = edges_tbl, directed = TRUE)
}
