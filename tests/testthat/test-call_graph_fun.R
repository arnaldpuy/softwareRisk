make_fixture_dir <- function() {
  td <- file.path(tempdir(), paste0("cg_test_", as.integer(runif(1, 1, 1e8))))
  dir.create(td, showWarnings = FALSE)
  writeLines(c(
    "load_data <- function(x) x",
    "clean_data <- function(x) load_data(x)",
    "calc_scores <- function(x) if (length(x) > 0) mean(x) else 0",
    "helpers = function(x) vapply(x, calc_scores, numeric(1))",
    "compute_risk <- function(x, trim = if (TRUE) 1 else 2) calc_scores(clean_data(x))",
    "recurse <- function(n) if (n > 0) recurse(n - 1) else 0"
  ), file.path(td, "model.R"))
  td
}

test_that("call_graph_fun builds the expected graph from a directory", {
  td <- make_fixture_dir()
  on.exit(unlink(td, recursive = TRUE))

  g <- suppressMessages(call_graph_fun(dir = td))
  expect_s3_class(g, "tbl_graph")

  ig <- tidygraph::as.igraph(g)
  nm <- igraph::V(ig)$name
  expect_setequal(nm, c("load_data", "clean_data", "calc_scores",
                        "helpers", "compute_risk", "recurse"))

  edges <- igraph::as_data_frame(ig, what = "edges")
  # direct calls
  expect_true(any(edges$from == "clean_data" & edges$to == "load_data"))
  expect_true(any(edges$from == "compute_risk" & edges$to == "calc_scores"))
  expect_true(any(edges$from == "compute_risk" & edges$to == "clean_data"))
  # function passed as value to vapply() is treated as a call
  expect_true(any(edges$from == "helpers" & edges$to == "calc_scores"))
  # self-recursion dropped by default
  expect_false(any(edges$from == edges$to))
})

test_that("call_graph_fun computes cyclomatic complexity correctly", {
  td <- make_fixture_dir()
  on.exit(unlink(td, recursive = TRUE))

  g <- suppressMessages(call_graph_fun(dir = td))
  ig <- tidygraph::as.igraph(g)
  cyclo <- stats::setNames(igraph::vertex_attr(ig, "cyclo"),
                           igraph::V(ig)$name)

  expect_equal(unname(cyclo["load_data"]), 1)     # no decision points
  expect_equal(unname(cyclo["calc_scores"]), 2)   # one if
  expect_equal(unname(cyclo["compute_risk"]), 2)  # if in default argument
  expect_equal(unname(cyclo["recurse"]), 2)       # one if
})

test_that("call_graph_fun keeps self-loops and applies exclude on request", {
  td <- make_fixture_dir()
  on.exit(unlink(td, recursive = TRUE))

  g <- suppressMessages(call_graph_fun(dir = td, keep_self_loops = TRUE))
  edges <- igraph::as_data_frame(tidygraph::as.igraph(g), what = "edges")
  expect_true(any(edges$from == "recurse" & edges$to == "recurse"))

  g2 <- suppressMessages(call_graph_fun(dir = td, exclude = "recurse"))
  nm <- igraph::V(tidygraph::as.igraph(g2))$name
  expect_false("recurse" %in% nm)
})

test_that("call_graph_fun output feeds directly into all_paths_fun", {
  td <- make_fixture_dir()
  on.exit(unlink(td, recursive = TRUE))

  g <- suppressMessages(call_graph_fun(dir = td))
  out <- all_paths_fun(g)
  expect_true(nrow(out$paths) > 0)
  expect_true(all(is.finite(out$nodes$risk_score)))
})

test_that("call_graph_fun works on an installed package", {
  g <- suppressMessages(call_graph_fun(pkg = "softwareRisk"))
  nm <- igraph::V(tidygraph::as.igraph(g))$name
  expect_true(all(c("all_paths_fun", "uncertainty_fun", "gini_index_fun") %in% nm))

  edges <- igraph::as_data_frame(tidygraph::as.igraph(g), what = "edges")
  # all_paths_fun calls gini_index_fun and slope_fun
  expect_true(any(edges$from == "all_paths_fun" & edges$to == "gini_index_fun"))
  expect_true(any(edges$from == "all_paths_fun" & edges$to == "slope_fun"))
})

test_that("call_graph_fun validates its input", {
  expect_error(call_graph_fun(), "exactly one")
  expect_error(call_graph_fun(pkg = "stats", dir = "."), "exactly one")
  expect_error(call_graph_fun(dir = file.path(tempdir(), "does_not_exist_xyz")),
               "existing directory")
  expect_error(call_graph_fun(pkg = "not_an_installed_package_xyz"),
               "not installed")

  td <- file.path(tempdir(), "cg_empty_dir")
  dir.create(td, showWarnings = FALSE)
  on.exit(unlink(td, recursive = TRUE))
  expect_error(call_graph_fun(dir = td), "No `.R` files")
})

test_that("read_call_graph builds a valid graph from data.frames", {
  calls_df <- data.frame(from = c("a", "b"), to = c("b", "c"))
  cyclo_df <- data.frame(name = c("a", "b", "c"), cyclo = c(2, 5, 1))

  g <- read_call_graph(calls_df, cyclo_df)
  expect_s3_class(g, "tbl_graph")

  ig <- tidygraph::as.igraph(g)
  expect_setequal(igraph::V(ig)$name, c("a", "b", "c"))
  expect_equal(igraph::vertex_attr(ig, "cyclo"),
               c(2, 5, 1))

  out <- all_paths_fun(g)
  expect_equal(nrow(out$paths), 1L)  # single path a -> b -> c
})

test_that("read_call_graph reads csv files and honours custom column names", {
  td <- file.path(tempdir(), "cg_csv_test")
  dir.create(td, showWarnings = FALSE)
  on.exit(unlink(td, recursive = TRUE))

  e_csv <- file.path(td, "edges.csv")
  m_csv <- file.path(td, "metrics.csv")
  write.csv(data.frame(caller = "a", callee = "b"), e_csv, row.names = FALSE)
  write.csv(data.frame(fn = c("a", "b"), complexity = c(3, 4)), m_csv,
            row.names = FALSE)

  g <- read_call_graph(e_csv, m_csv,
                       from_col = "caller", to_col = "callee",
                       name_col = "fn", complexity_col = "complexity")
  ig <- tidygraph::as.igraph(g)
  expect_setequal(igraph::V(ig)$name, c("a", "b"))
  expect_equal(igraph::vertex_attr(ig, "complexity"), c(3, 4))
})

test_that("read_call_graph rejects malformed input", {
  calls_df <- data.frame(from = c("a", "b"), to = c("b", "c"))
  cyclo_df <- data.frame(name = c("a", "b", "c"), cyclo = c(2, 5, 1))

  expect_error(read_call_graph(data.frame(x = 1), cyclo_df),
               "missing required columns")
  expect_error(read_call_graph(calls_df, data.frame(x = 1)),
               "missing required columns")
  expect_error(read_call_graph(data.frame(from = "a", to = "zzz"), cyclo_df),
               "missing from `metrics`")
  expect_error(
    read_call_graph(calls_df,
                    data.frame(name = c("a", "a", "b", "c"),
                               cyclo = c(1, 1, 1, 1))),
    "duplicated"
  )
  expect_error(
    read_call_graph(calls_df,
                    data.frame(name = c("a", "b", "c"), cyclo = c(0, 1, 1))),
    ">= 1"
  )
  expect_error(
    read_call_graph(calls_df,
                    data.frame(name = c("a", "b", "c"),
                               cyclo = c("x", "y", "z"))),
    "numeric"
  )
  expect_error(read_call_graph("no_such_file.csv", cyclo_df), "not found")

  # self-loops are kept with a message
  expect_message(
    read_call_graph(data.frame(from = c("a", "a"), to = c("a", "b")),
                    data.frame(name = c("a", "b"), cyclo = c(1, 1))),
    "self-loop"
  )
})
