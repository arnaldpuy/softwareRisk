make_small_ua <- function() {
  data(synthetic_graph, envir = environment())
  ig <- tidygraph::as.igraph(synthetic_graph)
  sub_ig <- igraph::induced_subgraph(ig, 1:12)
  g <- tidygraph::as_tbl_graph(sub_ig)
  out <- all_paths_fun(g)
  uncertainty_fun(out, N = 2^4, order = "first")
}

test_that("sensitivity_plot_fun builds the aggregate and node-level views", {
  ua <- make_small_ua()

  p1 <- suppressMessages(sensitivity_plot_fun(ua))
  expect_s3_class(p1, "ggplot")
  b1 <- ggplot2::ggplot_build(p1)
  expect_true(length(b1$data) >= 1)

  p2 <- suppressMessages(
    sensitivity_plot_fun(ua, nodes = ua$nodes$name[1:3])
  )
  expect_s3_class(p2, "ggplot")
  b2 <- ggplot2::ggplot_build(p2)
  # faceted by node
  expect_equal(length(unique(b2$layout$layout$PANEL)), 3L)
})

test_that("sensitivity_plot_fun filters by index type", {
  ua <- make_small_ua()
  p <- suppressMessages(sensitivity_plot_fun(ua, index = "Ti"))
  b <- ggplot2::ggplot_build(p)
  # a single fill group remains
  expect_equal(length(unique(b$data[[1]]$fill)), 1L)
})

test_that("sensitivity_plot_fun validates its input", {
  ua <- make_small_ua()
  expect_error(sensitivity_plot_fun(list(a = 1)), "uncertainty_fun")
  expect_error(sensitivity_plot_fun(ua, nodes = "not_a_node"), "not found")
  expect_error(sensitivity_plot_fun(ua, nodes = character(0)), "character")
  expect_error(sensitivity_plot_fun(ua, index = "banana"))

  bad <- ua
  bad$nodes$sensitivity_analysis <- NULL
  expect_error(sensitivity_plot_fun(bad), "sensitivity_analysis")
})
