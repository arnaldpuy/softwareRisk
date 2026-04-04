test_that("uncertainty_fun returns expected structure", {
  data(synthetic_graph)

  # Use a small subgraph: first 10 nodes
  ig <- tidygraph::as.igraph(synthetic_graph)
  sub_ids <- 1:10
  sub_ig <- igraph::induced_subgraph(ig, sub_ids)
  sub_graph <- tidygraph::as_tbl_graph(sub_ig)

  out <- all_paths_fun(sub_graph, p = 1)

  skip_if(nrow(out$paths) == 0, "No paths in subgraph")

  results <- uncertainty_fun(all_paths_out = out, N = 2^2, order = "first")

  expect_type(results, "list")
  expect_named(results, c("nodes", "paths"))

  expect_true("name" %in% names(results$nodes))
  expect_true("uncertainty_analysis" %in% names(results$nodes))
  expect_true("sensitivity_analysis" %in% names(results$nodes))

  expect_true("path_id" %in% names(results$paths))
  expect_true("uncertainty_analysis" %in% names(results$paths))
  expect_true("gini_index" %in% names(results$paths))
  expect_true("risk_trend" %in% names(results$paths))

  # SA structure: each element is a sensobol sobol_indices list with a
  # $results data frame covering all four raw Sobol parameters
  si <- results$nodes$sensitivity_analysis[[1]]
  expect_true(is.list(si))
  expect_true("results" %in% names(si))
  expected_params <- c("alpha", "beta", "gamma", "p")
  expect_true(all(expected_params %in% unique(si$results$parameters)))

  # UA draws: length N for nodes and paths
  N <- 2^2
  expect_length(results$nodes$uncertainty_analysis[[1]], N)
  first_path_ua <- results$paths$uncertainty_analysis[[1]]
  expect_true(is.numeric(first_path_ua))
  expect_length(first_path_ua, N)

  # Sobol indices should be finite and within [-0.1, 1.1] (allow small numerics)
  si_vals <- si$results$original
  expect_true(all(is.finite(si_vals)))
  expect_true(all(si_vals >= -0.1 & si_vals <= 1.1))
})
