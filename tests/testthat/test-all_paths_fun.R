test_that("all_paths_fun returns expected structure on synthetic_graph", {
  data(synthetic_graph)
  out <- all_paths_fun(synthetic_graph, p = 1)

  expect_type(out, "list")
  expect_named(out, c("nodes", "paths"))

  expect_s3_class(out$nodes, "tbl_df")
  expect_s3_class(out$paths, "tbl_df")

  expected_node_cols <- c("name", "cyclomatic_complexity", "indeg", "outdeg", "btw", "risk_score")
  expect_true(all(expected_node_cols %in% names(out$nodes)))

  expected_path_cols <- c("path_id", "path_nodes", "path_str", "hops",
                          "path_risk_score", "path_cc", "gini_node_risk",
                          "risk_slope", "risk_mean", "risk_sum")
  expect_true(all(expected_path_cols %in% names(out$paths)))

  expect_true(nrow(out$nodes) > 0)
  expect_true(nrow(out$paths) > 0)
  expect_true(all(out$nodes$risk_score >= 0 & out$nodes$risk_score <= 1))
})

test_that("all_paths_fun validates p range", {
  data(synthetic_graph)
  expect_error(all_paths_fun(synthetic_graph, p = 5), "\\[-1, 2\\]")
  expect_error(all_paths_fun(synthetic_graph, p = -2), "\\[-1, 2\\]")
})

test_that("all_paths_fun returns finite risk scores with a zero weight and p < 0", {
  data(synthetic_graph)

  # regression test: 0 * (0^p) = 0 * Inf = NaN before the eps guard
  out <- all_paths_fun(synthetic_graph, alpha = 0.5, beta = 0.5, gamma = 0, p = -1)
  expect_true(all(is.finite(out$nodes$risk_score)))
  expect_true(all(out$nodes$risk_score >= 0 & out$nodes$risk_score <= 1))
  expect_true(all(is.finite(out$paths$path_risk_score)))

  out2 <- all_paths_fun(synthetic_graph, alpha = 0, beta = 0.9, gamma = 0.1, p = -0.5)
  expect_true(all(is.finite(out2$nodes$risk_score)))
})

test_that("all_paths_fun produces no NaN gini values with default settings", {
  data(synthetic_graph)
  out <- all_paths_fun(synthetic_graph, p = 1)
  expect_true(all(is.finite(out$paths$gini_node_risk)))
})
