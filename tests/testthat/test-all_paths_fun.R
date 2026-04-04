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
