# diamond graph: e -> m1 -> s and e -> m2 -> s, two paths
make_diamond_out <- function() {
  g <- tidygraph::tbl_graph(
    nodes = data.frame(name = c("e", "m1", "m2", "s"),
                       cyclo = c(2, 10, 4, 6)),
    edges = data.frame(from = c(1, 1, 2, 3), to = c(2, 3, 4, 4)),
    directed = TRUE
  )
  all_paths_fun(g)
}

test_that("node_exposure_fun computes participation and risk load correctly", {
  out <- make_diamond_out()
  expect_equal(nrow(out$paths), 2L)

  ne <- node_exposure_fun(out)
  expect_s3_class(ne, "tbl_df")
  expect_equal(nrow(ne), 4L)

  n_paths <- stats::setNames(ne$n_paths, ne$name)
  expect_equal(unname(n_paths[c("e", "s")]), c(2L, 2L))
  expect_equal(unname(n_paths[c("m1", "m2")]), c(1L, 1L))

  # nodes on all paths carry the full risk load
  total_risk <- sum(out$paths$path_risk_score)
  rl <- stats::setNames(ne$risk_load, ne$name)
  expect_equal(unname(rl["e"]), total_risk)
  expect_equal(unname(rl["s"]), total_risk)
  expect_equal(unname(rl["m1"] + rl["m2"]), total_risk)

  # shares
  expect_equal(stats::setNames(ne$path_share, ne$name)[["e"]], 1)
  expect_equal(stats::setNames(ne$risk_load_share, ne$name)[["e"]], 1)

  # mean path risk consistency
  expect_equal(ne$mean_path_risk, ne$risk_load / ne$n_paths)
})

test_that("node_exposure_fun is sorted by risk load with coherent ranks", {
  data(synthetic_graph)
  out <- all_paths_fun(synthetic_graph)
  ne <- node_exposure_fun(out)

  expect_true(all(diff(ne$risk_load) <= 0))
  expect_equal(ne$rank_risk_load, seq_len(nrow(ne)))
  expect_true(all(ne$n_paths >= 1))
  expect_true(all(ne$path_share > 0 & ne$path_share <= 1))
})

test_that("node_exposure_fun returns an empty tibble when there are no paths", {
  # two-node cycle: no entry or sink nodes
  g <- tidygraph::tbl_graph(
    nodes = data.frame(name = c("a", "b"), cyclo = c(1, 2)),
    edges = data.frame(from = c(1, 2), to = c(2, 1)),
    directed = TRUE
  )
  out <- all_paths_fun(g)
  ne <- node_exposure_fun(out)
  expect_equal(nrow(ne), 0L)
  expect_true(all(c("name", "risk_score", "n_paths", "risk_load",
                    "rank_risk_load") %in% names(ne)))
})

test_that("node_exposure_fun validates its input", {
  expect_error(node_exposure_fun(list(a = 1)), "all_paths_fun")
  out <- make_diamond_out()
  out$paths$path_nodes <- NULL
  expect_error(node_exposure_fun(out), "missing required columns")
})
