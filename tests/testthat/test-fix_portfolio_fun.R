test_that("fix_portfolio_fun greedily reduces total path risk", {
  data(synthetic_graph)
  out <- all_paths_fun(synthetic_graph)

  res <- fix_portfolio_fun(out, budget = 5)
  expect_named(res, c("portfolio", "plot"))
  expect_s3_class(res$portfolio, "tbl_df")
  expect_s3_class(res$plot, "ggplot")
  expect_equal(nrow(res$portfolio), 5L)

  # objective decreases monotonically and deltas are positive
  expect_true(all(diff(res$portfolio$objective_after) < 0))
  expect_true(all(res$portfolio$delta > 0))
  expect_true(all(diff(res$portfolio$cum_reduction) > 0))
  expect_equal(anyDuplicated(res$portfolio$node), 0L)

  # first step is verifiable: fixing the selected node and recomputing the
  # total path risk by hand must reproduce objective_after
  v <- res$portfolio$node[1]
  r_map <- stats::setNames(out$nodes$risk_score, out$nodes$name)
  total_after <- sum(vapply(out$paths$path_nodes, function(nodes_vec) {
    r <- r_map[nodes_vec]
    r[nodes_vec == v] <- 0
    1 - prod(1 - r)
  }, numeric(1)))
  expect_equal(res$portfolio$objective_after[1], total_after)

  # ...and it is the best single fix
  candidates <- unique(unlist(out$paths$path_nodes))
  totals <- vapply(candidates, function(u) {
    sum(vapply(out$paths$path_nodes, function(nodes_vec) {
      r <- r_map[nodes_vec]
      r[nodes_vec == u] <- 0
      1 - prod(1 - r)
    }, numeric(1)))
  }, numeric(1))
  expect_equal(total_after, min(totals))
})

test_that("fix_portfolio_fun supports the max objective", {
  data(synthetic_graph)
  out <- all_paths_fun(synthetic_graph)

  res <- suppressMessages(fix_portfolio_fun(out, budget = 3, objective = "max"))
  expect_true(all(res$portfolio$objective_after <=
                    max(out$paths$path_risk_score)))
  expect_true(all(diff(res$portfolio$objective_after) < 0))
})

test_that("fix_portfolio_fun stops early when the budget exceeds useful fixes", {
  # single path with two nodes: after fixing both, nothing improves
  g <- tidygraph::tbl_graph(
    nodes = data.frame(name = c("a", "b"), cyclo = c(5, 3)),
    edges = data.frame(from = 1, to = 2),
    directed = TRUE
  )
  out <- all_paths_fun(g)
  res <- suppressMessages(fix_portfolio_fun(out, budget = 10))
  expect_true(nrow(res$portfolio) <= 2L)
  # risk fully removed
  expect_equal(res$portfolio$objective_after[nrow(res$portfolio)], 0,
               tolerance = 1e-12)
})

test_that("fix_portfolio_fun validates its input", {
  data(synthetic_graph)
  out <- all_paths_fun(synthetic_graph)
  expect_error(fix_portfolio_fun(list(a = 1)), "all_paths_fun")
  expect_error(fix_portfolio_fun(out, budget = 0), "positive")
  expect_error(fix_portfolio_fun(out, objective = "banana"))

  g <- tidygraph::tbl_graph(
    nodes = data.frame(name = c("a", "b"), cyclo = c(1, 2)),
    edges = data.frame(from = c(1, 2), to = c(2, 1)),
    directed = TRUE
  )
  out_empty <- all_paths_fun(g)
  expect_error(fix_portfolio_fun(out_empty), "no paths")
})
