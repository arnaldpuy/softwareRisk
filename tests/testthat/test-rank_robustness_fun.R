# hand-built uncertainty output with three paths and three draws
make_fake_ua <- function() {
  list(
    nodes = tibble::tibble(
      name = c("a", "b"),
      uncertainty_analysis = list(c(0.9, 0.2, 0.9), c(0.1, 0.8, 0.1)),
      sensitivity_analysis = list(NULL, NULL)
    ),
    paths = tibble::tibble(
      path_id = 1:3,
      path_str = c("p1", "p2", "p3"),
      hops = c(1L, 1L, 1L),
      uncertainty_analysis = list(
        c(0.9, 0.1, 0.9),   # top in draws 1 and 3
        c(0.5, 0.9, 0.5),   # top in draw 2
        c(0.1, 0.5, 0.1)    # never top
      ),
      gini_index = list(0, 0, 0),
      risk_trend = list(0, 0, 0)
    )
  )
}

test_that("rank_robustness_fun computes exact top-k probabilities", {
  rr <- rank_robustness_fun(make_fake_ua(), top_k = 1, what = "paths")

  expect_named(rr, c("summary", "consensus_correlation", "top_k", "what",
                     "n_draws"))
  expect_equal(rr$n_draws, 3L)

  s <- rr$summary
  prob <- stats::setNames(s$top_k_prob, s$id)
  expect_equal(unname(prob[c("1", "2", "3")]), c(2 / 3, 1 / 3, 0))

  # path 3 is always ranked last
  expect_equal(s$rank_median[s$id == 3], 3)
  # sorted by decreasing probability
  expect_true(all(diff(s$top_k_prob) <= 0))
  expect_true("path_str" %in% names(s))
})

test_that("rank_robustness_fun works on nodes", {
  rr <- rank_robustness_fun(make_fake_ua(), top_k = 1, what = "nodes")
  s <- rr$summary
  prob <- stats::setNames(s$top_k_prob, s$id)
  expect_equal(unname(prob[c("a", "b")]), c(2 / 3, 1 / 3))
})

test_that("rank_robustness_fun clamps top_k with a message", {
  expect_message(rr <- rank_robustness_fun(make_fake_ua(), top_k = 100),
                 "using top_k = 3")
  expect_equal(rr$top_k, 3L)
  expect_equal(rr$summary$top_k_prob, rep(1, 3))
})

test_that("rank_robustness_fun treats NA draws as never in the top-k", {
  fake <- make_fake_ua()
  fake$paths$uncertainty_analysis[[1]] <- c(NA_real_, NA_real_, NA_real_)
  rr <- rank_robustness_fun(fake, top_k = 1)
  expect_equal(rr$summary$top_k_prob[rr$summary$id == 1], 0)
})

test_that("rank_robustness_fun on real uncertainty_fun output is coherent", {
  data(synthetic_graph)
  out <- all_paths_fun(synthetic_graph)
  ua <- uncertainty_fun(out, N = 2^4, order = "first")

  rr <- rank_robustness_fun(ua, top_k = 10)
  expect_equal(nrow(rr$summary), nrow(out$paths))
  # on average exactly top_k items are in the top-k per draw
  expect_equal(sum(rr$summary$top_k_prob), 10, tolerance = 0.5)
  expect_true(rr$consensus_correlation > 0 && rr$consensus_correlation <= 1)
})

test_that("rank_robustness_fun validates its input", {
  expect_error(rank_robustness_fun(list(a = 1)), "uncertainty_fun")
  expect_error(rank_robustness_fun(make_fake_ua(), top_k = 0), "positive")

  fake <- make_fake_ua()
  fake$paths$uncertainty_analysis[[1]] <- c(0.1, 0.2)  # wrong length
  expect_error(rank_robustness_fun(fake), "same length")

  fake2 <- make_fake_ua()
  fake2$paths <- fake2$paths[0, ]
  expect_error(rank_robustness_fun(fake2), "no rows")
})

test_that("rank_robustness_plot returns a ggplot", {
  rr <- rank_robustness_fun(make_fake_ua(), top_k = 1)
  p <- rank_robustness_plot(rr, top_n = 2)
  expect_s3_class(p, "ggplot")
  built <- ggplot2::ggplot_build(p)
  expect_equal(nrow(built$data[[2]]), 2L)  # top_n points

  expect_error(rank_robustness_plot(list(a = 1)), "rank_robustness_fun")
  expect_error(rank_robustness_plot(rr, top_n = 0), "positive")
})
