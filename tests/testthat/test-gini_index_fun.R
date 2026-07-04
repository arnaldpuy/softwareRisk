test_that("gini_index_fun returns 0 for equal values", {
  expect_equal(gini_index_fun(c(5, 5, 5, 5)), 0)
})

test_that("gini_index_fun returns a positive value for unequal values", {
  g <- gini_index_fun(c(1, 2, 3, 4))
  expect_true(g > 0 && g < 1)
})

test_that("gini_index_fun returns 0 for length-1 input", {
  expect_equal(gini_index_fun(42), 0)
})

test_that("gini_index_fun removes non-finite values", {
  g <- gini_index_fun(c(NA, 1, 2, Inf, 3))
  expect_true(is.numeric(g) && length(g) == 1)
})

test_that("gini_index_fun returns 0 (not NaN) for all-zero input", {
  expect_equal(gini_index_fun(c(0, 0, 0)), 0)
})
