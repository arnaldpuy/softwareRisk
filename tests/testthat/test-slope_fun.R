test_that("slope_fun returns correct slope for a perfect linear sequence", {
  expect_equal(slope_fun(c(1, 2, 3, 4)), 1)
})

test_that("slope_fun returns correct slope for a decreasing sequence", {
  expect_equal(slope_fun(c(4, 3, 2, 1)), -1)
})

test_that("slope_fun returns 0 for a constant sequence", {
  expect_equal(slope_fun(c(5, 5, 5, 5)), 0)
})

test_that("slope_fun returns 0 for length-1 input", {
  expect_equal(slope_fun(42), 0)
})

test_that("slope_fun returns 0 for empty input", {
  expect_equal(slope_fun(numeric(0)), 0)
})

test_that("slope_fun removes non-finite values but keeps original positions", {
  # finite values 1, 2, 3 sit at positions 2, 3, 5: OLS slope on those
  # positions is 9/14, not 1 (which a compressed index would give)
  expect_equal(slope_fun(c(NA, 1, 2, Inf, 3)), 9 / 14)
})

test_that("slope_fun does not distort the trend when values are dropped", {
  # values at positions 1 and 3 rise by 1 per position
  expect_equal(slope_fun(c(1, NA, 3)), 1)
})
