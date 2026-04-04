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

test_that("slope_fun removes non-finite values", {
  expect_equal(slope_fun(c(NA, 1, 2, Inf, 3)), 1)
})
