context("predict")

test_that("predict with newx = NULL matches predict with newx = X", {

  set.seed(1)
  n <- 1000
  p <- 10

  X <- matrix(
    data = rnorm(n = n * p), nrow = n, ncol = p
  )

  beta <- matrix(data = rnorm(10), nrow = 10, ncol = 1)

  y <- rnorm(n) + X %*% beta

  out_n <- suppressWarnings(
    fit_mr_ash(X, y, init_mr_ash(X, y, intercept = F), verbose = "none")
  )

  preds_new <- predict(out_n, X)
  preds_def <- predict(out_n)
  test <- preds_new - preds_def

  expect_equal(test, rep(0, n), scale = 1, tolerance = 1e-8)
})
