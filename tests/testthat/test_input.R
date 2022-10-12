# Tests here check that the flags used in function mr_ash are functional,
# and that function is not sensitive to X scaling.
context("input")

# Testing that multiplying X by a constant does not affect output phi
test_that("Scaling X result in same phi", {
  # Simulate X and y
  set.seed(1)
  n        <- 200
  p        <- 400
  pve      <- 0.2
  s        <- 10
  data     <- simulate_regression_data(n = n, p = p, pve = pve, s = s)
  epstol   <- 1e-12
  X.scaled <- data$X * runif(1, -10, 10) + epstol

  # fit mr.ash (X, y), (X.scaled, y)
  fit0 <- init_mr_ash(data$X, data$y)
  capture.output(fit.Xy <- fit_mr_ash(data$X, data$y, fit0 = fit0,
                                      standardize = TRUE))
  capture.output(fit.Xsy <- fit_mr_ash(X.scaled, data$y, fit0 = fit0,
                                       standardize = TRUE))

  # Check that phi values are invariant wrt X scaling
  expect_equal(fit.Xy$phi,fit.Xsy$phi,scale = 1,tolerance = 1e-8)
})

# When data set is standardized, setting "standardize == TRUE"
# should have no effect on the output... right?
test_that("Standardized X yield the same result regardless of
          'standardize' flag", {

  # Simulate X and y
  set.seed(1)
  n        <- 200
  p        <- 400
  pve      <- 0.2
  s        <- 10
  data     <- simulate_regression_data(n = n, p = p, pve = pve, s = s)
  X.std    <- scale(data$X, center = TRUE, scale = TRUE)

  # fit mr.ash (X, y), (X.std, y)
  capture.output(fit.Xsy <- fit_mr_ash(X.std, data$y, standardize = T))
  capture.output(fit.Xy  <- fit_mr_ash(X.std, data$y))

  # Check that phi values are invariant wrt standardize flag
  expect_equal(fit.Xy$phi,fit.Xsy$phi,scale = 1,tolerance = 1e-6)
  # range(abs(fit.Xy$phi - fit.Xsy$phi))
})


# Check intercept flag: if the intercept flag is FALSE, fit$intercept == 0
# Although this test is kind of useless once we know function doesn't throw
# an error...
test_that("Check length of beta w/ and w/o intercept flat", {

  # Simulate X and y
  set.seed(1)
  n        <- 200
  p        <- 400
  pve      <- 0.2
  s        <- 10
  data     <- simulate_regression_data(n = n, p = p, pve = pve, s = s)

  # fit mr.ash (X, y), intercept == FALSE
  capture.output(fit.Xy <- fit_mr_ash(data$X, data$y, intercept = F))
  expect_equal(fit.Xy$intercept, 0)

})
