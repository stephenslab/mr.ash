context("input")

# Testing that multiplying X by a constant does not affect output phi
test_that("Scaling X result in same phi", {
  # Simulate X and y
  set.seed(1)
  n        <- 200
  p        <- 400
  pve      <- 0.2
  s        <- 10
  data     <- simulate_data(n, p, pve, s)
  epstol   <- 1e-12
  X.scaled <- data$X *  runif(1, -10, 10) + epstol
  
  # fit mr.ash (X, y), (X.scaled, y), (X, y.scaled), (X.scaled, y.scaled)
  capture.output(fit.Xy   <- mr_ash(data$X, data$y))
  capture.output(fit.Xsy  <- mr_ash(X.scaled, data$y))

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
  data     <- simulate_data(n, p, pve, s)
  X.std    <- scale(data$X, center = TRUE, scale = TRUE)
  
  # fit mr.ash (X, y), (X.std, y)
  capture.output(fit.Xsy <- mr_ash_test(X.std, data$y, standardize = T))
  capture.output(fit.Xy <- mr_ash_test(X.std, data$y))
  
  # Check that phi values are invariant wrt standardize flag
  expect_equal(fit.Xy$phi,fit.Xsy$phi,scale = 1,tolerance = 1e-8)
  # range(abs(fit.Xy$phi - fit.Xsy$phi))
})
