# Tests here intend to check the fit of function "mr_ash". Most of them are
# simple sanity checks.

context("fit")

test_that("re-running mr.ash after solution has converged yields same fit",{

  # Simulate a small regression data set with n = 200 samples and
  # p = 400 predictors.
  set.seed(1)
  n          <- 200
  p          <- 400
  X          <- matrix(rnorm(n*p),n,p)
  beta       <- double(p)
  beta[1:10] <- 1:10
  y          <- drop(X %*% beta + rnorm(n))

  # Fit the mr.ash model, then fit a second time in which the fit is
  # initialized to the estimates returned from first mr.ash call.
  capture.output(fit1 <- fit_mr_ash(X,y,control = list(convtol = 1e-14)))

  fit0_2 <- init_mr_ash(X, y, b = fit1$b,prior.weights = fit1$prior$weights,
                        prior.sd = fit1$prior$sd,resid.sd = fit1$resid.sd)
  capture.output(fit2 <- fit_mr_ash(X,y,fit0_2,
                                control = list(convtol = 1e-14)))

  # The mr.ash model fits, aside from a few bookkeeping details,
  # should be almost the same.
  fit1$progress <- NULL
  fit2$progress <- NULL
  expect_equal(fit1,fit2,scale = 1,tolerance = 1e-5)
})

# ELBO output non-decreasing
test_that("non-decreasing ELBO", {

  # Simulate X and y
  set.seed(1)
  n     <- 200
  p     <- 400
  pve   <- 0.2
  s     <- 10
  data  <- simulate_regression_data(n = n, p = p, pve = pve, s = s)

  # fit mr.ash
  capture.output(fit <- fit_mr_ash(data$X, data$y))

  # ELBO should be non-decreasing
  expect_true(all(fit$elbo == cummax(fit$elbo)))
})


# lfsr calculation should yield output between 0 and 1
test_that("lfsr between 0 and 1", {
  # Simulate X and y
  set.seed(2)
  n     <- 200
  p     <- 400
  pve   <- 0.2
  s     <- 10
  data  <- simulate_regression_data(n = n, p = p, pve = pve, s = s)

  # fit mr.ash
  capture.output(fit <- fit_mr_ash(data$X, data$y))
  test <- sum(fit$lfsr < 0 | fit$lfsr > 1)
  expect_equal(test, 0, scale = 1, tolerance = 1e-8)
})


# lfsr calculation should yield output between 0 and 1
test_that("computations scale appropriately with X and y", {
  # Simulate X and y
  set.seed(2)
  n     <- 200
  p     <- 400
  pve   <- 0.2
  s     <- 10
  data  <- simulate_regression_data(n = n, p = p, pve = pve, s = s)

  # fit mr.ash
  capture.output(fit <- fit_mr_ash(data$X, data$y))
  capture.output(fit2 <- fit_mr_ash(data$X*10, data$y*2))

  expect_equal(fit2$b/fit$b, rep(2/10,p)) # beta scales with y/X
  expect_equal(fit2$resid.var/fit$resid.var, 2^2) # scales with y^2
  expect_equal(fit2$prior$var/fit$prior$var,0.01) # with 1/X^2

  expect_equal(fit2$phi/fit$phi, 1) # this is very idealized! Probably won't happen
  expect_equal(fit2$m/fit$m, 2/10) # m scales with b
  expect_equal(fit2$s2/fit$s2, (2/10)^2) # scales with b^2

})


