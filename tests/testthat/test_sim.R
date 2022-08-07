context("simulation")

test_that("pve from predicted data matches input pve", {

  set.seed(1)
  sim <- simulate_regression_data(
    n = 100, p = 100, s = 25, pve = .25, ncov = 5, intercept = 1
  )

  a = as.numeric(t(sim$beta) %*% cov(sim$X) %*% sim$beta)
  r = a / (a + 1)

  expect_equal(r, .25, scale = 1, tolerance = 1e-6)
})
