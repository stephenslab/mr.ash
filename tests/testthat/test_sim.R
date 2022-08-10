context("simulation")

test_that("pve from simulated data matches input pve", {
  set.seed(1)
  exp_pve_vec <- c(0, 0.1, 0.5, 0.9)
  num_sims    <- length(exp_pve_vec)
  sim_pve_vec <- numeric(num_sims)
  for (i in 1:num_sims) {
    sim <- simulate_regression_data(
      n = 10000, p = 100, s = 10, pve = exp_pve_vec[i], intercept = -
    )
    sim_pve_vec[i] <- with(sim,var(drop(X %*% b))/var(y))
  }
  expect_equal(exp_pve_vec, sim_pve_vec, scale = 1, tolerance = 0.01)
})

test_that("pve from simulated data is 1 when pve is 1 and sigma is 0", {
  set.seed(1)
  sim <- simulate_regression_data(
    n = 100, p = 100, s = 10, pve = 1, sigma = 0, intercept = -1
  )
  r <- with(sim,var(drop(X %*% b))/var(y))
  expect_equal(r, 1, scale = 1, tolerance = 1e-4)
})

