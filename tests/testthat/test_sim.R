context("simulation")

test_that("pve from simulated data matches input pve", {

  set.seed(1)

  exp_pve_vec <- c(0, .1, .5, .9)
  num_sims <- length(exp_pve_vec)
  sim_pve_vec <- numeric(num_sims)

  for (i in 1:num_sims) {

    sim <- simulate_regression_data(
      n = 100, p = 100, s = 25, pve = exp_pve_vec[i], ncov = 5, intercept = 1
    )


    a = drop(t(sim$b) %*% cov(sim$X) %*% sim$b)
    r = a / (a + 1)
    sim_pve_vec[i] <- r

  }

  expect_equal(exp_pve_vec, sim_pve_vec, scale = 1, tolerance = 1e-6)


})

test_that("pve from simulated data is 1 when pve is 1 and sigma is 0", {

  sim <- simulate_regression_data(
    n = 100, p = 100, s = 90, pve = 1, sigma = 0, ncov = 5, intercept = 1
  )


  a = drop(t(sim$b) %*% cov(sim$X) %*% sim$b)
  r = a / (a + 1)

  expect_equal(r, 1, scale = 1, tolerance = 1e-2)

})

