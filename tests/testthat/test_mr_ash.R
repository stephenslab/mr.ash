# Tests here aim at checking that the package "mr.ash" has the same outputs
# as "mr.ash.alpha". 

context("mr_ash")

# We want to test that mr.ash outputs the same result as mr.ash.alpha
test_that("check identical fit", {
  
  # Simulate X and y
  set.seed(1)
  n     <- 200
  p     <- 400
  pve   <- 0.2
  s     <- 10
  data  <- simulate_data(n, p, pve, s)
  
  # fitting mr.ash.alpha and convirting classes
  capture.output(fit.alpha <- mr.ash.alpha::mr.ash(data$X, data$y))
  fit.alpha$beta   <- drop(fit.alpha$beta)
  fit.alpha$pi     <- drop(fit.alpha$pi)
  fit.alpha$varobj <- -fit.alpha$varobj
  fit.alpha$sa2    <- fit.alpha$data$sa2
  fit.alpha$data   <- NULL
  # attr(fit.alpha, "class") <- "list"
  
  # fitting mr.ash
  capture.output(fit <- mr_ash(data$X, data$y))
  
  # Conform mr.ash output to mr.ash.alpha output
  fit$elbo     <- NULL
  fit$varobj   <- fit$progress$elbo
  fit$iter     <- max(fit$progress$iter)
  fit$progress <- NULL
  fit          <- fit[names(fit.alpha)]
  # attr(fit, "class") <- "list"
  
  # Skipping sa2 because of ongoing debate in slack channel
  #fit$data$sa2       <- NULL
  #fit.alpha$data$sa2 <- NULL
  
  # Check equivalence
  expect_equal(fit, fit.alpha, scale = 1,tolerance = 1e-8, 
               check.attributes = F, use.names = F)
})


# Test that mr.ash outputs the same phi, m, and s2 values as mr.ash.alpha
test_that("check identical posterior", {
  # Simulate X and y
  set.seed(1)
  n     <- 200
  p     <- 400
  pve   <- 0.2
  s     <- 10
  data  <- simulate_data(n, p, pve, s)
  
  # fitting mr.ash.alpha
  capture.output(fit.alpha <- mr.ash.alpha::mr.ash(data$X, data$y))
  post.alpha <- mr.ash.alpha::get.full.posterior(fit.alpha)
  
  # fitting mr.ash
  capture.output(fit <- mr_ash(data$X, data$y))
  post <- list(phi = fit$phi, m = fit$m, s2 = fit$s2)
  
  # Check same phi value
  expect_equal(post, post.alpha, scale = 1,tolerance = 1e-8)
})
