context("mr_ash")

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
  capture.output(fit1 <- mr_ash(X,y,control = list(convtol = 1e-14)))
  capture.output(fit2 <- mr_ash(X,y,beta.init = fit1$beta,pi = fit1$pi,
                                sa2 = fit1$sa2,sigma2 = fit1$sigma2,
                                control = list(convtol = 1e-14)))

  # The mr.ash model fits, aside from a few bookkeeping details,
  # should be almost the same.
  fit1$progress <- NULL
  fit2$progress <- NULL
  expect_equal(fit1,fit2,scale = 1,tolerance = 1e-8)
})

# We want to test that mr.ash outputs the same result as mr.ash.alpha
test_that("equal phi values", {
  # skip("Leah is working on this test")
  
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
  attr(fit.alpha, "class") <- "list"
  
  # fitting mr.ash
  capture.output(fit <- mr_ash(data$X, data$y))
  
  # Conform mr.ash output to mr.ash.alpha output
  fit$elbo     <- NULL
  fit$varobj   <- fit$progress$elbo
  fit$iter     <- max(fit$progress$iter)
  fit$progress <- NULL
  fit          <- fit[names(fit.alpha)]
  attr(fit, "class") <- "list"
  
  # Skipping sa2 because of ongoing debate in slack channel
  fit$data$sa2       <- NULL
  fit.alpha$data$sa2 <- NULL
  
  # Check equivalence
  expect_equal(fit, fit.alpha, scale = 1,tolerance = 1e-8)
})


# ELBO output non-decreasing
test_that("non-decreasing ELBO", {
  
  # Simulate X and y
  set.seed(1)
  n     <- 200
  p     <- 400
  pve   <- 0.2
  s     <- 10
  data  <- simulate_data(n, p, pve, s)
  
  # fit mr.ash
  capture.output(fit <- mr_ash(data$X, data$y))
  
  # ELBO should be non-decreasing
  expect_true(all(fit$elbo == cummax(fit$elbo)))
})
