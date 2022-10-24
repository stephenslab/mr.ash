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
  data  <- simulate_regression_data(n = n, p = p, pve = pve, s = s)

  w = colSums(data$X ^ 2)
  sa2 = (2^((0:19) / 20) - 1)^2
  sa2 = sa2 / median(w) * n
  K = length(sa2)
  pi = rep(1,K)/K

  # fitting mr.ash.alpha and convirting classes
  capture.output(fit.alpha <- mr.ash.alpha::mr.ash(
    data$X, data$y, sa2 = sa2, pi = pi)
  )
  fit.alpha$beta   <- drop(fit.alpha$beta)
  fit.alpha$pi     <- drop(fit.alpha$pi)
  fit.alpha$varobj <- -fit.alpha$varobj
  fit.alpha$sa2    <- fit.alpha$data$sa2
  fit.alpha$data   <- NULL
  # attr(fit.alpha, "class") <- "list"

  # fitting mr.ash
  capture.output(fit <- fit_mr_ash(
    data$X, data$y, init_mr_ash(
      data$X, data$y, init.method = "null", prior.sd = sa2, prior.weights = pi
      )
    )
  )

  # Conform mr.ash output to mr.ash.alpha output
  fit$elbo     <- NULL
  fit$varobj   <- fit$progress$elbo
  fit$iter     <- max(fit$progress$iter)
  fit$progress <- NULL
  fit$pi <- fit$posterior.weights
  fit$sa2 <- fit$prior$sd
  fit$sigma2 <- fit$resid.sd
  fit$beta <- fit$b
  fit$intercept <- fit$b_0
  fit          <- fit[names(fit.alpha)]
  # attr(fit, "class") <- "list"

  # Skipping sa2 because of ongoing debate in slack channel
  #fit$data$sa2       <- NULL
  #fit.alpha$data$sa2 <- NULL

  fit$beta <- drop(fit$beta)

  # Check equivalence
  #browser()
  expect_equal(fit, fit.alpha, tolerance = 1e-4,
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
  data  <- simulate_regression_data(n = n, p = p, pve = pve, s = s)

  w = colSums(data$X ^ 2)
  sa2 = (2^((0:19) / 20) - 1)^2
  sa2 = sa2 / median(w) * n
  K = length(sa2)
  pi = rep(1,K)/K

  # fitting mr.ash.alpha
  capture.output(fit.alpha <- mr.ash.alpha::mr.ash(
    data$X, data$y, sa2 = sa2, pi = pi)
  )
  post.alpha <- mr.ash.alpha::get.full.posterior(fit.alpha)

  # fitting mr.ash
  capture.output(fit <- fit_mr_ash(
    data$X, data$y, init_mr_ash(
      data$X, data$y, init.method = "null", prior.sd = sa2, prior.weights = pi
    )
  )
  )
  post <- list(phi = fit$phi, m = fit$m, s2 = fit$s2)

  # make dimnames match
  colnames(post.alpha$phi) <- colnames(post$phi)
  colnames(post.alpha$m) <- colnames(post$m)
  colnames(post.alpha$s2) <- colnames(post$s2)

  # Check same phi value
  expect_equal(post, post.alpha, scale = 1,tolerance = 1e-6)
})
