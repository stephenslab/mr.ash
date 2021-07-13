context("mr_ash")

test_that("re-running mr.ash after solution has converged yields same fit",{
  set.seed(1)

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
                                sa2 = fit1$data$sa2,sigma2 = fit1$sigma2,
                                control = list(convtol = 1e-14)))

  # The mr.ash model fits, aside from a few bookkeeping details,
  # should be almost the same.
  fit1$progress <- NULL
  fit2$progress <- NULL
  expect_equal(fit1,fit2,scale = 1,tolerance = 1e-8)
})

# We want to test that mr.ash outputs the same result as mr.ash.alpha
# Not sure if we can directly compare mr.ash output objects? Or do we
# want to compare each component (get.full.posterior results?)
test_that("equal phi values", {
  skip("Leah is working on this test")
  f <- list.files(path = "/project2/mstephens/lwang19/mr.ash/data/", pattern = "\\.rds$", full.names = T)
  fs <- lapply(f, readRDS)
  
  for (i in 1:length(fs)) {
    X <- fs[[i]]$X
    y <- fs[[i]]$y
    
    fit.alpha <- mr.ash.alpha::mr.ash(X, y)
    fit.alpha$beta <- drop(fit.alpha$beta)
    fit.alpha$pi <- drop(fit.alpha$pi)
    fit <- mr_ash(X, y, verbose = "none")
    
    # expect_equal(fit$beta, fit.alpha$beta)
    expect_identical(fit, fit.alpha) 
    # expect_identical compares the entirety of two R objects
  }
})
