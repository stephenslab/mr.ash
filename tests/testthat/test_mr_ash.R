context("mr_ash")

test_that("re-running mr.ash after solution has converged yields same result",{
  set.seed(1)
  
}))

# We want to test that mr.ash outputs the same result as mr.ash.alpha
# Not sure if we can directly compare mr.ash output objects? Or do we
# want to compare each component (get.full.posterior results?)
test_that("equal phi values", {
  skip()
  f <- list.files(path = "/project2/mstephens/lwang19/mr.ash/data/", pattern = "\\.rds$", full.names = T)
  fs <- lapply(f, readRDS)
  
  for (i in 1:length(fs)) {
    X <- fs[[i]]$X
    y <- fs[[i]]$y
    
    fit.alpha <- mr.ash.alpha::mr.ash(X, y)
    fit.alpha$beta <- drop(fit.alpha$beta)
    fit.alpha$pi <- drop(fit.alpha$pi)
    fit <- mr.ash::mr.ash(X, y, verbose = "none")
    
    # expect_equal(fit$beta, fit.alpha$beta)
    expect_identical(fit, fit.alpha) 
    # expect_identical compares the entirety of two R objects
  }
})
