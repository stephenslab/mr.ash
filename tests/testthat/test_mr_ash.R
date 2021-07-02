# We want to test that mr.ash outputs the same result as mr.ash.alpha
# Not sure if we can directly compare mr.ash output objects? Or do we
# want to compare each component (get.full.posterior results?)

library(mr.ash.alpha)
context("Equivalence to mr.ash.alpha")

test_that("equal phi values", {
  f <- list.files(path = "/project2/mstephens/lwang19/mr.ash/data/", pattern = "\\.rds$", full.names = T)
  fs <- lapply(f, readRDS)
  
  for (i in 1:length(fs)) {
    X <- fs[[i]]$X
    y <- fs[[i]]$y
    
    fit.alpha <- mr.ash.alpha::mr.ash(X, y)
    post.alpha <- mr.ash.alpha::get.full.posterior(fit.alpha)
    phi.alpha <- post.alpha$phi
    
    fit <- mr.ash::mr.ash(X, y)
    post <- mr.ash::get.full.posterior(fit)
    phi <- post$phi
    
    expect_equal(phi, phi.alpha) # expect_equal compare numerical equivalence
    #expect_identical(fit, fit.alpha) 
    # expect_identical compares the entirety of two R objects
  }
})