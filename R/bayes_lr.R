# For each variable (column of X), compute the least-squares estimate
# of b (bhat), and its variance (shat). If y is NULL, only the
# variance (shat) is returned.
simple_lr <- function (X, y = NULL, se = 1) {

  # The first two lines are not very memory efficient, and could be
  # improved.
  X    <- scale(X, center = TRUE, scale = FALSE)
  xx   <- colSums(X^2)
  y    <- y - mean(y)
  bhat <- drop(y %*% X)/xx
  shat <- se/xx
  return(list(bhat = bhat, shat = shat))
}

# Fit a univariate linear regression model y ~ x*b + e, e ~ N(0,se),
# separately for each column of X, in which b ~ N(0,s0). Here, se and
# s0 are *variances* (not standard deviations). X and y are assumed to
# be "centered" so that the mean of y and each column of X is zero,
# and xx = colSums(X^2). The summary statistics are bhat =
# sum(x*y)/sum(x*x) and shat = se/sum(x*x), where
bayes_lr_ridge <- function (bhat, shat, s0 = 1) {

  # Compute the posterior mean (mu1) and variance (s1) assuming a
  # normal prior with zero mean and variance s0.
  s1  <- s0/(1 + s0/shat)
  mu1 <- s1/shat * bhat

  # Compute the log-Bayes factor.
  lbf <- log(shat/(s0 + shat))/2 + mu1^2/(2*s1)

  # Return the least-squares estimate of b (bhat), its variance
  # (shat), the posterior mean (mu1) and variance (s1), and the
  # log-Bayes factor (lbf).
  return(list(bhat = bhat,
              shat = shat,
              mu1  = mu1,
              s1   = s1,
              lbf  = lbf))
}
