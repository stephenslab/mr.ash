# Fit a univariate linear regression model y ~ x*b + e, e ~ N(0,se),
# separately for each column of X, in which b ~ N(0,s0). Here, se and
# s0 are variances.
bayes_lr_ridge <- function (X, y, se = 1, s0 = 1) {

  # These two lines are not very memory efficient, and could be
  # improved.
  X  <- scale(X, center = TRUE, scale = FALSE)
  xx <- colSums(X^2)

  # Compute the least-squares estimate of b, and its variance.
  bhat <- drop(y %*% X)/xx
  shat <- se/xx
  
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

function () {
    w <- colSums(X^2)
    S   <- outer(1/w, sa2, '+') * sigma2
    Phi <- -b^2/S/2 - log(S)/2
    Phi <- exp(Phi - apply(Phi,1,max))
    Phi <- Phi / rowSums(Phi)
    pi  <- colMeans(Phi)
}
