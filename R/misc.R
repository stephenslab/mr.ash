# Return TRUE if x is a finite numeric matrix with no missing entries.
is.numeric.matrix <- function (x)
  is.matrix(x) && is.numeric(x) && !anyNA(x) && all(is.finite(x))

# Return TRUE if x is a finite numeric vector with no missing entries.
is.numeric.vector <- function (x)
  is.numeric(x) && !anyNA(x) && all(is.finite(x))

# Return TRUE if x is a finite numeric scalar with no missing entries.
is.scalar <- function (x)
  is.numeric.vector(x) & length(x) == 1

# Remove covariate effects Regresses Z out from X and y; that is, X
# and y are projected into the space orthogonal to Z.
#'
#' @importFrom Matrix forceSymmetric
#'
remove_covariate <- function (X, y, Z, standardize = FALSE, intercept = TRUE) {

  # check if Z is null and intercept = FALSE
  if (is.null(Z) && !intercept && !standardize) {
    return(list(X = X, y = y, Z = Z,
                ZtZiZX = rep(0,dim(X)[2]), ZtZiZy = 0))
  } else if (is.null(Z) && !intercept) {

    # standardize but don't center
    return(list(X = scale(X, center = FALSE), y = y, Z = Z,
                ZtZiZX = rep(0,dim(X)[2]), ZtZiZy = 0))

  }

  # redefine y
  y = c(as.double(y))
  n = length(y)

  # add intercept if intercept = TRUE
  if (intercept) {
    if (is.null(Z))
      Z <- matrix(1,n,1)
    else
      Z <- cbind(1,Z)
  }

  if (ncol(Z) == 1) {
    ZtZ         = forceSymmetric(crossprod(Z))       # (Z^T Z) symmetric
    ZtZiZy      = as.vector(solve(ZtZ,c(y %*% Z)))   # (Z^T Z)^{-1} Z^T y
    ZtZiZX      = as.matrix(solve(ZtZ,t(Z) %*% X))   # (Z^T Z)^{-1} Z^T X
    X           = scale(X, center = intercept, scale = standardize)
    alpha       = mean(y)
    y           = y - alpha

  } else {
    ZtZ         = forceSymmetric(crossprod(Z))       # (Z^T Z) symmetric
    ZtZiZy      = as.vector(solve(ZtZ,c(y %*% Z)))   # (Z^T Z)^{-1} Z^T y
    ZtZiZX      = as.matrix(solve(ZtZ,t(Z) %*% X))   # (Z^T Z)^{-1} Z^T X
    y     = y - c(Z %*% ZtZiZy)
    X     = X - Z %*% ZtZiZX
  }

  return(list(X = X, y = y, Z = Z,
              ZtZiZX = ZtZiZX, ZtZiZy = ZtZiZy))
}

var.n <- function(x) {
  a <- x - mean(x)
  return (sum(a^2) / length(a))
}

# computelfsrmix function from Peter
computelfsrmix <- function (phi, mu, s) {
  # Inputs:
  # phi: p x K matrix of posterior assignment probabilities
  # mu:  p x K matrix of posterior means for each mixture component
  # s:   p x K matrix of posterior variances
  # Return: a vector of length p

  # Get the number of variables (p) and the number of mixture
  # components (k).
  p <- nrow(phi)
  k <- ncol(phi)

  # For each variable, get the posterior probability that the
  # regression coefficient is exactly zero.
  p0 <- phi[,1]

  # For each variable, get the posterior probability that the
  # regression coefficient is negative.
  if (k == 2) {
    pn <- phi[,2] * pnorm(0,mu[,2],sqrt(s[,2]))
  } else {
    pn <- rowSums(phi[,-1] * pnorm(0,mu[,-1],sqrt(s[,-1])))
  }

  # Compute the local false sign rate (LFSR) following the formula
  # given in the Biostatistics paper, "False discovery rates: a new
  # deal".
  lfsr     <- rep(0,p)
  b        <- pn > 0.5*(1 - p0)
  lfsr[b]  <- 1 - pn[b]
  lfsr[!b] <- p0[!b] + pn[!b]
  lfsr[lfsr < 0] <- 0 # avoid negatives

  return(lfsr)
}


# This function extracts mean, sd, and lfsr from fit object
get_posterior_summary <- function(fit) {
  posterior <- list()

  # Separate the posterior summary
  posterior$m    <- fit$m
  posterior$s2   <- fit$s2
  posterior$lfsr <- fit$lfsr

  # return as posterior object
  class(posterior) <- c("posterior", "list")
  return(posterior)
}
