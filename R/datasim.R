#' @title Simulate Regression Data
#'
#' @description Simulate regression data from a mr ash model.
#'
#' @details Data are generated from the model \deqn{y | X, \beta, \sigma^2, \beta_0 ~
#' N(\beta_0 + X \beta + Zu, \sigma I_n),} where the regression coefficients
#' \eqn{\beta} are first generated as \deqn{\beta | s ~ ((p - s) / p) \delta_0 +
#' (s / p) N(0, 1),} and then scaled to achieve the desired proportion
#' of variance in y explained by X (\code{pve}).
#'
#' @param n number of samples. Must be at least 2.
#'
#' @param p number of regression variables. Must be at least 2.
#'
#' @param s number of non-zero effects
#'
#' @param sigma (optional) standard deviation of the residual. Defaults to 1.
#'
#' @param pve (optional) proportion of variance in Y explained by X. Defaults to
#' .5.
#'
#' @param standardize_X (optional) whether or not to standardize the columns of
#' X. Defaults to \code{TRUE}.
#'
#' @param intercept (optional) intercept of regression equation. Defaults to 0.
#' \eqn{\beta_0} below.
#'
#' @param ncov (optional) number of covariates in Z. Defaults to 0
#'
#' @return A list with components
#' \itemize{
#'   \item y - \code{n} vector containing regression response.
#'   \item X - \code{n} by \code{p} matrix of regression variables X.
#'   \item Z - \code{n} by \code{ncov} matrix of covariates Z.
#'   \item u - \code{ncov} vector containing covariate coefficients u.
#'   \item beta - \code{p} vector contained regression coefficients \eqn{\beta}.
#' }
#'
#' @importFrom stats rnorm sd cor
#'
#' @export
#'
#' @examples
#'
#' # simulate 100 x 100 regression data with 75% sparsity
#' simulate_regression_data(100, 100, 25)
#'
simulate_regression_data <- function (
  n,
  p,
  s,
  sigma = 1,
  pve = .5,
  standardize_X = TRUE,
  intercept = 0,
  ncov = 0
) {

  # Argument checking
  if (!(is.scalar(n) && n >= 2 && is.int(n)))
    stop("Input argument \"n\" should be an integer equal to 2 or more")
  if (!(is.scalar(p) && p >= 2 && is.int(p)))
    stop("Input argument \"p\" should be an integer equal to 2 or more")
  if (!(is.scalar(s) && s >= 1 && s <= p && is.int(s)))
    stop("Input argument \"s\" should be an integer between 1 and p (inclusive)")
  if (!(is.scalar(sigma) && sigma >= 0))
    stop("Input argument \"sigma\" should be a scalar greater than 0")
  if (!(is.scalar(pve) && pve >= 0 && pve < 1))
    stop("Input argument \"pve\" should be a scalar between 0 (inclusive) and 1 (exclusive)")
  if (!is.scalar(intercept))
    stop("Input argument \"intercept\" should be a scalar")
  if (!(is.scalar(ncov) && ncov >= 0 && is.int(ncov)))
    stop("Input argument \"pve\" should be a scalar between 0 and 1 (inclusive)")


  # simulate data matrix X:
  X <- matrix(rnorm(n * p), n, p)
  if (standardize_X)
    X <- scale(X)

  # check if s <= p
  if (s > p) {
    stop("number of effects more than variables")
  }

  # generate effect variables
  beta.idx = sample(p, s)
  beta = rep(0, p)

  # generate effects
  if(s > 0){
    beta.values = rnorm(s)
    beta[beta.idx] = beta.values
  }

  # Adjust the effects so that we control for the proportion of variance
  # explained (pve). That is, we adjust beta so that r = a/(a+1), where we
  # define a = beta'*cov(X)*beta. Here, sb is the variance of the (nonzero)
  # effects.
  sb   <- pve/(1-pve)/var(c(X %*% beta))
  beta <- sqrt(sb * sigma) * beta

  # Generate the covariate data (Z), and the linear effects of the
  # covariates (u).
  if (ncov > 0) {
    Z <- matrix(
      data = rnorm(n * ncov),
      nrow = n,
      ncol = ncov
    )
    u <- rnorm(ncov)
  } else {
    Z <- NULL
    u <- NULL
  }

  # Generate the quantitative trait measurements.
  y <- intercept + X %*% beta + sqrt(sigma)*rnorm(n)
  if (ncov > 0)
    y <- y + Z %*% u
  y <- c(y)

  # give names to outputs
  rnames <- paste0("s", 1:n)
  cnames_X <- paste0("v", 1:p)
  names(y) <- rnames
  rownames(X) <- rnames
  colnames(X) <- cnames_X
  names(beta) <- cnames_X

  if (ncov > 0) {

    cnames_Z <- paste0("c", 1:ncov)
    names(u) <- cnames_Z
    colnames(Z) <- cnames_Z
    rownames(Z) <- rnames

  }

  return(
    list(
      y = y, X = X, Z = Z, u = u, beta = beta
    )
  )

}
