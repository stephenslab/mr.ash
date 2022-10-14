#' @title Simulate Data from Multiple Regression Model
#'
#' @description Simulate data from a multiple regression model.
#'
#' @details Data are generated from the linear regressionmodel \deqn{y
#' \sim N(b_0 + Xb + Zu, \sigma^2I),} in which the non-zero elements
#' of vector \eqn{b} are drawn from the standard normal, and scaled to
#' achieve the desired proportion of variance in y explained by X (see
#' input argument \dQuote{pve}).
#'
#' @param n The number of samples. Should be 2 or greater.
#'
#' @param p The number of regression variables. Should be 2 or greater.
#'
#' @param s A non-negative number specifying the number of non-zero
#'   regression coefficients.
#'
#' @param sigma The standard deviation of the residual.
#'
#' @param pve Assuming \code{ncov = 0}, this is the proportion of
#'   variance in y explained by X. When \code{ncov > 0}, this is the
#'   proportion of variance in y explained by X after removing the
#'   linear effects of the covariates Z.
#'
#' @param center_X If \code{center_X = TRUE}, the columns of X are
#'   centered before simulating y so that each column as a mean of
#'   zero. Note that if \code{standardize_X = TRUE} this argument must
#'   also be set to \code{TRUE}.
#'
#' @param standardize_X If \code{standardize_X = TRUE}, the columns of
#'   X are standardized before simulating y so that each column has a
#'   mean of zero and a standard deviation of 1.
#'
#' @param intercept The intercept \eqn{b_0} in the regression model.
#'
#' @param ncov The number of covariates (columns of Z).
#'
#' @return A list with the following components:
#'
#' \item{y}{Vector of length n containing the regression responses.}
#'
#' \item{X}{\code{n} by \code{p} matrix of regression variables X.}
#'
#' \item{Z}{\code{n} by \code{ncov} matrix of covariates Z.}
#'
#' \item{u}{Vector of length \code{ncov} the containing covariate
#'   coefficients, u.}
#'
#' \item{b}{Vector of length \code{p} containing the regression
#'   coefficients \eqn{b}.}
#'
#' \item{intercept}{The intercept in the regression model.}
#'
#' \item{sigma}{The standard deviation of the residual.}
#'
#' @importFrom stats rnorm sd cor
#'
#' @export
#'
#' @examples
#' set.seed(1)
#' 
#' # Simulate 100 x 100 regression data with 75% sparsity.
#' simulate_regression_data(n = 100, p = 100, s = 25)
#'
#' # simulate 100 x 100 regression data with # 50% sparsity and 100%
#' # variance explained.
#' simulate_regression_data(n = 100, p = 100, s = 50, pve = 1, sigma = 0)
#'
#' # Simulate 100 x 100 regression data with 50% sparsity and 10
#' # covariates.
#' simulate_regression_data(n = 100, p = 100, s = 50, ncov = 10)
#'
simulate_regression_data <- function (
  n,
  p,
  s,
  sigma = 1,
  pve = 0.5,
  center_X = TRUE,
  standardize_X = FALSE,
  intercept = 0,
  ncov = 0) {

  # Argument checking.
  if (!(is.scalar(n) && n >= 2))
    stop("Input argument n should be an integer equal to 2 or more")
  if (!(is.scalar(p) && p >= 2))
    stop("Input argument p should be an integer equal to 2 or more")
  if (!(is.scalar(s) && s >= 0 && s <= p))
    stop("Input argument s should be an integer between 0 and p ",
         "(inclusive)")
  if (!(is.scalar(sigma) && sigma >= 0))
    stop("Input argument sigma should be a scalar greater than 0")
  if (!(is.scalar(pve) && pve >= 0 && pve <= 1))
    stop("Input argument pve should be a scalar between 0 ",
         "and 1 (inclusive)")
  if (!is.scalar(intercept))
    stop("Input argument intercept should be a scalar")
  if (!(is.scalar(ncov) && ncov >= 0))
    stop("Input argument ncov should be an greater than or equal to 0 ",
         "(inclusive)")
  if (!center_X && standardize_X)
    stop("Input argument center_X should be TRUE when standardize_X is TRUE")
  if (pve == 1 && sigma != 0)
    stop("If pve = 1, sigma should be zero")
  if (s == 0 && pve > 0)
    stop("If s = 0, then pve should also be zero")
  
  # Take floor of all integer arguments (in case a non-integer value
  # was given).
  n    <- floor(n)
  p    <- floor(p)
  s    <- floor(s)
  ncov <- floor(ncov)

  # Check that s <= p.
  if (s > p)
    stop("Input s should be no greater than p")

  # Simulate data matrix.
  X <- matrix(rnorm(n*p), n, p)
  X <- scale(X,center = center_X, scale = standardize_X)
  
  # Generate the regression coefficients.
  b.idx <- sample(p, s)
  b <- rep(0, p)
  if (s > 0) {
    b.values <- rnorm(s)
    b[b.idx] <- b.values

    # Adjust the effects so that we control for the proportion of
    # variance explained (pve). That is, we adjust b so that r =
    # a/(a+1), where we define a = b'*cov(X)*b. Here, sb is the
    # variance of the (nonzero) effects.
    if (pve < 1) {
      sb <- sqrt(pve/(1-pve)/var(drop(X %*% b)))
      b  <- sb * sigma * b
    }
  }

  # Generate the covariate data (Z), and the linear effects of the
  # covariates (u).
  if (ncov > 0) {
    Z <- matrix(rnorm(n*ncov),n,ncov)
    u <- rnorm(ncov)
  } else {
    Z <- as.numeric(NA)
    u <- as.numeric(NA)
  }

  # Generate the quantitative trait measurements.
  y <- intercept + X %*% b + sigma*rnorm(n)
  if (ncov > 0)
    y <- y + Z %*% u
  y <- drop(y)

  # Add row and column names to the outputs.
  rnames      <- paste0("s", 1:n)
  cnames_X    <- paste0("v", 1:p)
  names(y)    <- rnames
  rownames(X) <- rnames
  colnames(X) <- cnames_X
  names(b)    <- cnames_X
  if (ncov > 0) {
    cnames_Z    <- paste0("c", 1:ncov)
    names(u)    <- cnames_Z
    colnames(Z) <- cnames_Z
    rownames(Z) <- rnames
  }

  # Create the final output.
  return(list(y = y, X = X, Z = Z, u = u, b = b,
              intercept = intercept, sigma = sigma))
}
