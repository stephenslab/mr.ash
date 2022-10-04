#' @rdname fit_mr_ash
#'
#' @param X The data matrix, a numeric matrix of dimension n x p; each
#'   column is a single predictor, and each row is an observation
#'   vector. Here, n is the number of samples and p is the number of
#'   predictors.
#'
#' @param y The observed outcomes, a numeric vector of length n.
#'
#' @param b Optional input argument specifying the initial estimate of
#'   the regression coefficients. It should be numeric vector of length
#'   p.
#'
#' @param prior.sd Optional input argument specifying the standard
#'   deviations of the mixture components in the mixture-of-normals
#'   prior.
#'
#' @param prior.weights Optio initial estimate of the mixture proportions
#'   \eqn{\pi_1, \ldots, \pi_K}. If \code{pi} is \code{NULL}, the
#'   mixture weights are initialized to \code{rep(1/K,K)}, where
#'   \code{K = length(sa2)}.
#'
#' @param resid.sd The initial estimate of the residual variance,
#'   \eqn{\sigma^2}. If \code{sigma2 = NULL}, the residual variance is
#'   initialized to the empirical variance of the residuals based on the
#'   initial estimates of the regression coefficients, \code{b},
#'   after removing linear effects of the intercept and any covariances.
#'
#' @param init.method method used to initialize \code{b} if not
#'   provided.  When set to \code{"glmnet"}, an L1 penalized regression
#'   of \code{y} on \code{X} is run (using cross-validation to select
#'   the magnitude of regularization) and the resulting regression
#'   coefficients are used as the initial values of \code{b}. It set to
#'   \code{"null"}, all regression coefficients are initialized to
#'   \code{0}. See \code{\link[glment]{cv.glmnet}} for more details.
#'
#' @return A \code{mr.ash} object.
#'
#' @seealso \code{\link{fit_mr_ash}}
#' 
#' @examples
#' dat <- simulate_regression_data(n = 100, p = 100, s = 50)
#'
#' # initialize with default values for all parameters
#' fit0_def <- init_mr_ash(dat$X, dat$y)
#'
#' # initialize b with all 0s
#' fit0_null <- init_mr_ash(dat$X, dat$y, init.method = "null")
#'
#' # specify custom mixture distribution for b
#' fit0_custom_mixt <- init_mr_ash(
#'   dat$X, dat$y, sa2 = (2^((0:19) / 20) - 1)^2, pi = rep(1/20, 20)
#' )
#'
#' @export
#'
init_mr_ash <- function (
  X, y, b, prior.sd, prior.weights, resid.sd,
  init.method = c("glmnet", "null")) {

  # Check and process input argument X.
  if (!(is.matrix(X) & is.numeric(X)))
    stop("Input argument X should be a numeric matrix")
  if (any(is.infinite(X)) | anyNA(X))
    stop("All entries of X should be finite and non-missing")
  if (is.integer(X))
    storage.mode(X) <- "double"
  n <- nrow(X)
  p <- ncol(X)
  if (n < 2 | p < 2)
    stop("Input matrix X should have at least 2 rows and 2 columns")

  # Check and process input argument y.
  if (!is.numeric(y))
    stop("X and y must be numeric")
  if (any(is.infinite(y)) | anyNA(y))
    stop("All entries of y should be finite and non-missing")
  y <- as.vector(y,mode = "double")
  if (length(y) != n)
    stop("The length of y should be the same as nrow(X)")

  # Process input argument init.method.
  init.method <- match.arg(init.method)

  # Check and process optional input b. If not provided, initialize
  # the coefficients using the chosen init.method.
  if (!missing(b)) {
    if (!is.numeric(b))
      stop("b must be a numeric vector")
    if (!length(b) == p)
      stop("b must have the same length as number of columns in X")
  } else {
    if (init.method == "null")
      b <- rep(0,p)
    else if (init.method == "glmnet") {
      lasso_cv_fit <- glmnet::cv.glmnet(X, y)

      # extract coefficients, excluding intercept
      b <- as.vector(coef(lasso_cv_fit, s = "lambda.min"))[-1]
    }
  }

  browser()
  
  # Check and process optional input prior.sd. 
  # TO DO
  
  # Check and process optional input prior.weights.
  # TO DO

  # Check and process optioinal input resid.sd
  # TO DO
  
  # check that the dimensions of X and y match

  if (!missing(sigma2)) {
    if (!is.numeric(sigma2))
      stop("sigma2 must be numeric")
    if (length(sigma2) != 1)
      stop("sigma2 must be a single number")
    if (sigma2 <= 0)
      stop("sigma2 must be greater than 0")
  } else {
    # initialize r
    r <- drop(y - X %*% b)
    sigma2 <- c(var.n(r))
  }

  if ((missing(pi) && !missing(sa2)) || (!missing(pi) && missing(sa2))) {
    stop("Either both pi and sa2 should be specified",
         "or neither should be specified")
  } else if (!missing(pi) && !missing(sa2)) {

    if(length(pi) != length(sa2))
      stop("pi and sa2 must be of the same length")

    if (any(sa2 < 0))
      stop("all the mixture component variances must be non-negative.")
    if (sa2[1] != 0)
      stop("the first mixture component variance sa2[1] must be 0.")
    if (!all(sort(sa2) == sa2))
      stop("sa2 must be sorted")

    if (!is.numeric(pi))
      stop("pi must be a numeric vector")
    if (!missing(sa2) && length(sa2) != length(pi))
      stop("pi and sa2 must be of the same length")
    if (any(pi < 0))
      stop("all elements of pi must be greater than or equal to 0")
    if (any(pi == 0))
      warning("some element(s) of pi are 0.",
              " This mixture component will not change with model training")

    pi <- pi/sum(pi)

  } else {

    d <- diag(crossprod(X))
    bhat <-	drop(t(X) %*% y)/d
    sehat <- sqrt(sigma2/d)
    sigmaamax <- 2*sqrt(max(bhat^2-sehat^2))
    sa2 <- seq(from = 0, to = sigmaamax, length.out = 20)

    w <- colSums(X^2)
    S   <- outer(1/w, sa2, '+') * sigma2
    Phi <- -b^2/S/2 - log(S)/2
    Phi <- exp(Phi - apply(Phi,1,max))
    Phi <- Phi / rowSums(Phi)
    pi  <- colMeans(Phi)

  }

  # TO DO: Add names to all the outputs.
  
  # Prepare the final output.
  fit0 <- list(
    sa2 = sa2, b = b, pi = pi, sigma2 = sigma2, progress = data.frame()
  )
  class(fit0) <- c("mr.ash","list")
  return(fit0)
}
