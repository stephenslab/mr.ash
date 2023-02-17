#' @rdname fit_mr_ash
#'
#' @param X The data matrix, a numeric matrix of dimension n x p; each
#'   column is a single predictor, and each row is an observation
#'   vector. Here, n is the number of samples and p is the number of
#'   predictors.
#'
#' @param y The observed outcomes, a numeric vector of length n.
#'
#' @param intercept Should intercept be estimated (\code{intercept =
#'   TRUE}) or set to zero (\code{intercept = FALSE})? Please note that
#'   \code{intercept = FALSE} is generally not recommended and could
#'   produce unexpected results.
#'
#' @param standardize When \code{standardize = TRUE},
#'   \dQuote{standardize} the variables. (Internally, the columns of X
#'   are divided by the the standard deviations so that each column has
#'   a standard deviation of one.) Note that the coefficients are always
#'   returned on the original scale.
#'
#' @param b Optional input argument specifying the initial estimate of
#'   the regression coefficients. It should be numeric vector of length
#'   p.
#'
#' @param prior.sd Optional input argument specifying the standard
#'   deviations of the mixture components in the mixture-of-normals
#'   prior.
#'
#' @param prior.weights Optional initial estimate of the prior mixture
#'   proportions. This vector will automatically be normalized so that
#'   the entries represent proportions (that is, they add up to 1).
#'
#' @param resid.sd Initial estimate of the residual standard deviation.
#'
#' @param init.method Method used to initialize the estimates of the
#'   regression coefficients.  When \code{init.method = "glmnet"}, the
#'   estimates are initialized using \code{\link[glmnet]{cv.glmnet}}.
#'   (Note that cv.glmnet may give slightly different results depending
#'   on the state of the random number generator, so use \code{set.seed}
#'   to guarantee a predictable initialization.) When \code{init.ethod =
#'   "null"}, the estimates are initialized to zero.
#'
#' @param s The value of the glmnet penalty parameter at which the
#'   coeffients are extracted (relevant for \code{init.method =
#'   "glmnet"} only).
#'
#' @param \dots Additional arguments passed to
#'   \code{\link[glmnet]{cv.glmnet}} (relevant for \code{init.method =
#'   "glmnet"} only).
#'
#' @examples
#' set.seed(1)
#' dat <- simulate_regression_data(n = 400, p = 100, s = 20)
#' X <- dat$X
#' y <- dat$y
#'
#' # Initialize the coefficients using glmnet.
#' fit0_glmnet <- init_mr_ash(X, y, init.method = "glmnet")
#'
#' # Initialize the coefficients to be all zero.
#' fit0_null <- init_mr_ash(X, y, init.method = "null")
#'
#' # Randomly initialize the coefficients.
#' fit0_rand <- init_mr_ash(X, y, b = rnorm(100))
#'
#' # Specify a custom mixture prior.
#' fit0_custom <- init_mr_ash(X, y, prior.sd = (2^((0:19)/20) - 1),
#'                            prior.weights = rep(1,20))
#'
#' @importFrom stats sd
#' @importFrom matrixStats colSds
#'
#' @export
#'
init_mr_ash <- function (
  X, y,
  intercept = TRUE, standardize = FALSE,
  b, prior.sd, prior.weights, resid.sd,
  init.method = c("glmnet", "null"),
  s = "lambda.1se", ...) {

  # Check and process input argument X.
  if (!is.numeric.matrix(X))
    stop("Input argument X should be a numeric matrix, and all entries of ",
         "X should be finite and non-missing")
  if (is.integer(X))
    storage.mode(X) <- "double"
  n <- nrow(X)
  p <- ncol(X)
  if (n < 2 | p < 2)
    stop("Input matrix X should have at least 2 rows and 2 columns")
  if (any(colSds(X) <= 0))
    stop("Input matrix X cannot have any zero-variance columns")

  # Check and process input argument y.
  if (!is.numeric.vector(y))
    stop("Input argument y should be a vector and all entries should be ",
         "finite and non-missing")
  y <- as.vector(y, mode = "double")
  if (length(y) != n)
    stop("The length of y should be the same as nrow(X)")

  # Process input argument init.method.
  init.method <- match.arg(init.method)

  # Check optional inputs "intercept" and "standardize".
  if (!is.trueorfalse(intercept))
    stop("Input argument \"intercept\" should be TRUE or FALSE")
  if (!is.trueorfalse(standardize))
    stop("Input argument \"intercept\" should be TRUE or FALSE")
  X  <- scale(X, center = intercept, scale = standardize)
  y  <- drop(scale(y, center = intercept, scale = FALSE))
  mx <- attr(X,"scaled:center")
  my <- attr(y,"scaled:center")
  sx <- attr(X,"scaled:scale")

  # Check and process optional input b. If not provided, initialize
  # the coefficients using the chosen init.method.
  if (!missing(b)) {
    if (!is.numeric.vector(b))
      stop("Input argument b should be a vector and all entries should be ",
           "finite and non-missing")
    b <- as.vector(b, mode = "double")
    if (!length(b) == p)
      stop("Input argument b should have one entry for each column of X")
  } else {
    if (init.method == "null")
      b <- rep(0, p)
    else if (init.method == "glmnet")
      b <- init_coef_glmnet(X, y, s, ...)
  }

  # Check and process optioinal input resid.sd.
  if (!missing(resid.sd)) {
    if (!is.scalar(resid.sd) && resid.sd < 0)
      stop("Input argument resid.sd should be a number greater than zero")
  } else {

    # Initialize the residual s.d. to the MLE (assuming the
    # coefficients, b, are known).
    resid.sd <- sqrt((n-1)/n) * sd(y - X %*% b)
  }
  resid.sd <- as.vector(resid.sd, mode = "double")

  # Check and process optional input prior.sd.
  if (!missing(prior.weights) & missing(prior.sd))
    stop("If prior.weights is provided then prior.sd should also be provided")
  if (!missing(prior.sd)) {
    if (!is.numeric.vector(prior.sd))
      stop("Input argument prior.sd should be a vector and all entries ",
           "should be finite and non-missing")
    if (any(prior.sd < 0))
      stop("All entries of prior.sd should be non-negative")
    if (prior.sd[1] != 0)
      stop("The first entry of prior.sd should be zero")
    if (!all(diff(prior.sd) > 0))
      stop("The entries of prior.sd should be increasing")
  } else {

    # Set the standard deviations of the mixture components in an
    # automated way based on the data.
    prior.sd <- init_prior_sd(X, y)
  }
  prior.sd <- as.vector(prior.sd, mode = "double")

  # Check and process optional input prior.weights.
  if (!missing(prior.weights)) {
    if (!is.numeric.vector(prior.weights))
      stop("Input argument prior.weights should be a vector and all entries ",
           "should be finite and non-missing")
    if(length(prior.weights) != length(prior.sd))
      stop("prior.sd and prior.weights should have the same number of ",
           "elements")
    if (any(prior.weights < 0))
      stop("All entries of prior.sd should be non-negative")
    if (any(prior.weights == 0))
      warning("Mixture components with weights initialized to zero will ",
              "never be used")
  } else
    prior.weights <- init_prior_weights(X, b, resid.sd^2, prior.sd^2)
  prior.weights <- prior.weights / sum(prior.weights)
  prior.weights <- as.vector(prior.weights, mode = "double")
  k <- length(prior.weights)

  # Return the coefficients on the original scale, and compute the
  # initial estimate of the intercept (also on the original scale).
  if (standardize)
    b <- b/sx
  if (intercept)
    b0 <- my - sum(mx * b)
  else
    b0 <- 0

  # Prepare the final output.
  names(b)             <- colnames(X)
  names(prior.sd)      <- paste0("k",1:k)
  names(prior.weights) <- paste0("k",1:k)
  fit <- list(intercept   = intercept,
              standardize = standardize,
              b0          = b0,
              b           = b,
              resid.sd    = resid.sd,
              prior       = list(sd = prior.sd, weights = prior.weights),
              progress    = NULL)
  class(fit) <- c("mr.ash","list")
  return(fit)
}

# Initialize the posterior mean estimates of the regression
# coefficients using glmnet.
#
#' @importFrom glmnet cv.glmnet
#' @importFrom glmnet coef.glmnet
init_coef_glmnet <- function (X, y, s, ...) {
  fit <- cv.glmnet(X, y, intercept = FALSE, standardize = FALSE, ...)
  return(as.numeric(coef.glmnet(fit, s = s))[-1])
}

# Get a reasonable setting for the standard deviations of the mixture
# components in the mixture-of-normals prior based on the data (X, y).
# Input se is an estimate of the residual *variance*, and n is the
# number of standard deviations to return. This code is adapted from
# the autoselect.mixsd function in the ashr package.
init_prior_sd <- function (X, y, n = 20) {
  res <- simple_lr(X, y)
  smax <- 2*max(res$bhat)
  return(seq(0, smax, length.out = n))
}

# Initialize the mixture weights in the mixture-of-normals prior by
# findinging the mixture weights that "best fit" b if we treat b as
# the "true" posterior mean of the regression coefficients. Here, se
# and s0 are *variances* (not standard deviations); specifically, se
# is an estimate of the variance of the residual and s0 specifies the
# variances in the mixture-of-normals prior.
init_prior_weights <- function (X, b, se = 1, s0) {

  # Get the number of variables (n) and the number of mixture
  # components (k).
  n <- ncol(X)
  k <- length(s0)

  # Compute the variances of the least-squares estimates.
  shat <- simple_lr(X, se = se)$shat

  # Compute the Bayes factors separately for each mixture component.
  # The log-Bayes factors are stored in an n x k matrix. If s0 = 1,
  # then the Bayes factor is equal to 1 by definition (i.e., the
  # log-Bayes factor is 0).
  lbf <- matrix(0, n, k)
  for (i in 1:k)
    if (s0[i] > 0)
      lbf[,i] <- bayes_lr_ridge(b, shat, s0[i])$lbf

  # Compute the "responsibilities"; that is, the probability that the
  # effect for the ith variable was drawn from the kth component.
  P <- exp(lbf - apply(lbf, 1, max))
  P <- P/rowSums(P)

  # Compute the M-step update for the mixture weights.
  return(colMeans(P))
}
