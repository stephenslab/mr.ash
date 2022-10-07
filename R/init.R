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
#' @param prior.weights Optional initial estimate of the prior mixture
#'   proportions.
#'
#' @param resid.sd Initial estimate of the residual standard deviation.
#'
#' @param init.method Method used to initialize the estimates of the
#'   regression coefficients.  When \code{init.method = "glmnet"}, the
#'   estimates are initialized using
#'   \code{\link[glmnet]{cv.glmnet}}. When \code{init.ethod = "null"},
#'   the estimates are initialized to zero.
#'
#' @param s The value of the glmnet penalty parameter at which the
#'   coeffients are extracted (relevant for \code{init.method =
#'   "glmnet"} only).
#' 
#' @param \dots Additional arguments passed to
#'   \code{\link[glmnet]{cv.glmnet}} (relevant for \code{init.method =
#'   "glmnet"} only).
#' 
#' @return A \code{mr.ash} object.
#'
#' @seealso \code{\link{fit_mr_ash}}
#' 
#' @examples
#' dat <- simulate_regression_data(n = 400, p = 100, s = 20)
#' X <- dat$X
#' y <- dat$y
#'
#' # Initialize the coefficients using glmnet.
#' fit0_glmnet <- init_mr_ash(X, y)
#'
#' # Initialize the coefficients to zero.
#' fit0_null <- init_mr_ash(X, y, init.method = "null")
#'
#' # Randomly initialize the coefficients.
#' fit0_rand <- init_mr_ash(X, y, b = rnorm(100))
#' 
#' # specify custom mixture distribution for b
#' fit0_custom_mixt <- init_mr_ash(
#'   dat$X, dat$y, sa2 = (2^((0:19) / 20) - 1)^2, pi = rep(1/20, 20)
#' )
#'
#' @importFrom stats sd
#' 
#' @export
#'
init_mr_ash <- function (
  X, y, b, prior.sd, prior.weights, resid.sd,
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

  # Check and process input argument y.
  if (!is.numeric.vector(y))
    stop("Input argument y should be a vector and all entries should be ",
         "finite and non-missing")
  y <- as.vector(y,mode = "double")
  if (length(y) != n)
    stop("The length of y should be the same as nrow(X)")

  # Process input argument init.method.
  init.method <- match.arg(init.method)

  # Check and process optional input b. If not provided, initialize
  # the coefficients using the chosen init.method.
  if (!missing(b)) {
    if (!is.numeric.vector(b))
      stop("Input argument b should be a vector and all entries should be ",
           "finite and non-missing")
    b <- as.vector(b,mode = "double")
    if (!length(b) == p)
      stop("Input argument b should have one entry for each column of X")
  } else {
    if (init.method == "null")
      b <- rep(0,p)
    else if (init.method == "glmnet")
      b <- init_coef_glmnet(X,y,s,...)
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
  resid.sd <- as.vector(resid.sd,mode = "double")
  
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
    prior.sd <- init_resid_sd(X,y,resid.sd^2)
  }
  prior.sd <- as.vector(prior.sd,mode = "double")
    
  # Prepare the final output.
  names(b) <- colnames(X)
  fit <- list(b         = b,
              resid.sd  = resid.sd,
              prior     = list(sd = prior.sd),
              progress  = NULL)
  class(fit) <- c("mr.ash","list")
  return(fit)

  # Check and process optional input prior.weights.
  # TO DO

  if (TRUE) {
    if(length(pi) != length(sa2))
      stop("pi and sa2 must be of the same length")
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

    w <- colSums(X^2)
    S   <- outer(1/w, sa2, '+') * sigma2
    Phi <- -b^2/S/2 - log(S)/2
    Phi <- exp(Phi - apply(Phi,1,max))
    Phi <- Phi / rowSums(Phi)
    pi  <- colMeans(Phi)
  }

  # Prepare the final output.
}

# Initialize the posterior mean estimates of the regression
# coefficients using glmnet.
#
#' @importFrom glmnet cv.glmnet
#' @importFrom glmnet coef.glmnet
init_coef_glmnet <- function (X, y, s, ...) {
  fit <- cv.glmnet(X, y, ...)
  return(drop(coef.glmnet(fit, s = s))[-1])
}

# Get a reasonable setting for the standard deviations of the mixture
# components in the mixture-of-normals prior based on the data (X, y).
# Input se is an estimate of the residual variance, and n is the
# number of standard deviations to return. This code is based on the
# autoselect.mixsd function from the ashr package.
init_resid_sd <- function (X, y, se = 1, n = 20) {
  res <- bayes_lr_ridge(X, y, se)
  smax <- with(res,2*sqrt(max(bhat^2 - shat)))
  return(seq(0, smax, length.out = n))
}
