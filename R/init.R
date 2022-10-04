#' @rdname fit_mr_ash
#'
#' @param X The input matrix, of dimension (n,p); each column is a
#'   single predictor; and each row is an observation vector. Here, n is
#'   the number of samples and p is the number of predictors. The matrix
#'   cannot be sparse.
#'
#' @param y The observed continuously-valued responses, a vector of
#'   length n.
#'
#' @param sa2 The vector of prior mixture component variances. The
#'   variances should be in increasing order, starting at zero; that is,
#'   \code{sort(sa2)} should be the same as \code{sa2}. When \code{sa2}
#'   is \code{NULL}, the default setting is used, \code{sa2[k] =
#'   (2^(0.05*(k-1)) - 1)^2}, for \code{k = 1:20}. For this default
#'   setting, \code{sa2[1] = 0}, and \code{sa2[20]} is roughly 1.
#'
#' @param b The initial estimate of the (approximate)
#'   posterior mean regression coefficients. This should be \code{NULL},
#'   or a vector of length p. When \code{b} is \code{NULL}, the
#'   posterior mean coefficients are all initially set to zero.
#'
#' @param pi The initial estimate of the mixture proportions
#'   \eqn{\pi_1, \ldots, \pi_K}. If \code{pi} is \code{NULL}, the
#'   mixture weights are initialized to \code{rep(1/K,K)}, where
#'   \code{K = length(sa2)}.
#'
#' @param sigma2 The initial estimate of the residual variance,
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
#' @return A \code{mr.ash} object
#'
#' @examples
#'
#' # initialization examples
#' data <- simulate_regression_data(n = 100, p = 100, s = 50)
#'
#' # initialize with default values for all parameters
#' fit0_def <- init_mr_ash(data$X, data$y)
#'
#' # initialize b with all 0s
#' fit0_null <- init_mr_ash(data$X, data$y, init.method = "null")
#'
#' # specify custom mixture distribution for b
#' fit0_custom_mixt <- init_mr_ash(
#'   data$X, data$y, sa2 = (2^((0:19) / 20) - 1)^2, pi = rep(1/20, 20)
#' )
#'
#' @export
#'
init_mr_ash <- function (
  X, y, sa2, b, pi, sigma2, init.method = c("glmnet", "null")
) {

  # get sizes
  n <- nrow(X)
  p <- ncol(X)

  init.method <- match.arg(init.method)

  if (!missing(b)) {
    if (!is.numeric(b))
      stop("b must be a numeric vector")
    if (!length(b) == p)
      stop("b must have the same length as number of columns in X")
  }
  else {
    if (init.method == "null") {
      b <- rep(0,p)
    } else if (init.method == "glmnet") {
      lasso_cv_fit <- glmnet::cv.glmnet(X, y)

      # extract coefficients, excluding intercept
      b <- as.vector(coef(lasso_cv_fit, s = "lambda.min"))[-1]
    }
  }

  # check that the dimensions of X and y match
  if (length(y) != n) {
    stop("The length of y must match the number of rows of X")
  }

  if (!is.numeric(X) || !is.numeric(y)) {
    stop("X and y must be numeric")
  }

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

    pi <- pi / sum(pi)

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

 fit0 <- list(
   sa2 = sa2, b = b, pi = pi, sigma2 = sigma2, progress = data.frame()
 )

 class(fit0) <- c("mr.ash","list")
 return(fit0)

}
