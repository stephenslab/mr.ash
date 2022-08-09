#' @rdname fit_mr_ash
#'
#' @param sa2 The vector of prior mixture component variances. The
#'   variances should be in increasing order, starting at zero; that is,
#'   \code{sort(sa2)} should be the same as \code{sa2}. When \code{sa2}
#'   is \code{NULL}, the default setting is used, \code{sa2[k] =
#'   (2^(0.05*(k-1)) - 1)^2}, for \code{k = 1:20}. For this default
#'   setting, \code{sa2[1] = 0}, and \code{sa2[20]} is roughly 1.
#'
#' @param beta.init The initial estimate of the (approximate)
#'   posterior mean regression coefficients. This should be \code{NULL},
#'   or a vector of length p. When \code{beta.init} is \code{NULL}, the
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
#'   initial estimates of the regression coefficients, \code{beta.init},
#'   after removing linear effects of the intercept and any covariances.
#'
#' @return A \code{mr.ash} object
#'
#' @export
#'
init_mr_ash <- function (
  sa2, beta.init, pi, sigma2
) {

  # additional check on pi if provided
  if (
    (!missing(pi) && !missing(sa2) && length(pi) != length(sa2)) ||
    (!missing(pi) && is.null(sa2) && length(pi) != 20)
  ) {

    stop("pi and sa2 must be of the same length.",
         " Either provide an sa2 with the same length as pi",
         " or provide a pi of length 20 to match the default length of sa2")

  }

  if (!missing(sa2)) {
    if (any(sa2 < 0))
      stop ("all the mixture component variances must be non-negative.")
    if (sa2[1] != 0)
      stop ("the first mixture component variance sa2[1] must be 0.")
    if (!all(sort(sa2) == sa2))
      stop ("sa2 must be sorted")
  } else {
    sa2 <- NULL
  }

  if (!missing(beta.init)) {
    if (!is.numeric(beta.init))
      stop("beta.init must be a numeric vector")
  }
  else {
    beta.init <- NULL
  }

  if (!missing(pi)) {

    if (!is.numeric(pi))
      stop("pi must be a numeric vector")
    if (!missing(sa2) && length(sa2) != length(pi))
      stop("pi and sa2 must be of the same length")

  } else {
    pi <- NULL
  }


  if (!missing(sigma2)) {

    if (!is.numeric(sigma2))
      stop("sigma2 must be numeric")
    if (length(sigma2) != 1)
      stop("sigma2 must be a single number")
    if (sigma2 <= 0)
      stop("sigma2 must be greater than 0")

  } else {
    sigma2 <- NULL
  }

 fit0 <- list(
   sa2 = sa2, beta.init = beta.init, pi = pi, sigma2 = sigma2
 )

 class(fit0) <- c("mr.ash","list")
 return(fit0)

}
