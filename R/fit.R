#' @rdname fit_mr_ash
#'
#' @title Multiple Regression with Adaptive Shrinkage
#'
#' @description Model fitting algorithms for Multiple Regression with
#'   Adaptive Shrinkage ("Mr.ASH"). Mr.ASH is a variational empirical
#'   Bayes (VEB) method for multiple linear regression. The fitting
#'  algorithms (locally) maximize the approximate marginal likelihood
#'   (the "evidence lower bound", or ELBO) via coordinate-wise updates.
#'
#' @details Mr.ASH is a statistical inference method for the following
#' multiple linear regression model: \deqn{y | X, \beta, \sigma^2 ~
#' N(X \beta, \sigma I_n),} in which the regression coefficients
#' \eqn{\beta} admit a mixture-of-normals prior, \deqn{\beta | \pi,
#' \sigma ~ g = \sum_{k=1}^K N(0, \sigma^2 \sigma_k^2).} Each mixture
#' component in the prior, \eqn{g}, is a normal density centered at
#' zero, with variance \eqn{\sigma^2 \sigma_k^2}.
#'
#' The fitting algorithm, if run for a large enough number of
#' iterations, will find an approximate posterior for the regression
#' coefficients, denoted by \eqn{q(\beta)}, residual variance
#' parameter \eqn{\sigma^2}, and prior mixture weights \eqn{\pi_1,
#' \ldots, \pi_K} maximizing the evidence lower bound, \deqn{F(q, \pi,
#' \sigma^2) = E_q log p(y | X, \beta, \sigma^2) - \sum_{j=1}^p
#' D_{KL}(q_j || g),} where \eqn{D_{KL}(q_j || g)} denotes the
#' Kullback-Leibler (KL) divergence, a measure of the "distance"
#' between (approximate) posterior \eqn{q_j(\beta_j)} and prior
#' \eqn{g(\beta_j)}. The fitting algorithm iteratively updates the
#' approximate posteriors \eqn{q_1, \ldots, q_p}, separately for each
#' \eqn{j = 1, \ldots, p} (in an order determined by
#' \code{update.order}), then separately updates the mixture weights
#' \eqn{\pi} and residual variance \eqn{\sigma^2}. This
#' coordinate-wise update scheme iterates until the convergence
#' criterion is met, or until the algorithm hits an upper bound on
#' the number of iterations (specified by \code{max.iter}).
#'
#' See \sQuote{References} for more details about the model and
#' algorithm.
#'
#' The \code{control} argument is a list in which any of the following
#' named components will override the default algorithm settings (as
#' defined by \code{mr_ash_control_default}):
#'
#' \describe{
#'
#' \item{\code{min.iter}}{The minimum number of outer loop iterations.}
#'
#' \item{\code{max.iter}}{The maximum number of outer loop iterations.}
#'
#' \item{\code{convtol}}{When \code{update.prior.weights = TRUE}, the
#'   outer-loop updates terminate when the largest change in the mixture
#'   weights is less than \code{convtol*K}; when \code{update.prior.weights =
#'   FALSE}, the outer-loop updates stop when the largest change in the
#'   estimates of the posterior mean coefficients is less than
#'   \code{convtol*K}.}
#'
#' \item{\code{update.prior.weights}}{If \code{update.prior.weights = TRUE}, the mixture
#'   proportions in the mixture-of-normals prior are estimated from the
#'   data.}
#'
#' \item{\code{update.resid.sd}}{If \code{update.resid.sd = TRUE}, the
#'   residual variance parameter \eqn{\sigma^2} is estimated from the
#'   data.}
#'
#' \item{\code{update.order}}{The order in which the coordinate
#'   ascent updates for estimating the posterior mean coefficients are
#'   performed. \code{update.order} can be \code{NULL}, \code{"random"},
#'   or any permutation of \eqn{(1,\ldots,p)}, where \code{p} is the
#'   number of columns in the input matrix \code{X}. When
#'   \code{update.order} is \code{NULL}, the coordinate ascent updates
#'   are performed in order in which they appear in \code{X}; this is
#'   equivalent to setting \code{update.order = 1:p}. When
#'   \code{update.order = "random"}, the coordinate ascent updates are
#'   performed in a randomly generated order, and this random ordering
#'   is different at each outer-loop iteration.}
#'
#' \item{\code{epstol}}{A small, positive number added to the
#'   likelihoods to avoid logarithms of zero.}}
#'
#' @param X The input matrix, of dimension (n,p); each column is a
#'   single predictor; and each row is an observation vector. Here, n is
#'   the number of samples and p is the number of predictors. The matrix
#'   cannot be sparse.
#'
#' @param y The observed continuously-valued responses, a vector of
#'   length n.
#'
#' @param fit0 Initialized \code{mr.ash} object resulting from a call to
#' \code{init_mr_ash}.
#'
#' @param standardize The logical flag for standardization of the
#'   columns of X variable, prior to the model fitting. The coefficients
#'   are always returned on the original scale.
#'
#' @param intercept When \code{intercept = TRUE}, an intercept is
#'   included in the regression model.
#'
#' @param control A list of parameters controlling the behaviour of
#'   the optimization algorithm. See \sQuote{Details}.
#'
#' @param verbose When \code{verbose = "detailed"}, detailed
#'   information about the algorithm's progress is printed to the
#'   console at each iteration; when \code{verbose = "progressbar"}, a
#'   plus (\dQuote{+}) is printed for each outer-loop iteration; and
#'   when \code{verbose = "none"}, no progress information is printed.
#'
#' @return A list object with the following elements:
#'
#' \item{intercept}{The estimated intercept.}
#'
#' \item{b}{Posterior mean estimates of the regression coefficients.}
#'
#' \item{resid.sd}{The estimated residual variance.}
#'
#' \item{posterior.weights}{A vector of containing the estimated mixture
#'   proportions.}
#'
#' \item{progress}{A list containing estimated parameters and objectives
#' over each outer-loop iteration.}
#'
#' \item{update.order}{The ordering used for performing the
#'   coordinate-wise updates. For \code{update.order = "random"}, the
#'   orderings for outer-loop iterations are provided in a vector of
#'   length \code{p*max.iter}, where \code{p} is the number of predictors.}
#'
#' \item{elbo}{Evidence lower bound of final fit.}
#'
#' \item{phi}{A p x K matrix containing the posterior assignment
#'   probabilities, where p is the number of predictors, and K is the
#'   number of mixture components. (Each row of \code{phi} should sum to
#'   1.)}
#'
#' \item{m}{A p x K matrix containing the posterior means conditional
#'   on assignment to each mixture component.}
#'
#' \item{s2}{A p x K matrix containing the posterior variances
#'   conditional on assignment to each mixture component.}
#'
#' \item{lfsr}{A vector of length p containing the local false
#'   discovery rate for each variable.}
#'
#' \item{prior}{list of values used in prior}
#'
#' \item{fitted}{fitted values for each row of \code{X}.}
#'
#'
#' @seealso \code{\link{predict.mr.ash}}
#'
#' @references
#' Y. Kim (2020), Bayesian shrinkage methods for high dimensional
#' regression. Ph.D. thesis, University of Chicago.
#'
#' @useDynLib mr.ash
#'
#' @importFrom Rcpp evalCpp
#' @importFrom utils modifyList
#' @importFrom stats var
#'
#' @examples
#' # Simulate a data set.
#' set.seed(1)
#' n <- 200
#' p <- 300
#' data <- simulate_regression_data(n, p, s = 250)
#'
#' ### fit Mr.ASH
#' fit.mr.ash <- fit_mr_ash(data$X,data$y)
#'
#' ### prediction routine
#' Xnew        = matrix(rnorm(n*p),n,p)
#' ynew        = Xnew %*% data$b + rnorm(n)
#' ypred       = predict(fit.mr.ash, Xnew)
#'
#' ### test error
#' rmse        = norm(ynew - ypred, '2') / sqrt(n)
#'
#' ### coefficients
#' bhat     = predict(fit.mr.ash, type = "coefficients")
#' # this equals c(fit.mr.ash$intercept, fit.mr.ash$b)
#'
#' @export
#'
fit_mr_ash <- function (
  X, y, fit0 = init_mr_ash(X, y), standardize = FALSE, intercept = TRUE,
  control = list(), verbose = c("progress", "detailed", "none")) {

  n <- nrow(X)
  p <- ncol(X)

  if (!is.null(standardize)) {

    if (!is.logical(standardize))
      stop("standardize must be set to either TRUE or FALSE")

  }

  if (!is.null(intercept)) {

    if (!is.logical(intercept))
      stop("intercept must be set to either TRUE or FALSE")

  }

  if(!is.list(control))
    stop("control must be a list")

  if (!inherits(fit0, "mr.ash"))
    stop("fit0 must be of class mr.ash")

  # Check and process input argument "control".
  control <- modifyList(mr_ash_control_default(),control,keep.null = TRUE)

  # Check and process input argument "verbose".
  verbose <- match.arg(verbose)

  # write off original X for fitted values
  og_data_X <- X

  # remove covariates
  res    <- remove_covariate(X,y,NULL,standardize,intercept)
  X      <- res$X
  y      <- res$y
  ZtZiZX <- res$ZtZiZX
  ZtZiZy <- res$ZtZiZy

  # standard beta

  if (standardize)
    beta <- drop(fit0$b) * attr(X,"scaled:scale")
  else
    beta <- drop(fit0$b)

  # initialize r
  r <- drop(y - X %*% beta)

  # precompute x_j^T x_j
  w <- colSums(X^2)

  K <- length(fit0$prior$sd)

  # run algorithm
  if (is.null(control$update.order))
    o <- rep(seq(0,p-1),control$max.iter)
  else if (is.numeric(control$update.order))
    o <- rep(control$update.order - 1,control$max.iter)
  else if (control$update.order == "random")
    o <- random_order(p,control$max.iter)
  method_q <- "sigma_dep_q"
  if (verbose != "none") {
    cat("Fitting mr.ash model (mr.ash 0.1-88).\n")
    cat(sprintf("number of samples: %d\n",n))
    cat(sprintf("number of variables: %d\n",p))
    cat(sprintf("number of mixture components: %d\n",K))
  }
  if (verbose == "detailed")
    cat("iter                elbo ||b-b'||   sigma2 w>0\n")
  out <- mr_ash_rcpp(X, y, w, fit0$prior$sd, fit0$prior$weights, beta, as.vector(r),
                     fit0$resid.sd, o, control$max.iter, control$min.iter,
                     control$convtol, control$epstol, method_q,
                     control$update.prior.weights, control$update.resid.sd,
                     switch(verbose, none = 0, progress = 1, detailed = 2))

  # polish return object
  out$progress <- data.frame(iter   = 1:control$max.iter,
                             elbo   = -out$varobj,
                             dbeta  = out$dbeta,
                             sigma2 = out$sigma2byiter,
                             w1     = out$w1)
  out$progress <- out$progress[1:out$iter,]
  out <- out[c("b","resid.sd","posterior.weights","progress")]
  out$elbo <- tail(out$progress$elbo,n = 1)
  out$intercept <- c(ZtZiZy - ZtZiZX %*% out$b)
  out$update.order <- o

  # rescale beta as needed
  if (standardize)
    out$b <- out$b / attr(X,"scaled:scale")

  ## warn if necessary
  if (control$update.prior.weights && out$posterior.weights[K] > 1/K)
    warning(sprintf(paste("The mixture proportion associated with the",
                          "largest prior variance is greater than %0.2e;",
                          "this indicates that the model fit could be",
                          "improved by using a larger setting of the",
                          "prior variance. Consider increasing the range",
                          "of the variances \"sa2\"."),1/K))

  # Add dimension names and prepare the final fit object.
  res <- get_full_posterior(
    X,y,w,out$b,out$posterior.weights,out$resid.sd,fit0$prior$sd
  )

  out$m <- res$m
  out$s2 <- res$s2
  out$phi <- res$phi
  out$lfsr <- res$lfsr
  out$fitted <- as.vector(og_data_X %*% out$b + out$intercept)
  out$b <- drop(out$b)
  out$posterior.weights <- drop(out$posterior.weights)

  out$prior <- list()
  out$prior$sd <- fit0$prior$sd
  out$prior$weights <- fit0$prior$weights

  names(out$b) <- colnames(X)
  class(out) <- c("mr.ash","list")
  return(out)
}

#' @title Approximation Posterior Expectations from Mr.ASH Fit
#'
#' @description Recover the parameters specifying the variational
#'   approximation to the posterior distribution of the regression
#'   coefficients. To streamline the model fitting implementation, and
#'   to reduce memory requirements, \code{\link{fit_mr_ash}} does not store
#'   all the parameters needed to specify the approximate posterior.
#'
#' @param fit A Mr.ASH fit obtained, for example, by running
#'   \code{\link{fit_mr_ash}}.
#'
#' @return A list object with the following elements:
#'
#' \item{phi}{A p x K matrix containing the posterior assignment
#'   probabilities, where p is the number of predictors, and K is the
#'   number of mixture components. (Each row of \code{phi} should sum to
#'   1.)}
#'
#' \item{m}{A p x K matrix containing the posterior means conditional
#'   on assignment to each mixture component.}
#'
#' \item{s2}{A p x K matrix containing the posterior variances
#'   conditional on assignment to each mixture component.}
#'
#' \item{lfsr}{A vector of length p containing the local false
#'   discovery rate for each variable}
#'
#' @examples
#' ## generate synthetic data
#' set.seed(1)
#' n           = 200
#' p           = 300
#' X           = matrix(rnorm(n*p),n,p)
#' beta        = double(p)
#' beta[1:10]  = 1:10
#' y           = X %*% beta + rnorm(n)
#'
#' # Fit mr.ash model.
#' fit <- fit_mr_ash(X, y)
#'
#' @export
#'
get_full_posterior <- function (X, y, w, beta, pi, sigma2, sa2) {

  # compute residual
  r = y - X %*% beta

  # compute bw and s2
  bw = as.vector((t(X) %*% r) + w * beta)
  s2 = sigma2 / outer(w, 1/sa2, '+')

  # compute m, phi
  m   = bw * s2/ sigma2
  phi = -log(1 + outer(w,sa2))/2 + m * (bw/2/sigma2)
  phi = c(pi) * t(exp(phi - apply(phi,1,max)))
  phi = t(phi) / colSums(phi)

  # compute lfsr
  lfsr <- computelfsrmix(phi, m, s2)
  return (list(phi = phi, m = m, s2 = s2, lfsr = lfsr))
}

#' @rdname fit_mr_ash
#'
#' @export
#'
mr_ash_control_default <- function()
  list(min.iter      = 1,
       max.iter      = 1000,
       convtol       = 1e-8,
       update.prior.weights     = TRUE,
       update.resid.sd = TRUE,
       update.order  = NULL,
       epstol        = 1e-12)
