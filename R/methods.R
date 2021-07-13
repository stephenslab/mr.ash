#' @title Predict Outcomes or Extract Coefficients from Mr.ASH Fit
#'
#' @description This function predicts outcomes (y) given the observed
#'   variables (X) and a Mr.ASH model; alternatively, retrieve the
#'   estimates of the regression coefficients.
#'
#' @param object A mr_ash fit, usually the result of calling
#'   \code{\link{mr_ash}}.
#'
#' @param newx The input matrix, of dimension (n,p); each column is a
#'   single predictor; and each row is an observation vector. Here, n is
#'   the number of samples and p is the number of predictors. When
#'   \code{newx} is \code{NULL}, the fitted values for the training data
#'   are provided.
#' 
#' @param type The type of output. For \code{type = "response"},
#'   predicted or fitted outcomes are returned; for \code{type =
#'   "coefficients"}, the estimated coefficients are returned.
#' 
#' @param ... Additional arguments passed to the default S3 method.
#'
#' @return For \code{type = "response"}, predicted or fitted outcomes
#' are returned; for \code{type = "coefficients"}, the estimated
#' coefficients are returned.
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
#' ## fit mr.ash model
#' fit.mr.ash  = mr_ash(X, y)
#' 
#' ## predict
#' Xnew        = matrix(rnorm(n*p),n,p)
#' ypred       = predict(fit.mr.ash, Xnew)
#' 
#' @importFrom stats predict
#' 
#' @export predict.mr.ash
#' 
#' @export
#' 
predict.mr.ash <- function (object, newx = NULL,
                            type=c("response","coefficients"),...) {
  
  type <- match.arg(type)
  if (type == "coefficients"){
    if(!missing(newx))
      stop("Do not supply newx when predicting coefficients")
    return(coef(object))
  }
  else if(missing(newx))
    return(object$fitted)
  else {
    if (!all(object$data$Z == 1))
      stop("predict.mr.ash is not implemented for covariates Z other than ",
           "intercept")
    return(drop(object$intercept + newx %*% coef(object)[-1]))
  }
}

#' @title Extract Regression Coefficients from Mr.ASH Fit
#'
#' @description Retrieve posterior mean estimates of the regression
#'   coefficients in a Mr.ASH model.
#' 
#' @param object A Mr.ASH fit, usually the result of calling
#'   \code{mr.ash}.
#'
#' @param ... Additional arguments passed to the default S3 method.
#' 
#' @return A p+1 vector. The first element gives the estimated
#'   intercept, and the remaining p elements are the estimated
#'   regression coefficients.
#'   
#' ## generate synthetic data
#' set.seed(1)
#' n           = 200
#' p           = 300
#' X           = matrix(rnorm(n*p),n,p)
#' beta        = double(p)
#' beta[1:10]  = 1:10
#' y           = X %*% beta + rnorm(n)
#' 
#' ## fit mr.ash model
#' fit.mr.ash  = mr_ash(X, y)
#' 
#' ## coefficient
#' coef.mr.ash = coef(fit.mr.ash)
#' intercept   = coef.mr.ash[1]
#' beta        = coef.mr.ash[-1]
#' 
#' @importFrom stats coef
#' 
#' @export coef.mr.ash
#' 
#' @export
#' 
coef.mr.ash <- function (object, ...)
  c(object$intercept,object$beta)
