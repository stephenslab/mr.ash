#' @title Simulate Regression Data
#'
#' @description Add descrption here.
#'
#' @importFrom stats rnorm
#' 
#' @export
#' 
simulate_regression_data <- function (n, p) {
  X <- matrix(rnorm(n*p),n,p)
  return(list(X = X))
}
