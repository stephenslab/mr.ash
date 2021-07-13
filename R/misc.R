# Remove covariate effects Regresses Z out from X and y; that is, X
# and y are projected into the space orthogonal to Z.
#' 
#' @importFrom Matrix forceSymmetric
#'
remove_covariate <- function (X, y, Z, standardize = FALSE, intercept = TRUE) {
  
  # check if Z is null and intercept = FALSE
  if (is.null(Z) & (intercept == FALSE)) {
    return(list(X = X, y = y, Z = Z,
                ZtZiZX = rep(0,dim(X)[2]), ZtZiZy = 0))
  }
  
  # redefine y
  y = c(as.double(y))
  n = length(y)
  
  # add intercept if intercept = TRUE
  if (intercept) {
    if (is.null(Z))
      Z <- matrix(1,n,1)
    else
      Z <- cbind(1,Z)
  }
  
  if (ncol(Z) == 1) {
    ZtZ         = forceSymmetric(crossprod(Z))       # (Z^T Z) symmetric
    ZtZiZy      = as.vector(solve(ZtZ,c(y %*% Z)))   # (Z^T Z)^{-1} Z^T y
    ZtZiZX      = as.matrix(solve(ZtZ,t(Z) %*% X))   # (Z^T Z)^{-1} Z^T X
    X           = scale(X, center = intercept, scale = standardize)
    alpha       = mean(y)
    y           = y - alpha
    
  } else {
    ZtZ         = forceSymmetric(crossprod(Z))       # (Z^T Z) symmetric
    ZtZiZy      = as.vector(solve(ZtZ,c(y %*% Z)))   # (Z^T Z)^{-1} Z^T y
    ZtZiZX      = as.matrix(solve(ZtZ,t(Z) %*% X))   # (Z^T Z)^{-1} Z^T X
    y     = y - c(Z %*% ZtZiZy)
    X     = X - Z %*% ZtZiZX
  }
  
  return(list(X = X, y = y, Z = Z,
              ZtZiZX = ZtZiZX, ZtZiZy = ZtZiZy))
}

var.n <- function(x) {
  a <- x - mean(x)
  return (sum(a^2) / length(a))
}
