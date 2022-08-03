#' @title Simulate Regression Data
#'
#' @description Simulate non-sparse data matrix X (and a corresponding normal y?)
#'
#' @param n number of samples
#' 
#' @param p number of variables
#' 
#' @param s number of non-zero effects
#' 
#' @param se standard deviation of the residual
#' 
#' @param pve proportion of variance in Y explained by X
#' 
#' @param misc FILL IN OTHER PARAMS
#' 
#' @return Y an n vector simulated gaussian (not centered or scaled)
#' @return Y_std an n vector simulated gaussian centered and scaled
#' @return beta a p vector of effects
#' @return sigma
#' @return sigma_std
#' @return mean_corX mean of correlations of effect variables (lower triangular entries of correlation matrix of effect variables)
#'
#' @importFrom stats rnorm
#' 
#' @export
#' 

# Simulate regression data matrix X
simulate_regression_data <- function (n, p) {
  
  # check input: [CHECK atomic vs. scalar]
  if (!(is.scalar(n) & all(n >= 2)))
    stop("Input argument \"n\" should be 2 or more")
  if (!(is.scalar(p) & all(p >= 2)))
    stop("Input argument \"p\" should be 2 or more")
  
  # simulate data matrix X: 
  X <- matrix(rnorm(n*p),n,p)
  return(list(X = X))
}

#' @importFrom stats sd
#' @importFrom stats cor
#
# Simulate a normal y from data matrix X [MODIFIED FROM simulate3.R]
sim_gaussian <- function(X, pve, s) {
  n = dim(X)[1]
  p = dim(X)[2]
  
  # check if s <= p
  if (s > p) {
    stop("number of effects more than variables")
  }
  
  # generate effect variables
  beta.idx = sample(p, s)
  beta = rep(0,p)
  
  # obtain correlation structure
  if(s > 0){
    #beta.values = rnorm(s, 0, sqrt(effect_sigma))
    beta.values = rnorm(s)
    beta[beta.idx] = beta.values
  }
  if (s==1){
    mean_corX = 1
  } else {
    effectX = X[,beta.idx]
    corX = cor(effectX)
    mean_corX = mean(abs(corX[lower.tri(corX)]))
  }
  
  # generate non-standardized sim.y and standardized Y 
  if(s==0){
    sigma = 1
    sim.y = rnorm(n, 0, 1)
    Y = (sim.y - mean(sim.y))/sd(sim.y)
  } else {
    y = X %*% beta
    sigma = sqrt(var(y)*(1-pve)/pve)
    epsilon = rnorm(n, mean = 0, sd = sigma)
    sim.y = y + epsilon
    Y = (sim.y - mean(sim.y))/sd(sim.y)
  }
  
  return(list(Y_std = Y, Y = sim.y, 
              sigma = sigma, sigma_std = sigma/sd(sim.y),
              beta = beta, mean_corX = mean_corX))
}


# Master function; Generates X and y with specified parameters
simulate_data <- function(n, p, pve, s) {
  sim_X <- simulate_regression_data(n, p)
  sim_y <- sim_gaussian(sim_X$X, pve, s)
  
  lst <- list(X = sim_X$X, y = sim_y$Y, sigma = sim_y$sigma, 
             beta = sim_y$beta, pve = pve, s = s)
  return(lst)
}
