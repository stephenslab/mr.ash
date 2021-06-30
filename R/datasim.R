#' @title Simulate Regression Data
#'
#' @description Simulate non-sparse data matrix X (and a corresponding normal y?)
#' 
#' @param n number of samples 
#' @param p number of variables
#' @param s number of non-zero effects
#' @param se standard deviation of the residual
#' @param pve proportion of variance in Y explained by X
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

# Original stand-in function
simulate_regression_data <- function (n, p) {
  
  # simulate data matrix X: 
  X <- matrix(rnorm(n*p),n,p)
  return(list(X = X))
}


# Simulate the gene matrix X. 
simulate_gene_data <- function(n, p) {
  
  # check input: [CHECK atomic vs. scalar]
  if (!(is.atomic(n) & all(n >= 2)))
    stop("Input argument \"n\" should be 2 or more")
  if (!(is.atomic(p) & all(p >= 2)))
    stop("Input argument \"p\" should be 2 or more")
  
  # simulate gene data set [BORROWED FROM varbvs.qtl.R]
  maf <- 0.05 + 0.45*runif(p)
  X   <- (runif(n*p) < maf) +
         (runif(n*p) < maf)
  X   <- matrix(as.double(X),n,p,byrow = TRUE)
  
  return(X)
}

# Simulate a normal y from data matrix X [MODIFIED FROM simulate3.R]
sim_gaussian <- function(X, pve, s) {
  n = dim(X)[1]
  p = dim(X)[2]
  
  # generate effect variables
  beta.idx = sample(p, effect_num)
  beta = rep(0,p)
  
  # obtain correlation structure
  if(effect_num > 0){
    beta.values = rnorm(effect_num, 0, sqrt(effect_sigma))
    beta[beta.idx] = beta.values
  }
  if (effect_num==1){
    mean_corX = 1
  } else {
    effectX = X[,beta.idx]
    corX = cor(effectX)
    mean_corX = mean(abs(corX[lower.tri(corX)]))
  }
  
  # generate non-standardized sim.y and standardized Y 
  if(effect_num==0){
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
              beta = beta, mean_corX = mean_corX)
}