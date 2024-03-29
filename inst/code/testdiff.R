# Test differences in mr.ash and mr.ash.alpha outputs

# setwd("/project2/mstephens/lwang19/mr.ash")
set.seed(1)

# Load data
data = readRDS("../datafiles/sim_data_1.rds") # 20 simulated datasets in data folder
X = data$X
y = data$y

# Fit using mr.ash.alpha.
fit.alpha <- mr.ash.alpha::mr.ash(X, y)
post.alpha <- mr.ash.alpha::get.full.posterior(fit.alpha)
phi.alpha <- post.alpha$phi

# Fit using mr.ash.
fit <- mr.ash::mr.ash(X,y)
post <- mr.ash::get.full.posterior(fit)
phi <- post$phi

# Identified differences through manual testing: 
identical(fit, fit.alpha) # FALSE
identical(phi, phi.alpha) # FALSE

## Class differences: 
alpha.b = fit.alpha$beta
class(alpha.b) # matrix
class(fit.alpha$pi) # matrix

fit.b = fit$beta
class(fit.b)  # numeric
class(fit$pi) # numeric

fit.alpha$beta <- drop(fit.alpha$beta)
fit.alpha$pi <- drop(fit.alpha$pi)
identical(fit, fit.alpha) # TRUE

## numerical differences (due to get.full.posterior):
phi.alpha[c(1:5), c(1:5)]
phi[c(1:5), c(1:5)]

phi.alpha[1,1]
phi[1,1]


