# Simulating test datasets

setwd("/project2/mstephens/lwang19/mr.ash/")
source("R/datasim.R")

n <- sample.int(100, size = 20)
p <- sample.int(1000, size = 20)
pve <- runif(20)
s <- sample.int(20, size = 20)

for (i in 1:length(n)) {
  sim_X <- simulate_regression_data(n[i], p[i])$X
  sim_y <- sim_gaussian(sim_X, pve, s[i])$Y
  saveRDS(list(X = sim_X, y = sim_y), file = paste0("data/sim_data_", i, ".rds"))
}

# Script generated 20 simulated datasets in mr.ash/R/data for testing purposes.
# Simulated data sets are saved on midway. Only sim_data_1 is uploaded to Github under
# mr.ash/inst/datafiles. 