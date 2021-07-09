# Simulating test datasets

setwd("/project2/mstephens/lwang19/mr.ash/")
source("R/datasim.R")
source("R/misc2.R")

set.seed(1)
n <- sample(2:400, size = 20)
p <- sample(2000, size = 20)
pve <- runif(20)
s <- sample(20, size = 20, replace = T)

for (i in 1:length(n)) {
  sim_data = simulate_data(n[i], p[i], pve[i], s[i])
  saveRDS(sim_data, file = paste0("data/sim_data_", i, ".rds"))
}

# Script generated 20 simulated datasets in mr.ash/R/data for testing purposes.
# Simulated data sets are saved on midway. Only sim_data_1 is uploaded to Github under
# mr.ash/inst/datafiles. 