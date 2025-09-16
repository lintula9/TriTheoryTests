# Source simulation
source("./second publication/intro_simulations.R")

# Stan model. ----
library(rstan)

# Define the data list.
stan_data <- list(
  Y       = scale(Y_case1),
  # Y_case2 = scale(Y_case2),
  # Y_case3 = scale(Y_case3),
  # Y_case4 = scale(Y_case4),
  N       = dim(Y_case1)[1],
  p       = dim(Y_case1)[2]
)


## Run the Stan model. ----
library(rstan)
stan_model       <- stan(
  file = list.files(path = ".", pattern = "intro_simulations.stan", 
                    recursive = TRUE, full.names = TRUE), 
  model_name = "basic_model",
  cores = 4,
  init = "random",
  data = stan_data,
)
