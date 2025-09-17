# Source simulation
source("./second publication/intro_simulations.R")

# Stan model. ----
library(rstan)

# Define the data list.
stan_data <- list(
  # Y       = scale(Y_case1), # Diagonal VAR(10) noise.
  Y = scale(Y_case2), # Diagonal VAR(1) noise.
  # Y = scale(Y_case3), # White noise.
  # Y = scale(Y_case4), # VMA(1) noise.
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
  iter = 4000,
  init = "random",
  data = stan_data, 
  include = T, 
  pars = c("Sigma", "A", "A_eps", "c")
)
