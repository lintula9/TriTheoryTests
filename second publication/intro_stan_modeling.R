# Source simulation and options.
source("./second publication/intro_simulations.R")
source("./stan_models/stan_options.R")

# Run the Stan model. ----
stan_file        <- list.files(path = ".", pattern = "intro_simulations.stan", 
                               recursive = TRUE, full.names = TRUE)
if(F){ # If modeling is needed to be reran.
  stan_model       <- cmdstan_model(stan_file, dir = "./Datas/")
  stan_fit <- stan_model$sample(
    data = stan_data,
    seed = 123,          # optional: set seed for reproducibility
    chains = 4,
    parallel_chains = 4, # same as cores in rstan
    iter_warmup = 1000,  # rstan's iter = 4000 means 2000 warmup + 2000 sampling (default)
    iter_sampling = 3000); gc()
  stan_fit$sampler_diagnostics()
  }
# Done?.
rm(list = ls()); gc()

## Posterior summaries. ----
# Fetch model. 
intro_csvs <- list.files("./Datas", pattern = "^intro_simulations-.*\\.csv$", full.names = TRUE)
stan_fit   <- cmdstanr::as_cmdstan_fit(intro_csvs) ; gc()

# Arbitrary statistic.
stan_fit$summary("A[1,1]", pr_lt_half = ~ mean(. <= 0.5)) ; gc()
# Plots.
mcmc_hist(stan_fit$draws("A")) ; gc()
mcmc_hist(stan_fit$draws("A_eps")) ; gc()
