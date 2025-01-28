# Original code was created by Eiko Fried, Faidra Papanikolaou & Sascha Epskamp. 
# This version of code was (heavily) edited in 29.11.2024 for reanalysis of Fried(2022)
# Mental Health and Social Contact During the COVID-19 Pandemic: 
# An Ecological Momentary Assessment Study.
# For the original code and data, see https://osf.io/kj5zh

# --------------- 1. Loading packages & Data ------------------------------
# List of required packages
setwd("C:/LocalData/lintusak/TriTheoryTests/")
required_packages <- c(
  "Matrix", "fastmatrix", "expm", "rstan",
  "qgraph", "tidyverse", 
  "ggplot2", "rstantools", "loo", "bayesplot",
  "shinystan", "cmdstanr"
)

# Function to check and install missing packages
install_if_missing <- function(packages) {
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg)
    }
    library(pkg, character.only = TRUE)
  }
}

install_if_missing(required_packages)

# --------------- 2. Estimation -----------------------------

# load 
load(file.path(getwd(), "/Fried_2022 data/clean_network.RData")); gc()
Data5b <- Data2

# Variables to investigate:
vars <- c("Q1","Q2","Q3","Q4","Q5","Q6","Q7","Q8","Q9",
          "Q10","Q11","Q12","Q13","Q14","Q15","Q16","Q17","Q18")

# Labels:
varLabs <- c("Relax","Irritable","Worry","Nervous","Future","Anhedonia",
             "Tired","Hungry","Alone","Angry","Social_offline","Social_online","Music",
             "Procrastinate","Outdoors","C19_occupied","C19_worry","Home")
names(Data5b)[names(Data5b) %in% vars] <- varLabs

varLabs <- c("Relax", "Worry",  
             "Nervous")
varLabs2 <- c("Relax", "Worry",  
              "Nervous", "Tired", "Hungry",
              "Alone", "Angry")


# Lagged variables
Data5b <- Data5b %>%
  group_by(id) %>%
  mutate(across(all_of(varLabs2), 
                ~ lag(.x), 
                .names = "{.col}_lag")) %>%
  ungroup()

# Exclude double beeps:
Data5b <- Data5b %>% group_by(id) %>%
  filter(beep != lag(beep)) %>% ungroup

# Exclude the IDs with lower than 50 obs count:
excluded_ids <- Data5b %>% group_by(id) %>% summarize(n = length(id)) %>% arrange(n, id) %>%
  filter(n < 50) %>% select(id) %>% unlist()

# Exclude the IDs with more than 7 days worth (half) of missing values in any variable
excluded_ids <- base::union(excluded_ids,
                            Data5b %>%
                              group_by(id) %>%
                              summarise(across(all_of(varLabs), ~ sum(is.na(.)))) %>%
                              filter(if_any(all_of(varLabs), ~ . > 28)) %>%
                              pull(id))
# Filter data.
Data5b <- Data5b %>% filter(!(id %in% excluded_ids))



# Truncate to 4.
Data5b <- Data5b %>% mutate(across(all_of(varLabs2), ~ pmin(.x, 4)))


# STAN Bayesian multilevel analysis. -----------------

  # Data:
{
  subjects <- unique(Data5b$id)
S <- length(subjects)
N <- nrow(Data5b)
K <- length(varLabs)  
subject <- Data5b$id
X <- matrix(nrow = N, ncol = K)
X[,1] <- Data5b$Relax
X[,2] <- Data5b$Worry
X[,3] <- Data5b$Nervous
X[is.na(X)] <- 1
start <- integer(S)
end <- integer(S)
beep <- Data5b$beep + 1
nbeeps <- length(unique(beep))
for (s in seq_along(subjects)) {
  subj_indices <- c()
  subj_indices <- which(Data5b$id == subjects[s])
  start[s] <- min(subj_indices)
  end[s] <- max(subj_indices)
}
missing_mask <- matrix(0, nrow = N, ncol = K)
missing_mask[, 1] <- as.integer(is.na(Data5b$Relax))
missing_mask[, 2] <- as.integer(is.na(Data5b$Worry))
missing_mask[, 3] <- as.integer(is.na(Data5b$Nervous))

# Vectorized form:
X_vector <- as.vector(X) # X as a vector so that vectorized operations can be used in R Stan.
X_k <- rep(1:K, each = N)
observed_indices = which(as.vector(missing_mask)==0) # For selecting only the obs vals.
N_obs = length(observed_indices)
X_observed = X_vector[observed_indices]
X_k_observed = X_k[observed_indices]

obsval_slice_start <- integer(0L)
obsval_slice_end <- integer(0L)
for (k in 1:K) {
  var_indices <- c()
  var_indices <- which(X_k_observed == k)
  obsval_slice_start[k] <- min(var_indices)
  obsval_slice_end[k] <- max(var_indices)
}


# Cutpoint settings:
cutpoint_prior_locations <- rowMeans(qnorm((sapply(Data5b %>% select(Relax, Worry, Nervous),
                                                   function(x) cumsum(prop.table(table(na.omit(x))))))))
cutpoint_prior_locations <- cutpoint_prior_locations[-length(cutpoint_prior_locations)]
cutpoint_count <- length(cutpoint_prior_locations)

# Settings for the X_star_zero, the t = 0 unobserved variables.
mu0 = rep(0, times = K)
Sigma0 = matrix(rep(0.5, times = K^2), ncol = K, nrow = K)
diag(Sigma0) <- 1
  }

# Prepare data list for Stan
stan_data <- list(
  S = S,
  K = K,
  N = N,
  subject = subject,
  X = X,
  missing_mask = missing_mask,
  start = start,
  end = end,
  cutpoint_prior_locations = cutpoint_prior_locations,
  cutpoint_count = cutpoint_count,
  mu0 = mu0,
  Sigma0 = Sigma0,
  beep = beep,
  nbeeps = nbeeps,
  X_vector = X_vector,
  X_k = X_k,
  observed_indices = observed_indices,
  N_obs = N_obs,
  X_observed = X_observed,
  obsval_slice_start = obsval_slice_start,
  obsval_slice_end = obsval_slice_end
)


  ## First analysis: 3 symptom Network ----
# Set of variables of interest, which we can set to pars argument, or later on extract.
inference_vars <- c(
  sapply(1:K, function(x) paste0("A[",x,",",1:K, "]")),
  sapply(1:K, function(x) paste0("B[",x,",",1:K, "]")),
  sapply(1:K, function(x) paste0("Omega[",x,",",1:K, "]")),
  sapply(1:K, function(x) paste0("cutpoints[",x,",",1:cutpoint_count, "]")),
  sapply(1:K, function(x) paste0("time_of_day_intercept[",x,",",1:nbeeps, "]")),
  sapply(1:K, function(x) paste0("X_star_innovation_sd[",x,"]"))
)
inference_vars_regex <- c("A", "Omega", "cutpoints","time_of_day_intercept")
# Run the Network model
stan_model_Net <- stan_model(file = "BayesianOrderedVAR.stan"); gc()
fit_Net <- sampling(stan_model_Net, data = stan_data, 
                iter = 1000, chains = 1, cores = 1, 
                control = list(adapt_delta = 0.99),
                init = function() {
                  list(
                    A = diag(rep(0, stan_data$K)),
                    L_corr = diag(rep(1, stan_data$K)),
                    subject_intercept_raw = matrix(0, stan_data$S, stan_data$K),
                    subject_intercept_sd = rep(1, stan_data$K),
                    X_star_innovation = matrix(0, stan_data$K, stan_data$N),
                    X_star_zero = matrix(rep(0, K*S), nrow = K, ncol = S),
                    cutpoints = matrix(rep(cutpoint_prior_locations,times = stan_data$K), ncol = stan_data$K, byrow = T),
                    time_of_day_intercept = replicate(stan_data$nbeeps, rep(0, stan_data$K), simplify = FALSE)
                  )},
                pars = inference_vars_regex
                ); gc()

# 2. Run MCMC with cmdstanr
mod_Net <- cmdstan_model("BayesianOrderedVAR_alpha.stan")
fit_Net <- mod_Net$sample(
  data = stan_data,
  seed = 123,                 # or your preferred seed
  refresh = 100,
  chains = 1,
  parallel_chains = 1,        # set >1 if you want parallel MCMC
  iter_warmup = 500,          # e.g. half of 1000
  iter_sampling = 500,        # total 1000 draws
  adapt_delta = 0.99, 
  init = function(chain_id) {
    list(
      A = diag(rep(0, stan_data$K)),
      L_corr = diag(rep(1, stan_data$K)),
      subject_intercept_raw = matrix(0, stan_data$S, stan_data$K),
      subject_intercept_sd = rep(1, stan_data$S),
      X_star_innovation = matrix(0, stan_data$K, stan_data$N),
      X_star_zero = matrix(0, stan_data$K, stan_data$S),
      cutpoints = matrix(
        rep(stan_data$cutpoint_prior_locations,
            times = stan_data$K),
        ncol  = stan_data$K,
        byrow = TRUE
      ),
      time_of_day_intercept = replicate(
        stan_data$nbeeps,
        rep(0, stan_data$K),
        simplify = FALSE
      )
    )
  }
)

# Possible conversion to matrix and saving.
fit_Net <- as.matrix(fit_Net); gc() # Rewrite, so that RAM is not occupied.
saveRDS(fit_Net, file = "Datas/BayesianVAR_3_symptoms_28_01.RDS"); gc()

# Diagnostics
plotnams <- c("A", "Omega","B", "cutpoints",
              "time_of_day_intercept")
for(i in 1:length(plotnams)){
  dev.new()
  print(  mcmc_trace(fit_Net, regex_pars  = plotnams[i] ))
  }; gc()
# Pairwise plots, for multicollinearit:
pairplotter <- function(var1,var2,fit = fit_Net) {
  mcmc_pairs(fit_Net, regex_pars = c(var1,var2))
}

# Extract parameters
A_mean <- matrix(0, ncol = K, nrow = K)
for(i in 1:K) for(j in 1:K){
  A_mean[i,j] <- mean(as.vector(fit_Net[,paste0("A[",i,",",j,"]")]))
  }
Omega_mean <- matrix(0, ncol = K, nrow = K)
for(i in 1:K) for(j in 1:K){
  Omega_mean[i,j] <- mean(as.vector(fit_Net[,paste0("Omega[",i,",",j,"]")]))
}

# Plot the Network:
par(mfrow=c(1,2))
qgraph(A_mean);qgraph(Omega_mean)
par(mfrow=c(1,1))



  ## Second analysis: More symptoms -----

# Data: OVERWRITTEN to not overflow the environment.
{
  subjects <- unique(Data5b$id)
  S <- length(subjects)
  N <- nrow(Data5b)
  K <- length(varLabs2)  
  subject <- Data5b$id
  X <- matrix(nrow = N, ncol = K)
  X[,1] <- Data5b$Relax
  X[,2] <- Data5b$Worry
  X[,3] <- Data5b$Nervous
  X[,4] <- Data5b$Tired
  X[,5] <- Data5b$Hungry
  X[,6] <- Data5b$Alone
  X[,7] <- Data5b$Angry
  X[is.na(X)] <- 1
  start <- integer(S)
  end <- integer(S)
  beep <- Data5b$beep + 1
  nbeeps <- length(unique(beep))
  for (s in seq_along(subjects)) {
    subj_indices <- c()
    subj_indices <- which(Data5b$id == subjects[s])
    start[s] <- min(subj_indices)
    end[s] <- max(subj_indices)
  }
  missing_mask <- matrix(0, nrow = N, ncol = K)
  missing_mask[, 1] <- as.integer(is.na(Data5b$Relax))
  missing_mask[, 2] <- as.integer(is.na(Data5b$Worry))
  missing_mask[, 3] <- as.integer(is.na(Data5b$Nervous))
  missing_mask[, 4] <- as.integer(is.na(Data5b$Tired))
  missing_mask[, 5] <- as.integer(is.na(Data5b$Hungry))
  missing_mask[, 6] <- as.integer(is.na(Data5b$Alone))
  missing_mask[, 7] <- as.integer(is.na(Data5b$Angry))
  
  # Vectorized form:
  X_vector <- as.vector(X) # X as a vector so that vectorized operations can be used in R Stan.
  X_k <- rep(1:K, each = N)
  observed_indices = which(as.vector(missing_mask)==0) # For selecting only the obs vals.
  N_obs = length(observed_indices)
  X_observed = X_vector[observed_indices]
  X_k_observed = X_k[observed_indices]
  
  obsval_slice_start <- integer(0L)
  obsval_slice_end <- integer(0L)
  for (k in 1:K) {
    var_indices <- c()
    var_indices <- which(X_k_observed == k)
    obsval_slice_start[k] <- min(var_indices)
    obsval_slice_end[k] <- max(var_indices)

  }
  }
  
# Cutpoint settings:
cutpoint_prior_locations <- rowMeans(qnorm((sapply(Data5b %>% select(Relax, Worry, Nervous),
                                                   function(x) cumsum(prop.table(table(na.omit(x))))))))
cutpoint_prior_locations <- cutpoint_prior_locations[-length(cutpoint_prior_locations)]
cutpoint_count <- length(cutpoint_prior_locations)

# Settings for the X_star_zero, the t = 0 unobserved variables.
mu0 = rep(0, times = K)
Sigma0 = matrix(rep(0.5, times = K^2), ncol = K, nrow = K)
diag(Sigma0) <- 1

# Prepare data list for Stan
stan_data <- list(
  S = S,
  K = K,
  N = N,
  subject = subject,
  X = X,
  missing_mask = missing_mask,
  start = start,
  end = end,
  cutpoint_prior_locations = cutpoint_prior_locations,
  cutpoint_count = cutpoint_count,
  mu0 = mu0,
  Sigma0 = Sigma0,
  beep = beep,
  nbeeps = nbeeps,
  X_vector = X_vector,
  X_k = X_k,
  observed_indices = observed_indices,
  N_obs = N_obs,
  X_observed = X_observed,
  obsval_slice_start = obsval_slice_start,
  obsval_slice_end = obsval_slice_end
)

# Set of variables of interest, which we can set to pars argument, or later on extract.
inference_vars <- c(
  sapply(1:K, function(x) paste0("A[",x,",",1:K, "]")),
# sapply(1:K, function(x) paste0("A[",x,",",1:K, "]")),
  sapply(1:K, function(x) paste0("Omega[",x,",",1:K, "]")),
  sapply(1:K, function(x) paste0("cutpoints[",x,",",1:cutpoint_count, "]")),
  sapply(1:K, function(x) paste0("time_of_day_intercept[",x,",",1:nbeeps, "]"))
)
# Nuisance set
nuisance_vars <- c(sapply(1:K, function(x) paste0("X_star_innovation[",x,",",1:N, "]")),
                   sapply(1:K, function(x) paste0("subject_innovation_sd[",x,",",1:S, "]")),
                   sapply(1:K, function(x) paste0("subject_intercept[",1:S,",",x, "]")),
                   sapply(1:K, function(x) paste0("subject_intercept_sd[",1:S,",",x, "]")),
                   sapply(1:K, function(x) paste0("subject_intercept_raw[",1:S,",",x, "]")),
                   sapply(1:K, function(x) paste0("X_star_zero[",x,",",1:S, "]")),
                   sapply(1:K, function(x) paste0("X_star[",x,",",1:N, "]")),
                   "lp__")

# Run the Network model
stan_model_Net <- stan_model(file = "BayesianOrderedVAR.stan")
fit_Net_7 <- sampling(stan_model_Net, data = stan_data, 
                    iter = 4000, chains = 8, cores = 8, 
                    control = list(adapt_delta = 0.95),
                    init = function() {
                      list(
                        A = diag(rep(0.1, stan_data$K)),
                        L_corr = diag(rep(1, stan_data$K)),
                        subject_intercept_raw = matrix(0, stan_data$S, stan_data$K),
                        subject_intercept_sd = rep(1, stan_data$K),
                        X_star_innovation = matrix(0, stan_data$K, stan_data$N),
                        X_star_zero = matrix(rep(0, K*S), nrow = K, ncol = S),
                        subject_innovation_sd = matrix(1, stan_data$K, stan_data$S),
                        cutpoints = replicate(stan_data$K, seq(-1, 1, length.out = stan_data$cutpoint_count), simplify = FALSE),
                        time_of_day_intercept = replicate(stan_data$nbeeps, rep(0, stan_data$K), simplify = FALSE)
                      )}); gc()

# Split the draws to inference set, likelihood set and nuisance set
fit_Net_7 <- as.matrix(fit_Net_7); gc() # Rewrite, so that RAM is not occupied.
saveRDS(as.matrix(fit_Net_7)[,dimnames(fit_Net_7)$parameters %in% nuisance_vars], "Datas/BayesOrderedVAR_FIT_7_NUISANCE_POSTERIOR.RDS", compress = T); gc()
# Likelihood set
saveRDS(fit_Net_7[, grepl("log_lik_7", dimnames(fit_Net_7)$parameters) ], "Datas/BayesOrderedVAR_FIT_7_LIKELIHOOD_POSTERIOR.RDS", compress = T); gc()
# Inference set:
fit_Net_7 <- fit_Net_7[, dimnames(fit_Net_7)$parameters %in% inference_vars ]; gc()
saveRDS(fit_Net_7, file ="Datas/BayesOrderedVAR_FIT_7_INFERENCE_POSTERIOR_27_01.RDS"); gc()
# Diagnostics
plotnams <- c("A","B", "Omega", "cutpoints",
              "time_of_day_intercept")
for(i in 1:length(plotnams)){
  dev.new()
  print(  mcmc_trace(fit_Net_7, regex_pars  = plotnams[i] ))
}

# Posterior plots
for(i in 1:length(plotnams)){
  dev.new()
  print(  mcmc_intervals(fit_Net_7, regex_pars  = plotnams[i] ))
}


# WAIC and LOO


# Extract parameters
posterior_samples <- rstan::extract(fit_Net_7, pars = "A", permuted = TRUE)

A_mean <- matrix(0, ncol = K, nrow = K)
for(i in 1:K) for(j in 1:K){
  A_mean[i,j] <- mean(as.vector(fit_Net_7[,paste0("A[",i,",",j,"]")]))
}
Omega_mean <- matrix(0, ncol = K, nrow = K)
for(i in 1:K) for(j in 1:K){
  Omega_mean[i,j] <- mean(as.vector(fit_Net_7[,paste0("Omega[",i,",",j,"]")]))
}

# Plot the Network:
par(mfrow=c(1,2))
qgraph(A_mean);qgraph(Omega_mean)
par(mfrow=c(1,1))




  ## DCF ------
# Run the DCF model
stan_model_DCF <- stan_model(file = "BayesianOrderedDCF.stan")
fit_dcf <- sampling(stan_model_DCF, data = stan_data, 
                iter = 2000, chains = 4, cores = 4, 
                control = list(adapt_delta = 0.80) )
saveRDS(fit_dcf, "Datas/BayesOrderedDCF_FIT.RDS", compress = T); gc()
# Diagnostics
print(fit_dcf, pars = c("psi", "Lambda", "cutpoints", paste0("eta[",1:10,"]"), 
                        paste0("eta_zero[",1:5,"]"),
                        "subject_intercept"))
dev.new(noRStudioGD = T)
traceplot(fit_dcf)
dev.new(noRstudioGD = T)
traceplot(fit_dcf, pars = "cutpoints")
dev.new(noRstudioGD = T)
traceplot(fit_dcf, pars = c("time_of_day_intercept", "eta_innovation"))
check_hmc_diagnostics(fit_dcf)

# WAIC and LOO
log_lik_matrix_dcf <- extract_log_lik(fit_dcf, parameter_name = "log_lik", merge_chains = FALSE)
waic_results_dcf <- waic(log_lik_matrix_dcf)
print(waic_results_dcf)
loo_results_dcf <- loo(log_lik_matrix_dcf)
print(loo_results_dcf)

