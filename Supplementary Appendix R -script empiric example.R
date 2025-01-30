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
  "shinystan", "cmdstanr", "posterior"
)

# Function to check and install missing packages
for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg)
    }
    library(pkg, character.only = TRUE)
  }

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

# Run MCMC with cmdstanr
mod_Net <- cmdstan_model("BayesianOrderedVAR_alpha.stan")
fit_Net <- mod_Net$sample(
  data = stan_data,
  seed = 123,                 # or your preferred seed
  refresh = 100,
  chains = 8,
  parallel_chains = 8, 
  num_cores = 8,       # set >1 if you want parallel MCMC
  iter_warmup = 2000,          # e.g. half of 1000
  iter_sampling = 2000,        # total 1000 draws
  adapt_delta = 0.95,
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

# Read data
if(F){
fit_Net <- as_cmdstan_fit(files = paste0("Datas/fit_net_3var_output-202501290753-",1:8,"-5996ed.csv") ); gc()
draws_df <- as_draws_df(fit_Net$draws(variables = c(inference_vars_regex, "lp__") )); gc(); rm(fit_Net); gc()
}
# Diagnostics and posterior distribution marignal plots. 
plotnams <- c("A", "Omega", "cutpoints","time_of_day_intercept");pdf(file = paste0("Datas/Bayespots_",format(Sys.time(), "%Y-%m-%d"), ".pdf"))
for(i in plotnams){print(  mcmc_trace(draws_df[ , grep(i, names(draws_df))]) );print(  mcmc_areas_ridges(draws_df[ , grep(i, names(draws_df))]) );
  }; gc(); dev.off()
# Extract posterior means as the VAR parameters:
A <- matrix( unlist(colMeans(draws_df[,grep("A", names(draws_df))])), ncol = K, nrow = K);
Z <- matrix( unlist(colMeans(draws_df[map_index,grep("Omega", names(draws_df))])), ncol = K, nrow = K)

par(mfrow=c(1,2)); pdf(file = paste0("Datas/Bayes_MAP_AZestimates",format(Sys.time(), "%Y-%m-%d"), ".pdf"))
qgraph(input = A);qgraph(input = Z);dev.off();par(mfrow=c(1,1));gc()

# Find the closest indistinguishable model, and compare
source("Supplementary Appendix R -script civ_find.R")
sqrt(mean(c(c(civ_find(A,Z)$A - A)^2, c(civ_find(A,Z)$Z[lower.tri(Z, diag = T)] - Z[lower.tri(Z, diag = T)])^2)))

# Randomly select 1000 posterior draws and obtain discrepancies
indices <- sample(1:nrow(draws_df),size = 1000L, replace = F)
estimated_var_samples <- draws_df[i,grep("A", names(draws_df)) | grep("Omega", names(draws_df))]; As <- grep("A", names(Alphas_Omegas)); Os <- grep("Omega", names(Alphas_Omegas)); 
civ_var_samples <- pbapply::pbsapply(1:length(indices),
                  FUN = function(i){
                    A_temp <- matrix( estimated_var_samples[i,As], ncol = K, nrow = K);
                    Z_temp <- matrix( draws_df[i,Os], ncol = K, nrow = K);
                    closest <- civ_find(A_temp,Z_temp); A_ <- fastmatrix::vec(closest$A); Z_ <- fastmatrix::vec(closest$Z); theta_0 <- c(A_,Z_)
                    return(theta_0)
                    }, simplify = "matrix"); gc()
civ_var_samples <- t(civ_var_samples)
discrepancy_samples <- estimated_var_samples - civ_var_samples
S_discrepancy <- cov(discrepancy_samples)

# Compute Mahalanobis squared distance:
civ_var <- civ_find(A,Z)
discrepancy <- c(fastmatrix::vec(A - civ_var$A), fastmatrix::vech( Z - civ_var$Z ))
D_squared <- t(discrepancy) %*% solve(S_discrepancy) %*% discrepancy

# Scenario 1.: Set discrepany = 0 as the null hypothesis, normal theory based inference:
pchisq(q = D_squared, df = length(discrepancy))

# Scenario 2.: Set discrepancy > epsilon as the null hypothesis, normal theory based inference:
# Accept, say, RMSEA < 0.08. Find non-centrality parameter, if RMSEA = 0.08, alpha = 0.05.
# RMSEA function
ncpar_find <- function(ncppar, 
                       chisq = D_squared, 
                       df = length(discrepancy), 
                       n = nrow(Data5b)) {
  sqrt(( pchisq(chisq, df = df, ncp = ncppar) - df/(n-1))/(df)) - 0.08};
ncpar <- optim(par = list(ncppar = 0), fn = ncpar_find, gr = numDeriv::grad(ncpar_find))
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

# Run MCMC with cmdstanr
mod_Net_7 <- cmdstan_model("BayesianOrderedVAR_alpha.stan")
fit_Net_7 <- mod_Net_7$sample(
  data = stan_data,
  seed = 123,                 # or your preferred seed
  refresh = 100,
  chains = 8,
  parallel_chains = 8, 
  num_cores = 8,       # set >1 if you want parallel MCMC
  iter_warmup = 2000,          # e.g. half of 1000
  iter_sampling = 2000,        # total 1000 draws
  adapt_delta = 0.95,
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
        nrow  = stan_data$K,
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
