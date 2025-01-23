# Original code was created by Eiko Fried, Faidra Papanikolaou & Sascha Epskamp. 
# This version of code was (heavily) edited in 29.11.2024 for reanalysis of Fried(2022)
# Mental Health and Social Contact During the COVID-19 Pandemic: 
# An Ecological Momentary Assessment Study.
# For the original code and data, see https://osf.io/kj5zh

# --------------- 1. Loading packages & Data ------------------------------
# List of required packages
required_packages <- c(
  "Matrix", "fastmatrix", "BVAR", 
  "brms", "expm", "rstan",
  "qgraph", "tidyverse", 
  "ggplot2", "blavaan", "lavaan",
  "rstantools", "loo", "bayesplot"
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

# Exclude the IDs with lower than 50 obs count:
excluded_ids <- Data5b %>% group_by(id) %>% summarize(n = length(id)) %>% arrange(n, id) %>%
  filter(n < 50) %>% select(id) %>% unlist()
Data5b <- Data5b %>% filter(!(id %in% excluded_ids))

# Exclude double beeps:
Data5b <- Data5b %>% group_by(id) %>%
  filter(beep != lag(beep)) %>% ungroup

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
M = sum(missing_mask == 1)
missing_positions <- which(X == 0, arr.ind = TRUE)
sorted_positions <- missing_positions[order(missing_positions[, 1]), ]
missing_idx <- sorted_positions[,1,drop=T]
missing_var <- sorted_positions[,2,drop=T]
cutpoint_prior_locations <- 1:max(Data5b$Relax, na.rm = T) - mean(1:max(Data5b$Relax, na.rm = T))
cutpoint_count <- length(cutpoint_prior_locations)
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
  M = M,
  missing_idx = missing_idx,
  missing_var = missing_var,
  cutpoint_prior_locations = cutpoint_prior_locations,
  cutpoint_count = cutpoint_count,
  mu0 = mu0,
  Sigma0 = Sigma0,
  beep = beep,
  nbeeps = nbeeps
)


  ## Network ----
# Run the Network model
stan_model_Net <- stan_model(file = "BayesianOrderedVAR.stan")
fit_Net <- sampling(stan_model_Net, data = stan_data, 
                iter = 1000, chains = 4, cores = 4, 
                control = list(adapt_delta = 0.80),
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
                  )}
                )
saveRDS(fit_Net, "Datas/BayesOrderedVAR_FIT.RDS", compress = T); gc()
# Diagnostics
print(fit_Net, pars = c("A","Omega", "subject_intercept_sd"))
dev.new(noRStudioGD = T)
traceplot(fit_Net)
traceplot(fit_Net, pars = paste0("X_star[1,",1:9,"]"))
traceplot(fit_Net, pars = "cutpoints")
traceplot(fit_Net, pars = "Omega")
check_hmc_diagnostics(fit_Net)


# WAIC and LOO
log_lik_matrix_Net <- extract_log_lik(fit, parameter_name = "log_lik", merge_chains = FALSE)
waic_results_Net <- waic(log_lik_matrix_Net)
print(waic_results_Net)
loo_results_Net <- loo(log_lik_matrix_Net)
print(loo_results_Net)




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

