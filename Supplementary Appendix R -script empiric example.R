# Original code was created by Eiko Fried, Faidra Papanikolaou & Sascha Epskamp. 
# This version of code was (heavily) edited in 29.11.2024 for reanalysis of Fried(2022)
# Mental Health and Social Contact During the COVID-19 Pandemic: 
# An Ecological Momentary Assessment Study.
# For the original code and data, see https://osf.io/kj5zh

# --------------- 1. Loading packages & Data ------------------------------
# List of required packages
required_packages <- c(
  "Matrix", "fastmatrix", "expm", "rstan",
  "qgraph", "tidyverse", 
  "ggplot2", "rstantools", "bayesplot", "cmdstanr", "posterior",
  "viridisLite", "dplyr")

# Function to check and install missing packages
for (pkg in required_packages) {
      if(!requireNamespace(pkg, quietly = T)) {install.packages(pkg, dependencies = T)}
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
# cutpoint_prior_locations <- rowMeans(qnorm((sapply(Data5b %>% select(Relax, Worry, Nervous),
#                                                    function(x) cumsum(prop.table(table(na.omit(x))))))))
# cutpoint_prior_locations <- cutpoint_prior_locations[-length(cutpoint_prior_locations)]
# Custom adjustment to cutpoint priors, based on the model analyses. This is to identify the scales.
cutpoint_prior_locations <- c(-0.6, 1, 2.2)
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


  #3. First analysis: 3 symptom Network ----
# Set of variables of interest, which we can set to pars argument, or later on extract.
inference_vars_regex <- c("A", "Omega", "cutpoints","time_of_day_intercept")
inference_vars_regex_alpha <- c("A_effective","A","psi","Lambda", "Omega","L_", "cutpoints", "time_of_day_effect",
                                "ref_time_of_day_effect", "specific_time_of_day_effect")
# Run MCMC with cmdstanr
mod_Net <- cmdstan_model("BayesianOrderedVAR_alpha.stan")
nchains = 8
fit_Net <- mod_Net$sample(
  data = stan_data,
  seed = 123,                 # or your preferred seed
  refresh = 100,
  chains = nchains,
  parallel_chains = nchains, 
  iter_warmup = 2000,          
  iter_sampling = 2000,        
  adapt_delta = 0.95,
  init = function(chain_id) {
    list(
      A = diag(rep(0, stan_data$K)),
      psi = 0.5,
      Lambda_rest = rep(0.5,times=K-1),
      L_rest_diag = rep(0.01, times = K-1),
      L_rest_offdiag = rep(0.01, times = (K*(K-1)%/%2)),
      L_corr = diag(rep(1, stan_data$K)),
      subject_intercept_raw = matrix(0, stan_data$S, stan_data$K),
      subject_intercept_sd = rep(1, stan_data$S),
      X_star_innovation = matrix(0, stan_data$K, stan_data$N),
      X_star_zero = matrix(0, stan_data$K, stan_data$S),
      cutpoints = matrix(rep(stan_data$cutpoint_prior_locations,
            times = stan_data$K),
        ncol  = stan_data$K,
        byrow = TRUE
      ),
      ref_time_of_day_effect = rep(0, times =nbeeps),
      specific_time_of_day_effect = 
        replicate(
        stan_data$nbeeps,
        rep(0, stan_data$K -1),
        simplify = FALSE
      )
    )
  }
); gc()
# Save output.
fit_Net$save_output_files("/LocalData/lintusak/TriTheoryTests/Datas/", basename = "3VAR_CF_relaxedpriors", timestamp = T, random = F)

# Read data
if(F){
  fit_Net    <- as_cmdstan_fit(files = paste0("Datas/3VAR_CF_relaxedpriors-202502031854-",1:nchains,".csv") ); gc()
  draws_data <- as_draws_df(fit_Net$draws(variables = c(inference_vars_regex_alpha) ), .nhcains = nchains ); gc(); rm(fit_Net); gc()
  }
# Diagnostics and posterior distribution marignal plots. 
plotnams <- inference_vars_regex_alpha;pdf(file = paste0("Datas/Bayespots_",format(Sys.time(), "%Y-%m-%d"), ".pdf"));color_scheme_set("viridis")
for(i in plotnams){print(  mcmc_trace(draws_data, regex_pars = i ));print(  mcmc_areas_ridges(draws_data[ , grep(i, names(draws_data))]) );
  }; gc(); dev.off()

# Extract posterior means as the VAR parameters:
A <- matrix( unlist(colMeans(draws_data[,grep("A", names(draws_data))])),     ncol = K, nrow = K);
Z <- matrix( unlist(colMeans(draws_data[,grep("Omega", names(draws_data))])), ncol = K, nrow = K)

# Obtain complete VAR(1) model samples:
estimated_var_samples <- draws_data[,c( grep("A", names(draws_data)) , grep("Omega", names(draws_data)) )]; As <- grep("A_", names(estimated_var_samples)); Os <- grep("Omega", names(estimated_var_samples)); 
var_samples           <- pbapply::pblapply(1:nrow(draws_data),
                                 FUN = function(i){
                                   A_temp <- matrix( unlist(estimated_var_samples[i,As]), ncol = sqrt(length(As)), nrow = sqrt(length(Os)));
                                   Z_temp <- matrix( unlist(estimated_var_samples[i,Os]), ncol = sqrt(length(As)), nrow = sqrt(length(Os)));
                                   return(list(A = A_temp, Z = Z_temp))
                                 }); gc()
# Source the methods.
source("Supplementary Appendix R -script civ_find.R")

# Figure 4 in main text ----
  # Compute parallel analysis imitation and RMSEA, for the posterior mean.
result_parallel <- civ_parallel(A, Z)

  # Compute credible intervals for eigenvalues, congruencies.
eigen_congurency <- pbapply::pblapply(var_samples, FUN = function(x){
              res  <- try(civ_parallel(x$A,x$Z))
            eigens <- t(abs(res$eigenvals))
      congruencies <- res$all_factor_congruencies[,1]
    return(list(eigens = eigens, congruencies = congruencies)) }); gc()
eigen_dat <- data.frame(Re(do.call(rbind,lapply(eigen_congurency, FUN = function(x) cbind( x$eigens, 1:11 ) ))))

upper <- as.matrix(eigen_dat %>% group_by(X4) %>% reframe( across(paste0( "X", 1:(length(eigen_dat)-1) ), ~ quantile(.x, c(.975))) ))
lower <- as.matrix(eigen_dat %>% group_by(X4) %>% reframe( across(paste0( "X", 1:(length(eigen_dat)-1) ), ~ quantile(.x, c(.025))) ))

tiff(filename = "Figure_4.tiff", 
        width = 14, 
     height   = 14, 
        units = "in", 
          res = 480)
par(mfrow     = c(2,2) )
par(mar       = c(4,4,2,2) )
matplot(t(result_parallel$eigenvals),
        type = "n",
        ylab = "Eigenvalues", 
        xlab = expression(paste("Increment in time ", Delta, "T")),
        xaxt = "n",
        main = "Eigenvalues over cross-covariance increments"); grid()
axis(side = 1, at = 1:(ncol(result_parallel$eigenvals)), labels = 0:(ncol(result_parallel$eigenvals) - 1))
for( i in 2:ncol(upper)) {
  polygon(x = c(upper[,1], rev(lower[,1])), y = c(upper[,i], rev(lower[,i])),
          col = adjustcolor(cividis(i-1), alpha.f = 0.15))
}
matplot(t(result_parallel$eigenvals), 
        type = "b",
        col  = cividis(6), add = T)

cong_dat <- data.frame(Re(do.call(rbind,lapply(eigen_congurency, FUN = function(x) cbind( x$congruencies, 1:6 ) ))))
upper_c  <- as.matrix(cong_dat %>% group_by(X2) %>% 
                        reframe( quantile(X1, 0.975) ))
lower_c  <- as.matrix(cong_dat %>% group_by(X2) %>% 
                        reframe( quantile(X1, 0.025) ))

matplot(result_parallel$all_factor_congruencies[,1], type = "n",
        ylim = c(0,1),
        ylab = "Congruency coefficient", 
        xlab = "Cross-covariance pair",
        xaxt = "n",
        main = "Congruency for cross-covariance pairs"); grid()
 # NOTE: this is omitted since it only brings clutter. polygon(x = c(upper_c[,1], rev(upper_c[,1])), y=c(upper_c[,2], rev(lower_c[,2])) )
axis(1, labels = paste0("(", 0:10,", ", 1:11,")"),
     at = 1:11 )
matplot(result_parallel$all_factor_congruencies[,1], type = "b",
        col = cividis(6), add = T )
qgraph( A, layout = "circle", 
        labels = varLabs, title = "Coefficient matrix", mar = c(5,5,5,5) )
qgraph( Z, layout = "circle", 
        labels = varLabs, title = "Innovation covariance", mar = c(5,5,5,5) )

dev.off();gc();par(mfrow = c(1,1) )

# Compute RMSEA, distribution
  # mean
plot(RMSEA_approx(A,Z,error_ratios = seq(0.001,0.1, length.out = 20),N = nrow(Data5b)))

rmseas <- pbapply::pbsapply(var_samples, FUN = function(x){
  res <- try(RMSEA_approx(x$A,x$Z, N = nrow(Data5b))$RMSEA_approximation)
  if(is.character(res)) {res <- NA; warning("Non-convergence")}
  return(res) }); gc()
saveRDS(rmseas, file = "Datas/rmseas_3vars.RDS");gc()

source("Yuan 2016.R"); Yuan_2016(  length(fastmatrix::vech(Z)) + length(A) - (1+3*ncol(A))
                                     ,4135) 
# Congruency
congruency <- pbapply::pbsapply(1:length(var_samples), FUN = function(x){
  return(try(Re(civ_parallel(var_samples[[x]]$A,var_samples[[x]]$Z)$min_factor_congruency)))
}); gc()
saveRDS(congruency, file = "Datas/congruency_3vars.RDS"); gc()

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
  cutpoint_prior_locations <- c(-0.6, 1, 2.2)
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
nchains = 8
fit_Net_7 <- mod_Net_7$sample(
  data = stan_data,
  seed = 123,                 # or your preferred seed
  refresh = 100,
  chains = nchains,
  parallel_chains = nchains, 
  iter_warmup = 2000,          # e.g. half of 1000
  iter_sampling = 2000,        # total 1000 draws
  adapt_delta = 0.95,
  init = function(chain_id) {
    list(
      A = diag(rep(0, stan_data$K)),
      psi = 0.5,
      Lambda_rest = rep(0.5,times=K-1),
      L_rest_diag = rep(0.01, times = K-1),
      L_rest_offdiag = rep(0.01, times = (K*(K-1)%/%2)),
      L_corr = diag(rep(1, stan_data$K)),
      subject_intercept_raw = matrix(0, stan_data$S, stan_data$K),
      subject_intercept_sd = rep(1, stan_data$S),
      X_star_innovation = matrix(0, stan_data$K, stan_data$N),
      X_star_zero = matrix(0, stan_data$K, stan_data$S),
      cutpoints = matrix(rep(stan_data$cutpoint_prior_locations,
                             times = stan_data$K),
                         nrow  = stan_data$K,
                         byrow = TRUE
      ),
      ref_time_of_day_effect = rep(0, times =nbeeps),
      specific_time_of_day_effect = 
        replicate(
          stan_data$nbeeps,
          rep(0, stan_data$K -1),
          simplify = FALSE
        )
    )
  }
); gc()
# Save output
fit_Net_7$save_output_files("Datas/",basename = "BVAR_7_variables_30_01")
# Read data
if(F){
    ## CHECK THAT THIS IS THE CORRECT FILE
  fit_Net_7    <- as_cmdstan_fit(files = paste0("Datas/HOMEPC_VAR_7-202502050525-",1:8,".csv") ); gc()
  draws_data_7 <- as_draws_df(fit_Net_7$draws(variables = c(inference_vars_regex_alpha, "lp__") )); gc(); rm(fit_Net_7); gc()
}
# Diagnostics and posterior distribution marignal plots. 
plotnams <- inference_vars_regex_alpha;pdf(file = paste0("Datas/Bayespots_7VAR_",format(Sys.time(), "%Y-%m-%d"), ".pdf"));color_scheme_set("viridis")
for(i in plotnams){print(  mcmc_trace(draws_data_7, regex_pars = i ));print(  mcmc_areas_ridges(draws_data_7[ , grep(i, names(draws_data_7))]) );
}; gc(); dev.off()

# Compute RMSEA, distribution
estimated_var_samples_7 <- draws_data_7[,c( grep("A", names(draws_data_7)) , grep("Omega", names(draws_data_7)) )]; As <- grep("A_", names(estimated_var_samples_7)); Os <- grep("Omega", names(estimated_var_samples_7)); gc()
var_samples_7 <- pbapply::pblapply(1:nrow(draws_data_7),
                                 FUN = function(i){
                                   A_temp <- matrix( unlist(estimated_var_samples_7[i,As]), ncol = sqrt(length(As)), nrow = sqrt(length(Os)));
                                   Z_temp <- matrix( unlist(estimated_var_samples_7[i,Os]), ncol = sqrt(length(As)), nrow = sqrt(length(Os)));
                                   return(list(A = A_temp, Z = Z_temp))
                                 }); gc()
rmseas_7 <- pbapply::pbsapply(var_samples_7, FUN = function(x){
  res <- try(civ_find(x$A,x$Z, cov.difference = T, N = nrow(Data5b))$RMSEA)
  if(is.character(res)) {res <- NA; warning("Non-convergence")}
  return(res) }); gc()
saveRDS(rmseas_7, file = "Datas/rmseas_7vars.RDS");gc()

props_explained_7 <- pbapply::pbsapply(1:length(var_samples_7), FUN = function(x){
  return(try(Re(civ_parallel(var_samples_7[[x]]$A,var_samples_7[[x]]$Z)$prop_explained[1])))
}); gc()
saveRDS(props_explained_7, file = "Datas/props_explained_7vars.RDS"); gc()

source("Yuan 2016.R"); Yuan_2016(  (2*7*(2*7-1) / 2) - (1+7+14),4135) 

  # Congruency
congruency_7 <- pbapply::pbsapply(1:length(var_samples_7), FUN = function(x){
  return(try(Re(civ_parallel(var_samples_7[[x]]$A,var_samples_7[[x]]$Z)$min_factor_congruency)))
}); gc()
saveRDS(congruency_7, file = "Datas/congruency_7vars.RDS"); gc()


