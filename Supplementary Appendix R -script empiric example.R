# Packages -----
# List of required packages
required_packages <- c(
  "Matrix", "fastmatrix", "expm", "rstan",
  "qgraph", "tidyverse", 
  "ggplot2", "rstantools", "bayesplot", "cmdstanr", "posterior",
  "viridisLite", "dplyr")

# Function to check and install missing packages
for (pkg in required_packages) {
  if(pkg == "cmdstanr" & !requireNamespace(pkg, quietly = T)) install.packages("cmdstanr", repos = c('https://stan-dev.r-universe.dev', getOption("repos"))) 
  else if(!requireNamespace(pkg, quietly = T)) {install.packages(pkg, dependencies = T)}
    library(pkg, character.only = TRUE)
}



# Load data ------

# load 
load(file.path("Fried_2022 data/clean_network.RData")); gc()
Data5b <- Data2

# Variables to investigate:
vars <- c("Q1","Q2","Q3","Q4","Q5","Q6","Q7","Q8","Q9",
          "Q10","Q11","Q12","Q13","Q14","Q15","Q16","Q17","Q18")

# Rename, and pick variables, which will be used for estimation.
varLabs <- c("Relax","Irritable","Worry","Nervous","Future","Anhedonia",
             "Tired","Hungry","Alone","Angry","Social_offline","Social_online","Music",
             "Procrastinate","Outdoors","C19_occupied","C19_worry","Home")
names(Data5b)[names(Data5b) %in% vars] <- varLabs

varLabs <- c("Relax", "Worry",  
             "Nervous")
varLabs2 <- c("Relax", "Worry",  
              "Nervous", "Tired", "Hungry",
              "Alone", "Angry")

# Exclude double beeps:
Data5b <- Data5b %>% group_by(id) %>%
  filter(beep != lag(beep)) %>% ungroup

# Exclude the IDs with lower than 50 obs count:
excluded_ids <- Data5b %>% group_by(id) %>% summarize(n = length(id)) %>% arrange(n, id) %>%
  filter(n < 50) %>% select(id) %>% unlist()
      # (In total, three IDS.)

# Exclude the IDs with more than 7 days worth (half) of missing values in any variable
excluded_ids <- base::union(excluded_ids,
                            Data5b %>%
                              group_by(id) %>%
                              summarise(across(all_of(varLabs), ~ sum(is.na(.)))) %>%
                              filter(if_any(all_of(varLabs), ~ . > 28)) %>%
                              pull(id))
      # (In total, four IDS.)

# Remove excluded ids from data.
Data5b <- Data5b %>% filter(!(id %in% excluded_ids))

# Truncate to 4, as asnwer category 5 has low cell counts.
Data5b <- Data5b %>% mutate(across(all_of(varLabs2), ~ pmin(.x, 4)))


# STAN Bayesian multilevel analysis. -----------------

  # Data preparations for Stan:
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
  # Missing mask is used to parameterize the missing values.
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

cutpoint_prior_locations <- c(-0.6, 1, 2.2)
cutpoint_count           <- length(cutpoint_prior_locations)

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
  # Note: Stan models took 12 - 24 hours to run on 
  # Intel(R) Core(TM) Ultra 7 165H, 1400 Mhz, 16 Core(s), 22 Logical Processor(s).

# Set of variables of interest, which we can set to pars argument, or later on extract.
inference_vars_regex <- c("A", "Omega", "cutpoints","time_of_day_intercept")
inference_vars_regex_alpha <- c("A_effective","A","psi","Lambda", "Omega","L_", "cutpoints", "time_of_day_effect",
                                "ref_time_of_day_effect", "specific_time_of_day_effect")
# Run MCMC with cmdstanr
nchains = 8
mod_Net <- cmdstan_model("BayesianOrderedVAR_alpha.stan")
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
fit_Net$save_output_files("/Datas/", basename = "3VAR_CF_relaxedpriors", timestamp = T, random = F)

# Read data
if(F){
    fit_Net    <- as_cmdstan_fit(files = paste0("Datas/3VAR_CF_relaxedpriors-202502031854-",1:nchains,".csv") ); gc()
  draws_data <- as_draws_df(fit_Net$draws(variables = c(inference_vars_regex_alpha) ), .nhcains = nchains ); gc(); rm(fit_Net); gc()
}

  # Analysis.
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
source("Supplementary Appendix R -script var_dcf_compare.R")

# Figure 4 in main text ----
  # Compute parallel analysis imitation and RMSEA, for the posterior mean.
result_parallel  <- var_dcf_compare(A, Z)

  # Compute credible intervals for eigenvalues, congruencies.
eigen_congurency <- pbapply::pblapply(var_samples, FUN = function(x){
              res  <- try(var_dcf_compare(x$A,x$Z))
            eigens <- t(abs(res$eigenvals))
      congruencies <- res$subsequent_pair_congruencies
    return(list(eigens = eigens, congruencies = congruencies)) }); gc()

    # How many explain above 0.90?
mean(unlist(lapply(eigen_congurency, function(x) any(0.90 > (x$eigens[ , 1] / (rowSums(x$eigens)))) )))
    # All.
    # How many congruencies above 0.90?
mean(unlist(lapply(eigen_congurency, function(x) any(0.90 > x$congruencies))))
    # All.

eigen_dat <- data.frame(
  Re(do.call(rbind,
             lapply(eigen_congurency, 
                    FUN = function(x) cbind( x$eigens, 1:11 ) ))))

upper <- as.matrix(eigen_dat %>% group_by(X4) %>% reframe( across(paste0( "X", 1:(length(eigen_dat)-1) ), ~ quantile(.x, c(.975))) ))
lower <- as.matrix(eigen_dat %>% group_by(X4) %>% reframe( across(paste0( "X", 1:(length(eigen_dat)-1) ), ~ quantile(.x, c(.025))) ))

tiff(filename = "Figure_3.tiff", 
     width    = 17, 
     height   = 19, 
     units    = "cm", 
     res      = 300,
     pointsize = 10)
par(mfrow     = c(2,2) )
par(mar       = c(4,4,2,0.5) )
colvec <- cividis(ncol(upper)-1) |> adjustcolor(alpha.f = 0.15)
matplot(t(result_parallel$eigenvals),
        type = "n",
        ylab = "Absolute value of eigenvalue", 
        xlab = expression(paste("Increment in time ", Delta, "T")),
        xaxt = "n",
        main = "Cross-covariance eigenvalues",
        font.main = 1); grid()
axis(side = 1, at = 1:(ncol(result_parallel$eigenvals)), labels = 0:(ncol(result_parallel$eigenvals) - 1))
for( i in 2:ncol(upper)) {
  polygon(x = c(upper[,1], rev(lower[,1])), y = c(upper[,i], rev(lower[,i])),
          col = colvec[i-1], border = F)
}
matplot(t(result_parallel$eigenvals), 
        type = "b",
        col  = cividis(ncol(upper)-1), add = T)


cong_dat <- data.frame(Re(do.call(rbind,lapply(eigen_congurency, FUN = function(x) cbind( x$congruencies, 1:10 ) ))))
upper_c  <- as.matrix(cong_dat %>% group_by(X2) %>% 
                        reframe( quantile(X1, 0.975) ))
lower_c  <- as.matrix(cong_dat %>% group_by(X2) %>% 
                        reframe( quantile(X1, 0.025) ))

matplot(result_parallel$subsequent_pair_congruencies, type = "n",
        ylim = c(0,1),
        ylab = "Congruency coefficient", 
        xlab = "Cross-covariance pair",
        xaxt = "n",
        main = "Largest eigenvector congruency",
        font.main = 1); grid()
polygon(x = c(upper_c[,1], rev(upper_c[,1])), y=c(upper_c[,2], rev(lower_c[,2])),
        col = colvec[1])
axis(1, labels = paste0("(", 0:10,", ", 1:11,")"),
     at = 1:11, cex.axis = 0.7 )
matplot(result_parallel$subsequent_pair_congruencies, type = "b",
        col = cividis(6), add = T )
qgraph( A, layout = "circle", 
        labels = varLabs, mar = c(2,2,7,2))
title("Coefficient matrix",
      font.main = 1,
      line     = -1)
qgraph( Z, layout = "circle", 
        labels = varLabs, mar = c(3,3,7,3))
title("Innovation covariance",
      font.main = 1,
      line     = -1)

dev.off();gc();par(mfrow = c(1,1) )

# Second analysis: More symptoms -----

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
nchains = 8
mod_Net_7 <- cmdstan_model("BayesianOrderedVAR_alpha.stan")
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


# Extract posterior means as the VAR parameters:
A_7 <- matrix( unlist(colMeans(draws_data_7[,grep("A", names(draws_data_7))])),     ncol = K, nrow = K);
Z_7 <- matrix( unlist(colMeans(draws_data_7[,grep("Omega", names(draws_data_7))])), ncol = K, nrow = K)

# Obtain complete VAR(1) model samples:
estimated_var_samples_7 <- draws_data_7[,c( grep("A", names(draws_data_7)) , grep("Omega", names(draws_data_7)) )]; As <- grep("A_", names(estimated_var_samples_7)); Os <- grep("Omega", names(estimated_var_samples_7)); 
var_samples_7           <- pbapply::pblapply(1:nrow(draws_data_7),
                                           FUN = function(i){
                                             A_temp <- matrix( unlist(estimated_var_samples_7[i,As]), ncol = sqrt(length(As)), nrow = sqrt(length(Os)));
                                             Z_temp <- matrix( unlist(estimated_var_samples_7[i,Os]), ncol = sqrt(length(As)), nrow = sqrt(length(Os)));
                                             return(list(A = A_temp, Z = Z_temp))
                                           }); gc()
# Source the methods.
source("Supplementary Appendix R -script var_dcf_compare.R")

# Figure 5 in main text ----
# Compute parallel analysis imitation and RMSEA, for the posterior mean.
result_parallel_7  <- var_dcf_compare(A_7, Z_7)

# Compute credible intervals for eigenvalues, congruencies.
eigen_congurency_7 <- pbapply::pblapply(var_samples_7, FUN = function(x){
  res  <- try(var_dcf_compare(x$A,x$Z))
  eigens <- t(abs(res$eigenvals))
  congruencies <- res$subsequent_pair_congruencies
  return(list(eigens = eigens, congruencies = congruencies)) }); gc()

# How many do not explain above 0.90?
mean(unlist(lapply(eigen_congurency_7, function(x) any(0.90 > (x$eigens[ , 1] / (rowSums(x$eigens)))) )))
# All.
# How many congruencies not above 0.90?
mean(unlist(lapply(eigen_congurency_7, function(x) any(0.90 > x$congruencies))))
# 0.133.

# Compute quantiles for eigenvalues.
eigen_distributions_7 <- pbapply::pblapply(var_samples_7, FUN = function(x){
  res  <- var_dcf_compare(x$A,x$Z)
  eigens <- t((res$eigenvals))
  return(list(eigens = eigens)) }); gc()

eigen_distributions_7 |> 
  lapply(FUN = function(x) cbind(Re(x$eigens), 0:10)) %>% 
  do.call(rbind, args = .) |> 
  as.data.frame() |> setNames(c(paste0("eigen",1:7),"Increment")) %>%
  tibble() %>%
  group_by(Increment) %>%
  summarise_all( .funs = function(x) quantile(x, c(0.05, 0.025, 0.01, 0.001)) ) |>
  print(n = 50)

eigen_dat_7 <- data.frame(Re(do.call(rbind,lapply(eigen_congurency_7, 
                                                FUN = function(x) cbind( x$eigens, 1:11 ) ))))

upper_7 <- as.matrix(eigen_dat_7 %>% group_by(X8) %>% reframe( across(paste0( "X", 1:(length(eigen_dat_7)-1) ), ~ quantile(.x, c(.975))) ))
lower_7 <- as.matrix(eigen_dat_7 %>% group_by(X8) %>% reframe( across(paste0( "X", 1:(length(eigen_dat_7)-1) ), ~ quantile(.x, c(.025))) ))

tiff(filename = "Figure_5.tiff", 
     width    = 17, 
     height   = 19, 
     units    = "cm", 
     res      = 300,
     pointsize = 10)
par(mfrow     = c(2,2) )
par(mar       = c(4,4,2,0.5) )
colvec <- cividis(ncol(upper_7)-1) |> adjustcolor(alpha.f = 0.15)
matplot(t(result_parallel_7$eigenvals),
        col  = colvec,
        type = "n",
        ylab = "Absolute value of eigenvalue", 
        xlab = expression(paste("Increment in time ", Delta, "T")),
        xaxt = "n",
        main = "Cross-covariance eigenvalues",
        font.main = 1); grid()
axis(side = 1, at = 1:(ncol(result_parallel_7$eigenvals)), labels = 0:(ncol(result_parallel_7$eigenvals) - 1))

matplot(t(result_parallel_7$eigenvals), 
        type = "b",
        col  = cividis(ncol(upper_7)-1), add = T)
for( i in 2:ncol(upper_7)) {
  polygon(x = c(upper_7[,1], rev(lower_7[,1])), y = c(upper_7[,i], rev(lower_7[,i])),
          col = colvec[i-1], 
          border = NA)
}

cong_dat_7 <- data.frame(Re(do.call(rbind,lapply(eigen_congurency_7, FUN = function(x) cbind( x$congruencies, 1:10 ) ))))
upper_c_7  <- as.matrix(cong_dat_7 %>% group_by(X2) %>% 
                        reframe( quantile(X1, 0.975) ))
lower_c_7  <- as.matrix(cong_dat_7 %>% group_by(X2) %>% 
                        reframe( quantile(X1, 0.025) ))
matplot(result_parallel_7$subsequent_pair_congruencies, type = "n",
        ylim = c(0,1),
        ylab = "Congruency coefficient", 
        xlab = "Cross-covariance pair",
        xaxt = "n",
        main = "Largest eigenvector congruency",
        font.main = 1); grid()
polygon(x = c(upper_c_7[,1], rev(upper_c_7[,1])), y=c(upper_c_7[,2], rev(lower_c_7[,2])),
        col = adjustcolor(cividis(1), alpha.f = 0.15), border = NA )
axis(1, labels = paste0("(", 0:10,", ", 1:11,")"),
     at = 1:11, cex.axis = 0.7 )
matplot(result_parallel_7$subsequent_pair_congruencies, type = "b",
        col = cividis(6), add = T )
qgraph( A_7, layout = "circle", 
        labels = varLabs2, mar = c(2,2,7,2) )
title("Coefficient matrix",
      font.main = 1,
      line     = -1)
qgraph( Z_7, layout = "circle", 
        labels = varLabs2, mar = c(3,3,7,3))
title("Innovation covariance",
      font.main = 1,
      line     = -1)
dev.off();gc();par(mfrow = c(1,1) )


# Save for future use.
if(F) {
  
  saveRDS(result_parallel_7, file = "parallel_7.RDS")
  saveRDS(eigen_congurency_7, file = "eigen_congurency_7.RDS")
  saveRDS(draws_data_7, file = "draws_data_7.RDS"); gc()
  
  result_parallel_7  <- readRDS("parallel_7.RDS")
  eigen_congurency_7 <- readRDS("eigen_congurency_7.RDS")
  draws_data_7       <- readRDS("draws_data_7.RDS")
  }