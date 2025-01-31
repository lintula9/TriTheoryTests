# Scratch or old

# ncp from RMSEA  
DF = 12

qchisq(p = 0.95, df = DF, ncp = par)













# Obtain the estimates:

# using brms
# Formula for VAR(1) using `brms`
formula <- bf(Relax ~ Relax_lag + Worry_lag + Nervous_lag + (1 | id)) +
  bf(Worry ~ Relax_lag + Worry_lag + Nervous_lag + (1 | id)) +
  bf(Nervous ~ Relax_lag + Worry_lag + Nervous_lag + (1 | id)) +
  set_rescor(TRUE)  # Enable estimation of the residual covariance matrix

# Fitting VAR(1) for the whole dataset (clustered data)
res_bayes_all <- brm(
  formula = formula,
  data = Data5b %>%
    na.omit(),
  family = gaussian(),
  chains = 4,
  iter = 2000,
  cores = 4,
  control = list(adapt_delta = 0.95)
)

if(F) {
  saveRDS(res_bayes_all, file = "BayesMLVAR.RDS")
  res_bayes_all <- readRDS("BayesMLVAR.RDS")
}

# Extract residual covariance matrix for individual or whole dataset
# Summary of the model to review residual correlations
summary(res_bayes_all)

# retrieve A
A <- matrix(summary(res_bayes_all)$fixed$Estimate[4:12], ncol = 3, nrow = 3)

# retrieve Z
Z_cor <- summary(res_bayes_all)$rescor_pars$Estimate
Z <- matrix(0,ncol = 3, nrow = 3)
diag(Z) <- 1
Z[lower.tri(Z, diag = F)] <- Z_cor
Z[upper.tri(Z, diag = F)] <- Z_cor

Z <- diag(summary(res_bayes_all)$spec_pars$Estimate) %*%
  Z %*% 
  diag(summary(res_bayes_all)$spec_pars$Estimate)


# Obtain the posterior draws for A matrices

draws <- as_draws_df(res_bayes_all)
A_draws <- draws %>% select(matches("^b_.*_lag$"))
A_array <- array(0, dim = c(nrow(A), ncol(A), nrow(A_draws) ))
for( i in 1:nrow(A_draws)) {
  A_array[ , , i ] <- matrix(as.numeric(A_draws[ i , ], ncol = 3, nrow = 3))
}

# Extract posterior draws for the Z matrix
# Convert posterior samples to a data frame

# Extract residual standard deviations and correlations using regex
residual_sd_draws <- draws %>% select(matches("^sigma"))
residual_cor_draws <- draws %>% select(matches("^rescor"))

# Compute posterior draws of the covariance matrix
Z_array <- array(0, dim = c(nrow(Z), ncol(Z), nrow(residual_sd_draws) ))
for( i in 1:nrow(residual_sd_draws)) {
  
  #Transform to covariance
  Z_cor_temp <- as.numeric(residual_cor_draws[ i , ])
  Z_temp <- matrix(0,ncol = 3, nrow = 3)
  diag(Z_temp) <- 1
  Z_temp[lower.tri(Z_temp, diag = F)] <- Z_cor_temp
  Z_temp[upper.tri(Z_temp, diag = F)] <- t(Z_temp)[upper.tri(Z_temp)]
  
  Z_temp <- diag(as.numeric(residual_sd_draws[  i , ])) %*%
    Z_temp %*% 
    diag(as.numeric(residual_sd_draws[  i , ]))
  
  #place.
  Z_array[ ,, i] <- Z_temp
  
}

## Inference -----

pb <- txtProgressBar(min = 0, max = nrow(draws), style = 3)  # Progress bar
mad_draws <- c()
mse_draws <- c()
rmse_draws <- c()
rmsea_stat <- c()
for( i in seq(1,nrow(draws),by = 10)) {
  setTxtProgressBar(pb, i)
  #Perform computations
  A <- A_array[ ,, i]
  Z <- Z_array[ ,, i]
  results <- civ_find3(A, Z)
  A_ <- results$A_result
  Z_ <- results$Z_result
  # Compute some summaries of the differences
  mad_draws[ i ] <- mean(abs(c(c(A_ - A),c(Z_ - Z))))
  mse_draws[i] <- mean((c(A_ - A)^2 + c(Z_ - Z)^2))
  rmse_draws[i] <- sqrt(mean((c(A_ - A)^2 + c(Z_ - Z)^2)))
  # Compute a statistic for RMSE(A), up to arbitrary T = 72 = maximum time points observed for a person.
  T_ = 72
  rmsea_stat[i] <- sum(sapply(  1:T_, 
                                function(k) sum( ( var_ccov(A,Z,Delta=k) - var_ccov(A_,Z_,Delta=k) )^2 ))) +
    sum( ( (var_ccov(A,Z,Delta=0) - var_ccov(A_,Z_,Delta=0))[tril(selector)] )^2 ) #FIX
}

# Close progress bar
close(pb)



if(F) {
  saveRDS(mad_draws, "mad_draws.RDS")
  saveRDS(mse_draws, "mse_draws.RDS")
  saveRDS(rmse_draws, "rmse_draws.RDS")
  saveRDS(rmsea_stat, "rmsea_stat.RDS")
  
  mad_draws <- readRDS("mad_draws.RDS")
  mse_draws <- readRDS("mse_draws.RDS")
  rmse_draws <- readRDS("rmse_draws.RDS")
  rmsea_stat <- readRDS("rmsea_stat.RDS")
  
}


# Compute RMSEA with 95%CI. This is heuristic, since the covariance in this case is arbitrarily large and
# as such there is no way to really compute df_...
N_ = nrow(Data5b)
df_ <- length(A) + sum(upper.tri(Z,diag = T)) - (length(varLabs) + 1) # factor loadings count plus CF autoregression.
RMSEA <- sqrt( ( rmsea_stat ) / ( df_ )  - 1/(N_-1)  )

psych::describe(cbind(mad_draws, mse_draws, RMSEA))

# MAD, MSE densities
ggplot() +
  geom_density(aes(x = mad_draws))
ggplot() +
  geom_density(aes(x = mse_draws))
ggplot() +
  geom_density(aes(x = RMSEA))




# --------------------- 3. Plotting ------

# Find the closest indistinguishable VAR. Note, that there are many.
#Method 1
dcf_var <- civ_find(A,Z, tol = 1e-10, n.iter = 3000)
dcf_var$A_result - A
dcf_var$Z_result - Z
dcf_var$Loadings
#Method 2
dcf_var2 <- civ_find2(A,Z, tol = 1e-10, n.iter = 3000)
dcf_var2$A_result - A
dcf_var2$Z_result - Z
dcf_var2$Loadings

# Plot the result

# Plot.
max_weight_Z <- max(c(dcf_var2$Z_result,Z))
max_weight_A <- max(c(A,dcf_var2$A_result))


par(mfrow = c(2, 2), oma = c(0, 0, 4, 0)) # Adjust oma for outer margin to accommodate the title
labels = varLabs
qgraph(Z, 
       title = "Contemporaneous covariance\nEstimated VAR",
       title.cex = 1.5,
       mar = c(4, 4, 6, 4),
       layout = "circle", 
       edge.width = 2,
       maximum = max_weight_Z,
       labels = labels # Use the expression labels
)
qgraph(A, 
       title = "Lagged effects\n",
       title.cex = 1.5,  
       mar = c(4, 4, 6, 4),
       maximum = max_weight_A, # Consistent scale for edges
       edge.width = 3,       # Adjust this value for larger edges
       layout = "circle",
       labels = labels       # Use the expression labels
)
qgraph(dcf_var2$Z_result, 
       title = "\nIndistinguishable VAR",
       title.cex = 1.5,  
       mar = c(4, 4, 6, 4),
       maximum = max_weight_Z, # Consistent scale for edges
       edge.width = 2,       # Adjust this value for larger edges
       layout = "circle",
       labels = labels       # Use the expression labels
)

qgraph(dcf_var2$A_result, 
       title = "",
       title.cex = 1.5,  
       mar = c(4, 4, 6, 4),
       maximum = max_weight_A, # Consistent scale for edges
       edge.width = 3,       # Adjust this value for larger edges
       layout = "circle",
       labels = labels       # Use the expression labels
)

mtext("Estimated and an indistinguishabel VAR(1) Network models.", outer = TRUE, cex = 1.5, font = 2)

par(mfrow = c(1,1))


# Simplest method of inspecting eigenvectors of A and Z:

which.min( eigen(Z)$vectors - eigen(A)$vectors ) 





# 4. Second analysis, repeat of first with more variables. --------------------


# Obtain the estimates:
# Formula for VAR(1) using `brms`
formulas <- lapply(varLabs2, function(var) {
  bf(as.formula(
    paste0(var, " ~ ", 
           paste(paste0(varLabs2, "_lag"), collapse = " + "), 
           " + (1 | id)"))
  )
})
# Combine the individual formulas and set residual correlation
formula <- Reduce(`+`, formulas) + set_rescor(TRUE)

# Fitting VAR(1) for the whole dataset (clustered data)
res_bayes_second <- brm(
  formula = formula,
  data = Data5b %>%
    na.omit(),
  family = gaussian(),
  chains = 4,
  iter = 2000,
  cores = 4,
  control = list(adapt_delta = 0.95)
)

if(F) {
  saveRDS(res_bayes_second, file = "BayesMLVAR_second.RDS")
  res_bayes_second <- readRDS("BayesMLVAR_second.RDS")
}

# Extract residual covariance matrix for individual or whole dataset
# Summary of the model to review residual correlations
summary(res_bayes_second)

# retrieve A
estimates <- summary(res_bayes_second)$fixed$Estimate
estimates[ -c(1:length(varLabs2)) ]
A <- matrix(estimates[ -c(1:length(varLabs2)) ], 
            ncol = length(varLabs2), nrow = length(varLabs2))

# retrieve Z
Z_cor <- summary(res_bayes_second)$rescor_pars$Estimate
Z <- matrix(0,ncol = nrow(A), nrow = ncol(A))
diag(Z) <- 1
Z[lower.tri(Z, diag = F)] <- Z_cor
Z[upper.tri(Z, diag = F)] <- Z_cor

Z <- diag(summary(res_bayes_second)$spec_pars$Estimate) %*%
  Z %*% 
  diag(summary(res_bayes_second)$spec_pars$Estimate)


# Obtain the posterior draws for A matrices

draws_second <- as_draws_df(res_bayes_second)
A_draws <- draws_second %>% select(matches("^b_.*_lag$"))
A_array <- array(0, dim = c(nrow(A), ncol(A), nrow(A_draws) ))
for( i in 1:nrow(A_draws)) {
  A_array[ , , i ] <- matrix(as.numeric(A_draws[ i , ], ncol = 3, nrow = 3))
}

# Extract posterior draws for the Z matrix
# Convert posterior samples to a data frame

# Extract residual standard deviations and correlations using regex
residual_sd_draws_second <- draws_second %>% select(matches("^sigma"))
residual_cor_draws_second <- draws_second %>% select(matches("^rescor"))

# Compute posterior draws of the covariance matrix
Z_array <- array(0, dim = c(nrow(Z), ncol(Z), nrow(residual_sd_draws_second) ))
for( i in 1:nrow(residual_sd_draws_second)) {
  
  #Transform to covariance
  Z_cor_temp <- as.numeric(residual_cor_draws_second[ i , ])
  Z_temp <- matrix(0,ncol = ncol(A), nrow = nrow(A))
  diag(Z_temp) <- 1
  Z_temp[lower.tri(Z_temp, diag = F)] <- Z_cor_temp
  Z_temp[upper.tri(Z_temp, diag = F)] <- t(Z_temp)[upper.tri(Z_temp)]
  
  Z_temp <- diag(as.numeric(residual_sd_draws_second[  i , ])) %*%
    Z_temp %*% 
    diag(as.numeric(residual_sd_draws_second[  i , ]))
  
  #place.
  Z_array[ ,, i] <- Z_temp
  
}

## Inference. -----

# Cross-covariance function
var_ccov <- function(A,Z,Delta) {
  
  return( (A %^% Delta) %*% matrix(Matrix::solve(diag(1,ncol=ncol(A)^2,nrow=nrow(A)^2) - Matrix::kronecker(A,A)) %*% 
                                     fastmatrix::vec(Z), ncol = ncol(A), nrow = nrow(A)))
}

pb <- txtProgressBar(min = 0, max = nrow(draws_second), style = 3)  # Progress bar
mad_draws_second <- c()
mse_draws_second <- c()
rmse_draws_second <- c()
rmsea_stat_second <- c()
N_ <- nrow(Data5b %>%
             na.omit())
for( i in seq(1,N_,by = 10)) {
  setTxtProgressBar(pb, i)
  #Perform computations
  A <- A_array[ ,, i]
  Z <- Z_array[ ,, i]
  results <- civ_find2(A, Z)
  A_ <- results$A_result
  Z_ <- results$Z_result
  # Compute some summaries of the differences
  mad_draws_second[ i ] <- mean(abs(c(c(A_ - A), c(Z_ - Z))))
  mse_draws_second[i] <- mean((c(A_ - A)^2 + c(Z_ - Z)^2))
  rmse_draws_second[i] <- sqrt(mean((c(A_ - A)^2 + c(Z_ - Z)^2)))
  # Compute a statistic for RMSE(A), up to arbitrary T = 72 = maximum of observation time length in data.
  T_ = 72
  rmsea_stat_second[i] <- sum(sapply(  1:T_, # This is actually simply the fit, F, statistic (for unweighted least squares).
                                       function(k) sum( ( var_ccov(A,Z,Delta=k) - var_ccov(A_,Z_,Delta=k) )^2 ) )) +
    sum( (fastmatrix::vech( var_ccov(A,Z,Delta=0) - var_ccov(A_,Z_,Delta=0) ))^2 ) # Only use the non-redundant elements of the within time point covariance.
}

if(F) {
  saveRDS(mad_draws_second, "mad_draws_second.RDS")
  saveRDS(mse_draws_second, "mse_draws_second.RDS")
  saveRDS(rmse_draws_second, "rmse_draws_second.RDS")
  saveRDS(rmsea_stat_second, "rmsea_stat_second.RDS")
  
  mad_draws_second <- readRDS("mad_draws_second.RDS")
  mse_draws_second <- readRDS("mse_draws_second.RDS")
  rmse_draws_second <- readRDS("rmse_draws_second.RDS")
  rmsea_stat_second <- readRDS("rmsea_stat_second.RDS")
  
}

# Close progress bar
close(pb)

# Compute RMSEA with 95%CI. 
df_ <- 49 + 28 - 7 - 1
df_ <- length(A) + sum(upper.tri(Z, diag = T)) - (length(varLabs2) + 1) # factor loadings count plus CF autoregression.
RMSEA_second <- sqrt(rmsea_stat_second/df_)
median(RMSEA_second, na.rm = T)

psych::describe(cbind(mad_draws_second, 
                      mse_draws_second, 
                      RMSEA_second))

# MAD, MSE densities
ggplot() +
  geom_density(aes(x = mad_draws_second))
ggplot() +
  geom_density(aes(x = mse_draws_second))
ggplot() +
  geom_density(aes(x = RMSEA_second))
