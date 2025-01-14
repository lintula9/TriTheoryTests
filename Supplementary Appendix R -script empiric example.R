# Original code was created by Eiko Fried, Faidra Papanikolaou & Sascha Epskamp. 
# This version of code was (heavily) edited in 29.11.2024 for reanalysis of Fried(2022)
# Mental Health and Social Contact During the COVID-19 Pandemic: 
# An Ecological Momentary Assessment Study.
# For the original code and data, see https://osf.io/kj5zh

# --------------- 1. Loading packages & Data ------------------------------
# List of required packages
required_packages <- c(
  "Matrix", "fastmatrix", "BVAR", "brms", "expm", "qgraph", "tidyverse"
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

# Install and load packages
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



# Select a plausible set of variables, which could be explained by a one CF model:
Data5b <- as_tibble(Data5b %>% select("Relax", "Worry", 
                            "Nervous", 
                            "id", "beep", "day", "conc"))
varLabs <- c("Relax", "Worry",  #Adjustment for below to work. 03.12.2024. Sakari Lintula
            "Nervous")

# Lagged variables
Data5b <- Data5b %>% 
  group_by(id) %>%
  mutate(Relax_lag = lag(Relax),
         Worry_lag = lag(Worry),
         Nervous_lag = lag(Nervous)) %>%
  ungroup()


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



# This is for the below sections transfer later!!!!!!!

# Cross-covariance function
var_ccov <- function(A,Z,Delta) {
    
  return( (A %^% Delta) %*% matrix(Matrix::solve(diag(1,ncol=ncol(A)^2,nrow=nrow(A)^2) - Matrix::kronecker(A,A)) %*% 
    fastmatrix::vec(Z), ncol = ncol(A), nrow = nrow(A)))
  }

pb <- txtProgressBar(min = 0, max = nrow(draws), style = 3)  # Progress bar
mad_draws <- c()
mse_draws <- c()
rmse_draws <- c()
for( i in 1:nrow(draws)) {
  setTxtProgressBar(pb, i)
  
  #Perform computations
  A <- A_array[ ,, i]
  Z <- Z_array[,, i]
  results <- civ_find2(A, Z)
  A_ <- results$A_result
  Z_ <- results$Z_result
  
  # Compute some summaries of the differneces
  mad_draws[ i ] <- mean(abs(c(c(A_ - A),c(Z_ - Z))))
  mse_draws[i] <- mean((c(A_ - A)^2 + c(Z_ - Z)^2))
  rmse_draws[i] <- sqrt(mean((c(A_ - A)^2 + c(Z_ - Z)^2)))
  
  # Compute RMSE of the differences in predicted covariance, up to arbitrary T.
  sqrt(mean(sapply(0:5, function(Delta) ((var_ccov(A,Z,Delta = Delta) - var_ccov(A,Z,Delta = Delta))^2  ))))
  
  }
# Close progress bar
close(pb)


# --------------------- 3. Obtain the closest indistinguishable model ------


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










# Scratch:
# Solve for the covariance
library(expm); library(Matrix)
Sigma_VAR = matrix(solve(diag(1, ncol = ncol(A)^2, nrow = nrow(A)^2) - fastmatrix::kronecker.prod(A)) %*% fastmatrix::vec(Z),
                   ncol = ncol(A), nrow = nrow(A))
K_VAR = function( Delta ) {
  
  A %^% Delta %*% Sigma_VAR
  }

eigen(Sigma_VAR) # Suggestive of there being one large component. 
# Capture Lambda as the normalized eigenvector of the largest eigenvalue.

L = eigen(Sigma_VAR)$vectors[,1] / c(sqrt(t(eigen(Sigma_VAR)$vectors[,1]) %*% eigen(Sigma_VAR)$vectors[,1]))
eigen(A)$vectors[,1] # The first vector of A is nearly in the same direction.
eigen(A)$vectors[,1] - (-1)*L 

  # Project A onto Lambda
A_0 =  L %*% t(L) %*% A %*% L %*% t(L)

  # Project the difference onto the orthogonal complement of Lambda
B_tilde = (A - A_0) %*% (diag(1, nrow = nrow(A), ncol = ncol(A)) - L %*% t(L))

# Create the VAR(1), indistinguishable from one dimensional D-CF(1) model.
A_tilde = A_0 + B_tilde
Z_tilde = L %*% t(L) %*% Z %*% L %*% t(L)# Project onto Lambda, as it must be proportional to LL^T


par(mfrow=c(2,2))
qgraph(Z, layout = "circle", labels = expression(X[1], X[2], X[3], X[4]))
title("VAR(1) and indistinguishable VAR(1).",outer = T)

qgraph(Z_tilde, layout = "circle", labels = expression(X[1], X[2], X[3], X[4]))
qgraph(A_tilde, layout = "circle", labels = expression(X[1], X[2], X[3], X[4]))
qgraph(A, layout = "circle", labels = expression(X[1], X[2], X[3], X[4]))
par(mfrow=c(1,1))

Sigma_DCF = matrix(solve(diag(1, ncol = ncol(A)^2, nrow = nrow(A)^2) - fastmatrix::kronecker.prod(A_tilde)) %*% fastmatrix::vec(Z_tilde),
                   ncol = ncol(A), nrow = nrow(A))
K_VAR = function( Delta ) {
  
  A_tilde %^% Delta %*% Sigma_DCF
}

