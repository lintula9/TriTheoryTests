# --- 0. Setup and Libraries ---
required_packages <- c("ggplot2", "dplyr", "cowplot", "tidyr")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)

set.seed(999) # New seed for rotation

# ==============================================================================
# 1. MODEL DEFINITION: ROTATED RANK-2 SYSTEM
# ==============================================================================

p <- 3
mu_base <- c(0, 0, 0)

# A. Base System (2 Active, 1 Silent)
# We decouple the 3rd dimension and give it 0 noise.
Theta_base <- matrix(c(
  0.5,  0.1, 0.0,
  -0.1, 0.8, 0.0,
  0.0,  0.0, 1.0   # 3rd dim decays independently
), nrow=p, ncol=p, byrow=TRUE)

Sigma_base <- matrix(c(
  1.0, 0.0, 0.0,
  0.5, 1.0, 0.0,
  0.0, 0.0, 0.0   # 3rd dim has NO noise
), nrow=p, ncol=p, byrow=TRUE)

# B. Create Random Rotation Matrix
# QR decomposition of random matrix gives a valid orthogonal rotation matrix
rnd_mat <- matrix(rnorm(9), 3, 3)
Rot <- qr.Q(qr(rnd_mat))

# C. Rotate Parameters into 3D Space
# X_new = R * X_old
Theta <- Rot %*% Theta_base %*% t(Rot)
Sigma <- Rot %*% Sigma_base
mu    <- Rot %*% mu_base

BBt <- Sigma %*% t(Sigma) # This is now Rank 2

# D. Calculate Theoretical C(0)
I_mat <- diag(p)
L_op <- kronecker(I_mat, Theta) + kronecker(Theta, I_mat)
vec_V <- solve(L_op, as.vector(BBt))
True_C0 <- matrix(vec_V, nrow=p, ncol=p)

# Check Eigenvalues (Expect 3rd one to be 0)
true_eigs <- sort(eigen(True_C0)$values, decreasing = TRUE)

cat("Theoretical True Eigenvalues (Expect 3rd is approx 0):\n")
print(round(true_eigs, 6))

# --- Define Measurement Noise (50/50 SNR against the NON-ZERO signals) ---
# Average variance of the signal part only
avg_signal_var <- mean(diag(True_C0))
noise_var <- avg_signal_var 
noise_sd <- sqrt(noise_var)

cat("Noise Variance Added:", round(noise_var, 3), "\n")

# ==============================================================================
# 2. SIMULATION FUNCTION
# ==============================================================================

mat_exp <- function(A) {
  e <- eigen(A)
  V <- e$vectors
  D <- diag(exp(e$values))
  return(Re(V %*% D %*% solve(V)))
}

run_simulation_batch <- function(N_seq, reps, dt, Theta, mu, True_V, noise_sd) {
  
  Phi <- mat_exp(-Theta * dt)
  Q <- True_V - Phi %*% True_V %*% t(Phi)
  
  # Robust Cholesky (Handles Rank Deficient Q gracefully)
  chol_Q <- tryCatch({
    # pivoting is crucial here as Q is rank 2!
    suppressWarnings(t(chol(Q, pivot = TRUE))) 
  }, error = function(e) {
    e_Q <- eigen(Q); val <- pmax(e_Q$values, 0)
    return(e_Q$vectors %*% diag(sqrt(val)) %*% t(e_Q$vectors))
  })
  
  if(is.matrix(chol_Q) && !is.null(attr(chol_Q, "pivot"))) {
    piv <- attr(chol_Q, "pivot")
    chol_Q <- chol_Q[order(piv), order(piv)]
  }
  
  results_list <- vector("list", length(N_seq) * reps)
  idx <- 1
  
  cat("Running Simulations...\n")
  pb <- txtProgressBar(min = 0, max = length(N_seq), style = 3)
  
  for (i in seq_along(N_seq)) {
    n <- N_seq[i]
    setTxtProgressBar(pb, i)
    
    for (r in 1:reps) {
      # 1. Generate Signal X
      X <- matrix(0, nrow=n, ncol=p)
      X[1, ] <- mu 
      Z <- matrix(rnorm(p * (n - 1)), nrow=p)
      curr_centered <- numeric(p)
      for (t in 1:(n - 1)) {
        curr_centered <- Phi %*% curr_centered + chol_Q %*% Z[, t]
        X[t+1, ] <- curr_centered + mu
      }
      
      # 2. Add Measurement Noise (Full Rank Noise!)
      Noise_Mat <- matrix(rnorm(n * p, sd = noise_sd), nrow = n, ncol = p)
      Y <- X + Noise_Mat
      
      # --- Metric Calculation ---
      
      # A. Clean Data (Reference)
      Est_C0 <- cov(X)
      diff_mat <- Est_C0 - True_V
      frob_error <- sqrt(sum(diff_mat^2))
      clean_eigs <- sort(eigen(Est_C0)$values, decreasing = TRUE)
      
      # B. Noisy Data (Standard Covariance)
      Est_C0_Noisy <- cov(Y)
      noisy_eigs <- sort(eigen(Est_C0_Noisy)$values, decreasing = TRUE)
      
      # C. Dynamic Factor (Lag-1 Estimator)
      # Uses Autocovariance at lag 1 to filter out noise
      Gam1 <- cov(Y[1:(n-1), ], Y[2:n, ])
      Est_C0_Dynamic <- (Gam1 + t(Gam1)) / 2
      # Sort eigenvalues magnitude (sometimes small negative values appear due to finite sample)
      dfm_vals <- eigen(Est_C0_Dynamic)$values
      dfm_eigs <- sort(dfm_vals, decreasing = TRUE)
      
      results_list[[idx]] <- data.frame(
        N = n, Rep = r, Error_Norm = frob_error,
        Clean_Eig1 = clean_eigs[1], Clean_Eig2 = clean_eigs[2], Clean_Eig3 = clean_eigs[3],
        Noisy_Eig1 = noisy_eigs[1], Noisy_Eig2 = noisy_eigs[2], Noisy_Eig3 = noisy_eigs[3],
        DFM_Eig1   = dfm_eigs[1],   DFM_Eig2   = dfm_eigs[2],   DFM_Eig3   = dfm_eigs[3]
      )
      idx <- idx + 1
    }
  }
  close(pb)
  return(do.call(rbind, results_list))
}

# ==============================================================================
# 3. EXECUTE & PLOT
# ==============================================================================

N_values <- seq(20, 1000, by = 20)
Replicates <- 50
raw_data <- run_simulation_batch(N_values, Replicates, 0.1, Theta, mu, True_C0, noise_sd)

# Process Summaries
get_summary <- function(data, prefix) {
  data %>%
    select(N, Rep, starts_with(prefix)) %>%
    pivot_longer(cols = starts_with(prefix), names_to = "Eigenvalue", values_to = "Value") %>%
    mutate(Eigenvalue = gsub(paste0(prefix, "_"), "", Eigenvalue)) %>%
    group_by(N, Eigenvalue) %>%
    summarise(
      Median = median(Value),
      Lower = quantile(Value, 0.10),
      Upper = quantile(Value, 0.90),
      .groups = "drop"
    )
}

err_summary <- raw_data %>% group_by(N) %>% 
  summarise(Median=median(Error_Norm), Lower=quantile(Error_Norm,0.1), Upper=quantile(Error_Norm,0.9))
clean_summary <- get_summary(raw_data, "Clean")
noisy_summary <- get_summary(raw_data, "Noisy")
dfm_summary   <- get_summary(raw_data, "DFM")

# Plots
true_lines <- data.frame(Eigenvalue = c("Eig1", "Eig2", "Eig3"), Value = true_eigs)
plot_cols <- hcl.colors(3) |> setNames(c("Eig1","Eig2","Eig3"))

# A. Error
p1 <- ggplot(err_summary, aes(x = N)) +
  geom_ribbon(aes(ymin=Lower, ymax=Upper), fill="steelblue", alpha=0.3) +
  geom_line(aes(y=Median), color="darkblue", size=1) +
  labs(title="A. Clean Estimation Error", y="Frobenius Norm", x=NULL) +
  theme_minimal() + theme(plot.title=element_text(face="bold"))

# B. Clean (Reference)
p2 <- ggplot(clean_summary, aes(x=N, group=Eigenvalue)) +
  geom_ribbon(aes(ymin=Lower, ymax=Upper, fill=Eigenvalue), alpha=0.2) +
  geom_line(aes(y=Median, color=Eigenvalue), size=1) +
  geom_hline(data=true_lines, aes(yintercept=Value, color=Eigenvalue), linetype="dashed", size=0.8) +
  labs(title="B. Clean Eigenvalues (True Rank 2)", y="Magnitude", x=NULL) +
  scale_color_manual(values=plot_cols) + scale_fill_manual(values=plot_cols) +
  theme_minimal() + theme(plot.title=element_text(face="bold"), legend.position="none")

# C. Noisy (Standard Cov)
p3 <- ggplot(noisy_summary, aes(x=N, group=Eigenvalue)) +
  geom_ribbon(aes(ymin=Lower, ymax=Upper, fill=Eigenvalue), alpha=0.2) +
  geom_line(aes(y=Median, color=Eigenvalue), size=1) +
  geom_hline(data=true_lines, aes(yintercept=Value, color=Eigenvalue), linetype="dashed", size=0.8) +
  labs(title="C. Noisy Estimator (Standard)", subtitle="3rd Eigenvalue artificially inflated", y="Magnitude", x="Sample Size N") +
  scale_color_manual(values=plot_cols) + scale_fill_manual(values=plot_cols) +
  theme_minimal() + theme(plot.title=element_text(face="bold"), legend.position="none")

# D. Dynamic Factor (Lag-1)
p4 <- ggplot(dfm_summary, aes(x=N, group=Eigenvalue)) +
  geom_ribbon(aes(ymin=Lower, ymax=Upper, fill=Eigenvalue), alpha=0.2) +
  geom_line(aes(y=Median, color=Eigenvalue), size=1) +
  geom_hline(data=true_lines, aes(yintercept=Value, color=Eigenvalue), linetype="dashed", size=0.8) +
  labs(title="D. Dynamic Model (Noise Removed)", subtitle="Correctly recovers 0 for 3rd Eigenvalue", y="Magnitude", x="Sample Size N") +
  scale_color_manual(values=plot_cols) + scale_fill_manual(values=plot_cols) +
  theme_minimal() + theme(plot.title=element_text(face="bold"), legend.position="bottom")

final_plot <- plot_grid(
  plot_grid(p1, p2, ncol=2, align="h"),
  plot_grid(p3, p4, ncol=2, align="h"),
  ncol=1, rel_heights=c(1, 1.2)
)

png(filename = "./Datas/Convergence_of_Eigenvalues.png",
    width = 16, height = 12, units = "in", family = "sans", res = 120)
print(final_plot)
dev.off()
