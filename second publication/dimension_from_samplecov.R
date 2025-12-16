# --- 0. Setup and Libraries ---
required_packages <- c("ggplot2", "dplyr", "cowplot", "tidyr")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)

set.seed(42)

# ==============================================================================
# 1. MODEL DEFINITION (Same as before)
# ==============================================================================

p <- 3
mu <- c(10, 0, -5)

# Mean Reversion Matrix (Theta)
Theta <- matrix(c(
  0.5,  0.1,  0.0,
  -0.1,  0.8,  0.1,
  0.0, -0.2,  1.0
), nrow=p, ncol=p, byrow=TRUE)

# Diffusion Matrix (Sigma)
Sigma <- matrix(c(
  1.0, 0.0, 0.0,
  0.5, 1.0, 0.0,
  0.2, 0.3, 1.0
), nrow=p, ncol=p, byrow=TRUE)

BBt <- Sigma %*% t(Sigma)

# Theoretical C(0) and Eigenvalues
I_mat <- diag(p)
L_op <- kronecker(I_mat, Theta) + kronecker(Theta, I_mat)
vec_V <- solve(L_op, as.vector(BBt))
True_C0 <- matrix(vec_V, nrow=p, ncol=p)
true_eigs <- sort(eigen(True_C0)$values, decreasing = TRUE)

# ==============================================================================
# 2. SIMULATION FUNCTION
# ==============================================================================

mat_exp <- function(A) {
  e <- eigen(A)
  V <- e$vectors
  D <- diag(exp(e$values))
  return(Re(V %*% D %*% solve(V)))
}

# Optimized for repeated runs
run_simulation_batch <- function(N_seq, reps, dt, Theta, mu, True_V) {
  
  # Pre-calculate transition parameters (fixed for all N if dt is fixed)
  Phi <- mat_exp(-Theta * dt)
  Q <- True_V - Phi %*% True_V %*% t(Phi)
  
  # Robust Cholesky
  chol_Q <- tryCatch({
    t(chol(Q, pivot = TRUE))
  }, error = function(e) {
    e_Q <- eigen(Q); val <- pmax(e_Q$values, 0)
    return(e_Q$vectors %*% diag(sqrt(val)) %*% t(e_Q$vectors))
  })
  if(is.matrix(chol_Q) && !is.null(attr(chol_Q, "pivot"))) {
    piv <- attr(chol_Q, "pivot")
    chol_Q <- chol_Q[order(piv), order(piv)]
  }
  
  # Storage list
  results_list <- vector("list", length(N_seq) * reps)
  idx <- 1
  
  cat("Running Simulations (Total batches:", length(N_seq), ")\n")
  pb <- txtProgressBar(min = 0, max = length(N_seq), style = 3)
  
  for (i in seq_along(N_seq)) {
    n <- N_seq[i]
    setTxtProgressBar(pb, i)
    
    for (r in 1:reps) {
      # 1. Generate Path
      # We can generate noise for all N at once, but matrix ops are fast enough
      X <- matrix(0, nrow=n, ncol=p)
      X[1, ] <- mu 
      Z <- matrix(rnorm(p * (n - 1)), nrow=p)
      
      curr_centered <- numeric(p)
      for (t in 1:(n - 1)) {
        curr_centered <- Phi %*% curr_centered + chol_Q %*% Z[, t]
        X[t+1, ] <- curr_centered + mu
      }
      
      # 2. Estimates
      Est_C0 <- cov(X)
      
      # Metric A: Error Norm
      diff_mat <- Est_C0 - True_V
      frob_error <- sqrt(sum(diff_mat^2))
      
      # Metric B: Eigenvalues
      est_eigs <- sort(eigen(Est_C0)$values, decreasing = TRUE)
      
      # Store
      results_list[[idx]] <- data.frame(
        N = n,
        Rep = r,
        Error_Norm = frob_error,
        Eig1 = est_eigs[1],
        Eig2 = est_eigs[2],
        Eig3 = est_eigs[3]
      )
      idx <- idx + 1
    }
  }
  close(pb)
  return(do.call(rbind, results_list))
}

# ==============================================================================
# 3. EXECUTE & PROCESS DATA
# ==============================================================================

# Parameters
N_values <- seq(20, 1000, by = 20) # N from 20 to 1000
Replicates <- 50                   # 50 Simulations per N

# Run
raw_data <- run_simulation_batch(N_values, Replicates, 0.1, Theta, mu, True_C0)

# --- Process Data for Plotting (Calculate Quantiles) ---

# 1. Error Plot Data
error_summary <- raw_data %>%
  group_by(N) %>%
  summarise(
    Median = median(Error_Norm),
    Lower = quantile(Error_Norm, 0.10),
    Upper = quantile(Error_Norm, 0.90)
  )

# 2. Eigenvalue Plot Data (Needs Pivot to Long format first)
eig_long <- raw_data %>%
  select(N, Rep, Eig1, Eig2, Eig3) %>%
  pivot_longer(cols = starts_with("Eig"), names_to = "Eigenvalue", values_to = "Value")

eig_summary <- eig_long %>%
  group_by(N, Eigenvalue) %>%
  summarise(
    Median = median(Value),
    Lower = quantile(Value, 0.10),
    Upper = quantile(Value, 0.90),
    .groups = "drop"
  )

# ==============================================================================
# 4. PLOTTING
# ==============================================================================

# --- Plot 1: Estimation Error (Ribbon) ---
p1 <- ggplot(error_summary, aes(x = N)) +
  # Shaded Confidence Band (10th-90th Percentile)
  geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "steelblue", alpha = 0.3) +
  # Median Line
  geom_line(aes(y = Median), color = "darkblue", size = 1) +
  labs(
    title = "A. Estimation Error Consistency",
    subtitle = "Median Error (Line) ± 10th-90th Percentile (Shaded)",
    x = "Sample Size N",
    y = "Error (Frobenius Norm)"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"))

# --- Plot 2: Eigenvalue Convergence (Ribbon) ---
# True Values Dataframe for dashed lines
true_lines <- data.frame(
  Eigenvalue = c("Eig1", "Eig2", "Eig3"),
  Value = true_eigs
)

p2 <- ggplot(eig_summary, aes(x = N, group = Eigenvalue)) +
  # Shaded Bands
  geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = Eigenvalue), alpha = 0.2) +
  # Median Lines
  geom_line(aes(y = Median, color = Eigenvalue), size = 1) +
  # True Value Dashed Lines
  geom_hline(data = true_lines, aes(yintercept = Value, color = Eigenvalue), 
             linetype = "dashed", size = 0.8) +
  labs(
    title = "B. Eigenvalue Spectrum",
    subtitle = "Shaded regions represent empirical variation (simulated)",
    x = "Sample Size N",
    y = "Eigenvalue Magnitude"
  ) +
  theme_minimal() +
  scale_fill_manual(values = c("red", "green4", "blue")) +
  scale_color_manual(values = c("red", "green4", "blue")) +
  theme(plot.title = element_text(face = "bold"), legend.position = "bottom")

# --- Combine ---
final_plot <- plot_grid(p1, p2, ncol = 2, align = "h")
tiff(filename = "./Datas/Convergence_of_Eigenvalues.tiff", 
     width = 16, height = 12, units = "in", family = "sans", res = 120)
print(final_plot)
dev.off()
