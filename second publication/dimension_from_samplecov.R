# --- 0. Setup and Libraries ---
if(!require("ggplot2")) install.packages("ggplot2")
if(!require("dplyr")) install.packages("dplyr")
if(!require("tidyr")) install.packages("tidyr")
if(!require("cowplot")) install.packages("cowplot")

library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)

set.seed(999) 

# ==============================================================================
# 1. MODEL DEFINITION
# ==============================================================================

p <- 3
mu_base <- c(0, 0, 0)

# A. Base System
Theta_base <- matrix(c(
  0.5,  0.1, 0.0,
  -0.1, 0.8, 0.0,
  0.0,  0.0, 1.0   
), nrow=p, ncol=p, byrow=TRUE)

# Base Diffusion (Target Stationary)
diag_vals_target <- c(1.0, 0.5, 0.0) 

# B. Rotation
rnd_mat <- matrix(rnorm(9), 3, 3)
Rot <- qr.Q(qr(rnd_mat))

# C. Rotated Parameters
Theta <- Rot %*% Theta_base %*% t(Rot)
mu    <- Rot %*% mu_base

# D. Theoretical Stationary Covariance (for reference)
Sigma_base <- diag(diag_vals_target)
BBt <- (Rot %*% Sigma_base) %*% t(Rot %*% Sigma_base)
I_mat <- diag(p)
L_op <- kronecker(I_mat, Theta) + kronecker(Theta, I_mat)
vec_V <- solve(L_op, as.vector(BBt))
True_C0 <- matrix(vec_V, nrow=p, ncol=p)

true_eigs <- sort(eigen(True_C0)$values, decreasing = TRUE)
avg_signal_var <- mean(diag(True_C0))
noise_sd <- sqrt(avg_signal_var) 

cat("Target Stationary Eigenvalues:", round(true_eigs, 4), "\n")

# ==============================================================================
# 2. SIMULATION FUNCTIONS
# ==============================================================================

mat_exp <- function(A) {
  e <- eigen(A)
  V <- e$vectors
  D <- diag(exp(e$values))
  return(Re(V %*% D %*% solve(V)))
}

sigmoid <- function(t, t0, k) {
  1 / (1 + exp(-k * (t - t0)))
}

# --- NEW: Generate Full Trajectory First ---
generate_full_trajectory <- function(n_steps, dt, Theta, mu, Rot, mode="stationary") {
  
  Phi <- mat_exp(-Theta * dt)
  X <- matrix(0, nrow=n_steps, ncol=p)
  X[1, ] <- mu 
  curr_centered <- numeric(p)
  
  # Parameters for Abrupt Transition
  # We place the transition exactly in the middle of the Fixed Horizon
  t_mid <- n_steps / 2 
  k_slope <- 0.1 # Very abrupt/steep transition
  
  # Pre-calc Stationary Cholesky
  Sigma_stat <- Rot %*% diag(diag_vals_target)
  Q_stat <- (Sigma_stat %*% t(Sigma_stat)) * dt 
  chol_Q_stat <- t(chol(Q_stat))
  
  for (t in 1:(n_steps - 1)) {
    
    if (mode == "stationary") {
      chol_Q <- chol_Q_stat
    } else {
      if (mode == "decay") {
        # Max * (1 - sigmoid)
        val_2 <- 0.5 * (1 - sigmoid(t, t_mid, k_slope))
      } else if (mode == "growth") {
        # Max * sigmoid
        val_2 <- 0.5 * sigmoid(t, t_mid, k_slope)
      }
      
      current_diag <- c(1.0, val_2, 0.0)
      Sig_t <- Rot %*% diag(current_diag)
      Q_t <- (Sig_t %*% t(Sig_t)) * dt
      
      chol_Q <- tryCatch({
        suppressWarnings(t(chol(Q_t, pivot=TRUE)))
      }, error = function(e) {
        e_Q <- eigen(Q_t); val <- pmax(e_Q$values, 0)
        return(e_Q$vectors %*% diag(sqrt(val)) %*% t(e_Q$vectors))
      })
      
      if(is.matrix(chol_Q) && !is.null(attr(chol_Q, "pivot"))) {
        piv <- attr(chol_Q, "pivot")
        chol_Q <- chol_Q[order(piv), order(piv)]
      }
    }
    
    Z <- rnorm(p)
    curr_centered <- Phi %*% curr_centered + chol_Q %*% Z
    X[t+1, ] <- curr_centered + mu
  }
  return(X)
}

# --- NEW: Subsampling Runner ---
run_subsampling_experiment <- function(N_seq, reps, dt, Theta, mu, Rot, noise_sd, mode="stationary") {
  
  # Fixed Horizon: Long enough to contain the whole history
  TOTAL_STEPS <- 2000 
  
  results_list <- vector("list", length(N_seq) * reps)
  idx <- 1
  
  cat(paste0("Running Fixed-Window Experiment: ", mode, "...\n"))
  pb <- txtProgressBar(min = 0, max = reps, style = 3)
  
  for (r in 1:reps) {
    setTxtProgressBar(pb, r)
    
    # 1. Generate ONE full history for this replicate
    #    (The "Physical Reality" over the long period)
    X_full <- generate_full_trajectory(TOTAL_STEPS, dt, Theta, mu, Rot, mode)
    
    # Generate Noise Matrix (Full size) - only added if needed
    Noise_Full <- matrix(rnorm(TOTAL_STEPS * p, sd = noise_sd), nrow = TOTAL_STEPS, ncol = p)
    Y_full <- X_full + Noise_Full
    
    # 2. Subsample this history for different N
    for (n_sample in N_seq) {
      
      # systematic sampling: spread n_sample points evenly across 1:TOTAL_STEPS
      # We randomize the start index slightly to satisfy "Random Equidistant"
      step_size <- floor(TOTAL_STEPS / n_sample)
      start_idx <- sample(1:step_size, 1)
      indices <- seq(start_idx, by=step_size, length.out=n_sample)
      indices <- indices[indices <= TOTAL_STEPS] # Safety clip
      
      # Extract Data
      X_sub <- X_full[indices, , drop=FALSE]
      Y_sub <- Y_full[indices, , drop=FALSE]
      
      # Calculate Eigenvalues
      if (mode == "stationary") {
        # A (Clean) and B (Noisy)
        eigs_clean <- sort(eigen(cov(X_sub))$values, decreasing = TRUE)
        eigs_noisy <- sort(eigen(cov(Y_sub))$values, decreasing = TRUE)
      } else {
        # C and D: Clean ONLY (No Measurement Noise requested)
        eigs_clean <- sort(eigen(cov(X_sub))$values, decreasing = TRUE)
        eigs_noisy <- eigs_clean # Placeholder
      }
      
      results_list[[idx]] <- data.frame(
        Mode = mode,
        N = n_sample, Rep = r,
        Clean_Eig1 = eigs_clean[1], Clean_Eig2 = eigs_clean[2], Clean_Eig3 = eigs_clean[3],
        Noisy_Eig1 = eigs_noisy[1], Noisy_Eig2 = eigs_noisy[2], Noisy_Eig3 = eigs_noisy[3]
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

# We go up to N=1000, which is half the total resolution of the simulated reality
N_values <- seq(50, 1000, by = 50) 
Replicates <- 30 
dt <- 0.1

data_stat   <- run_subsampling_experiment(N_values, Replicates, dt, Theta, mu, Rot, noise_sd, mode="stationary")
data_decay  <- run_subsampling_experiment(N_values, Replicates, dt, Theta, mu, Rot, noise_sd, mode="decay")
data_growth <- run_subsampling_experiment(N_values, Replicates, dt, Theta, mu, Rot, noise_sd, mode="growth")

all_data <- bind_rows(data_stat, data_decay, data_growth)

get_summary <- function(data, prefix) {
  data %>%
    select(Mode, N, Rep, starts_with(prefix)) %>%
    pivot_longer(cols = starts_with(prefix), names_to = "Eigenvalue", values_to = "Value") %>%
    mutate(Eigenvalue = gsub(paste0(prefix, "_"), "", Eigenvalue)) %>%
    group_by(Mode, N, Eigenvalue) %>%
    summarise(
      Median = median(Value),
      Lower = quantile(Value, 0.10),
      Upper = quantile(Value, 0.90),
      .groups = "drop"
    )
}

clean_summ <- get_summary(all_data, "Clean")
noisy_summ <- get_summary(all_data, "Noisy")

plot_A_data <- filter(clean_summ, Mode == "stationary")
plot_B_data <- filter(noisy_summ, Mode == "stationary")
plot_C_data <- filter(clean_summ, Mode == "decay") 
plot_D_data <- filter(clean_summ, Mode == "growth")

# Formatting
true_lines <- data.frame(Eigenvalue = c("Eig1", "Eig2", "Eig3"), Value = true_eigs)
plot_cols <- hcl.colors(3, palette = "Dark 2") |> setNames(c("Eig1","Eig2","Eig3"))

create_plot <- function(data, title, subtitle, show_legend=FALSE) {
  p <- ggplot(data, aes(x=N, group=Eigenvalue)) +
    geom_ribbon(aes(ymin=Lower, ymax=Upper, fill=Eigenvalue), alpha=0.15) +
    geom_line(aes(y=Median, color=Eigenvalue), size=1) +
    geom_hline(data=true_lines, aes(yintercept=Value, color=Eigenvalue), linetype="dashed", size=0.8) +
    labs(title=title, subtitle=subtitle, y="Eigenvalue", x=if(show_legend) "Sample Size N (Subsampling)" else NULL) +
    scale_color_manual(values=plot_cols) + scale_fill_manual(values=plot_cols) +
    theme_minimal() + 
    theme(plot.title=element_text(face="bold", size=11),
          plot.subtitle=element_text(size=9))
  
  if(!show_legend) {
    p <- p + theme(legend.position="none")
  } else {
    p <- p + theme(legend.position="bottom")
  }
  return(p)
}

pA <- create_plot(plot_A_data, 
                  "A. Stationary (Clean)", 
                  "Fixed Window Subsampling: Converges correctly.")

pB <- create_plot(plot_B_data, 
                  "B. Stationary (Noisy)", 
                  "Fixed Window Subsampling: Constant noise bias.")

# Note for C and D: The "True" dashed line represents the INITIAL Stationary state.
# We expect the solid line to settle roughly HALFWAY between the dashed line and zero
# because the window covers both states equally.

pC <- create_plot(plot_C_data, 
                  "C. Sigmoidal Decay (Fixed Window)", 
                  "Result: Flat line at AVERAGE of 'Before' and 'After'. No drift, just consistently 'wrong'.", 
                  show_legend=FALSE)

pD <- create_plot(plot_D_data, 
                  "D. Sigmoidal Growth (Fixed Window)", 
                  "Result: Flat line at AVERAGE of 'Before' and 'After'. Scarcity doesn't hide the history.", 
                  show_legend=TRUE)

final_plot <- plot_grid(
  plot_grid(pA, pB, ncol=2),
  plot_grid(pC, pD, ncol=2),
  ncol=1, rel_heights=c(1, 1.2)
)

png(filename = "./Datas/Fixed_Window_Subsampling.png",
    width = 14, height = 10, units = "in", family = "sans", res = 150)
print(final_plot)
dev.off()