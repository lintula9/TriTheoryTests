# Simulation overview. ----
# VAR(1) signals X^(i) via MTS::VARMAsim + noise eps^(i) via MTS
# Observed Y^(i) = X^(i) + eps^(i)
# Fit VAR(1) to Y^(i), compare estimated (Ahat, Sigmahat) to true (A, Sigma_u)

# Install packages. ----
pkgs <- c("MTS", "ggplot2", "gridExtra", "rstan")
for (pkg in pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg)}  else {
      library(pkg, character.only = T)}}

# Reproducibility seed. ----
set.seed(42)  # global reproducibility

# Parameters. ----
n_obs <- 1000        # kept observations
p     <- 3
L     <- 10          # "multiple lags" = 10
burn  <- 200         # burn-in dropped by skip=
SNR   <- 1.0         # signal-to-noise ratio per dimension (1.0 ≈ comparable scales)

# True VAR(1) for X. ----
A        <- matrix(c(0.50, 0.10, 0.00,
                     0.00, 0.40, 0.20,
                     0.00, 0.00, 0.30), nrow = p, byrow = TRUE)
Sigma_u  <- diag(p)

# Helpers. ----
hstack <- function(mats) do.call(cbind, mats)

scale_noise_to_SNR <- function(E, X, snr = 1.0) {
  # Scales each column of E so Var(E_j) ≈ Var(X_j) / snr
  E2 <- E
  for (j in 1:ncol(E)) {
    vx <- stats::var(X[, j]); ve <- stats::var(E[, j])
    if (is.finite(vx) && is.finite(ve) && ve > 0) {
      E2[, j] <- E[, j] * sqrt(vx / (snr * ve))
    }
  }
  E2
}

gen_VAR1 <- function(n_obs, A, Sigma_u, burn, seed) {
  set.seed(seed)
  sim <- MTS::VARMAsim(nobs = n_obs, arlags = 1, phi = A, sigma = Sigma_u, skip = burn)
  list(X = sim$series, U = sim$noises)   # return series and innovations
}

# Simulate signals X^(i) (independent paths). ----
sim1 <- gen_VAR1(n_obs, A, Sigma_u, burn, seed = 1001)
sim2 <- gen_VAR1(n_obs, A, Sigma_u, burn, seed = 1002)
sim3 <- gen_VAR1(n_obs, A, Sigma_u, burn, seed = 1003)
sim4 <- gen_VAR1(n_obs, A, Sigma_u, burn, seed = 1004)

X_case1 <- sim1$X; U_case1 <- sim1$U
X_case2 <- sim2$X; U_case2 <- sim2$U
X_case3 <- sim3$X; U_case3 <- sim3$U
X_case4 <- sim4$X; U_case4 <- sim4$U

# Simulate eps^(i) with MTS. ----
# Case 1: eps is diagonal VAR(10) (multi-lag autocorr, independent across dims)
base_phi <- 0.4 * (0.8)^(0:(L - 1))         # stable AR(10) profile
Phi_blocks_case1 <- lapply(1:L, function(k) diag(rep(base_phi[k], p)))
Phi_case1 <- hstack(Phi_blocks_case1)       # p x (p*L)
simE1 <- VARMAsim(nobs = n_obs, arlags = 1:L, phi = Phi_case1, sigma = diag(p), skip = 0)
e1 <- simE1$series

# Case 2: eps is diagonal VAR(1) (single-lag autocorr, independent across dims)
phi_vec <- c(0.5, 0.3, 0.6)
Phi_case2 <- diag(phi_vec)
simE2 <- VARMAsim(nobs = n_obs, arlags = 1, phi = Phi_case2, sigma = diag(p), skip = 0)
e2 <- simE2$series

# Case 3: eps is white (VAR(1) with phi = 0 => eps_t = innovations)
Phi_zero <- matrix(0, nrow = p, ncol = p)
simE3 <- VARMAsim(nobs = n_obs, arlags = 1, phi = Phi_zero, sigma = diag(p), skip = 0)
e3 <- simE3$series

# Case 4: eps is VMA(10) with dense cross-dim Theta_k (unstructured)
set.seed(123)  # reproducible Theta_k
Theta_blocks_case4 <- lapply(1:L, function(k) {
  M <- matrix(runif(p * p, -0.15, 0.15), nrow = p)
  (0.85)^(k - 1) * M
})
Theta_case4 <- hstack(Theta_blocks_case4)  # p x (p*L)
simE4 <- VARMAsim(nobs = n_obs, malags = 1:L, theta = Theta_case4, sigma = diag(p), skip = 0)
e4 <- simE4$series

# Build observed data Y = X + scaled noise ----
e1s <- scale_noise_to_SNR(e1, X_case1, SNR)
e2s <- scale_noise_to_SNR(e2, X_case2, SNR)
e3s <- scale_noise_to_SNR(e3, X_case3, SNR)
e4s <- scale_noise_to_SNR(e4, X_case4, SNR)

Y_case1 <- X_case1 + e1s
Y_case2 <- X_case2 + e2s
Y_case3 <- X_case3 + e3s
Y_case4 <- X_case4 + e4s

# Define the data list.
stan_data <- list(
  # Y       = scale(Y_case1), # Diagonal VAR(10) noise.
  Y = scale(Y_case2), # Diagonal VAR(1) noise.
  # Y = scale(Y_case3), # White noise.
  # Y = scale(Y_case4), # VMA(1) noise.
  N       = dim(Y_case1)[1],
  p       = dim(Y_case1)[2]
)