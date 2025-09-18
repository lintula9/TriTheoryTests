# Install packages. ----
pkgs <- c("MTS", "ggplot2", "gridExtra", "rstan")
for (pkg in pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg)}  else {
      library(pkg, character.only = T)}}

# Reproducibility seed. ----
set.seed(42)  # global reproducibility

# Parameters.
n_obs <- 10000       # kept observations
L     <- 10          # "multiple lags" = 10
burn  <- 200         # burn-in dropped by skip=
SNR   <- 1.0         # signal-to-noise ratio per dimension (1.0 ≈ comparable scales)

# Template VAR(1) coefficient matrix.
A               <- matrix(0, nrow = 20, ncol = 20)
A[upper.tri(A)] <- runif(n = upper.tri(A) |> sum(), -0.99, 0.99)

# Generate observed values at different time frequencies.
frequencies <- c(1,2,5,10,20,40,80,160)
frq_names   <- paste0("frequency_", frequencies)
num_var     <- 3:20
num_var_names <- paste0("num_var_", num_var)
result      <- lapply(num_var,
               FUN = \(p) {
                 true_simu <- MTS::VARMAsim(nobs  = n_obs, arlags = 1, 
                                            phi   = A[1:p,1:p], 
                                            sigma = diag(p), 
                                            skip  = burn)
                 lapply(frequencies,
               FUN = \(frq) {
                 true_simu$series[seq(from=1,to=nrow(true_simu$series),by=frq),]
                 }) |> setNames(frq_names) 
                 }) |> setNames(num_var_names)
               
# Fit models.
models      <- lapply(
  result, FUN = \(p) {
    lapply(p, FUN = \(frq){
      try(MTS::VARMA(da = frq, p = 1, q = 0))
    }) |> setNames(frq_names) 
    }) |> setNames(num_var_names)

# Eigen decomp each model implied within time point covariance matrix.
# Yule-Walker:
yule_walker <- function(model) {
  A <- model$phi; S <- model$sigma; p <- nrow(A)
  return(
    matrix(solve(diag(p) - fastmatrix::kronecker.prod(A))
                        %*% fastmatrix::vec(S),
                 ncol = p, nrow = p))}
eigen_decomps <- lapply(
  models, FUN = \(p) {
    lapply(p, FUN = \(frq){
      yule_walker(frq) |> eigen(only.values = T)
  }) |> setNames(frq_names) }) |> setNames(num_var_names)

# Define the data list.
stan_data <- list(
  # Y       = scale(Y_case1), # Diagonal VAR(10) noise.
  Y = scale(Y_case2), # Diagonal VAR(1) noise.
  # Y = scale(Y_case3), # White noise.
  # Y = scale(Y_case4), # VMA(1) noise.
  N       = dim(Y_case1)[1],
  p       = dim(Y_case1)[2]
)