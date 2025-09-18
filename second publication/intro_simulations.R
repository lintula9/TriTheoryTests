# Install packages. ----
pkgs <- c("MTS", "ggplot2", "gridExtra", "rstan")
for (pkg in pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg)}  else {
      library(pkg, character.only = T)}}
# Parameters.
n_obs <- 10000       # kept observations
burn  <- 200         # burn-in dropped by skip=
# Template VAR(1) coefficient matrix.
A               <- matrix(0, nrow = 20, ncol = 20)
A[upper.tri(A)] <- runif(n = upper.tri(A) |> sum(), -0.99, 0.99)
# Generate observed values at different time frequencies.
frequencies <- c(1,2,5,10,20,40,80,160)
frq_names   <- paste0("frequency_", frequencies)
num_var     <- 3:20
num_var_names <- paste0("num_var_", num_var)
result      <- lapply(
  num_var, FUN = \(p) {
    true_simu <- MTS::VARMAsim(
      nobs  = n_obs, arlags = 1, phi   = A[1:p,1:p], sigma = diag(p), skip=burn)
    lapply(frequencies,FUN = \(frq) {
      true_simu$series[seq(from=1,to=nrow(true_simu$series),by=frq),]
                 }) |> setNames(frq_names) }) |> setNames(num_var_names)
# Fit models.
models      <- lapply(
  result, FUN = \(p) {
    lapply(p, FUN = \(frq){
      try({MTS::VARMA(da = frq, p = 1, q = 0); invisible(gc())})
    }) |> setNames(frq_names) 
    }) |> setNames(num_var_names)
saveRDS(models, file = "./second publication/")
# Eigen decomp each model implied within time point covariance matrix.
# Yule-Walker:
yule_walker <- function(model) {
  if(class(model) == "try-error") {
    message("Try error found, returning 0 matrix.")
    return(matrix(0))} 
  else {
  A <- model$Phi; S <- model$Sigma; p <- nrow(A)
  return(matrix(solve(diag(p^2) - fastmatrix::kronecker.prod(A))
                        %*% fastmatrix::vec(S),
                 ncol = p, nrow = p))}}
eigen_decomps <- lapply(
  models, FUN = \(p) {
    lapply(p, FUN = \(frq){
      e <- (yule_walker(frq) |> eigen(only.values = T))$values
      # return(sum(round(e, digits = 8) != 0))
  }) |> setNames(frq_names) }) |> setNames(num_var_names)
# Create table out of freq vs. number of variables.
rank_df <- data.frame(eigen_decomps)
saveRDS(rank_df, file = "./second publication/model_dimensions.rds")