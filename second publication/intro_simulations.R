# Install & load -----------------------------------------------------------
pkgs <- c("MTS","ggplot2","gridExtra","rstan","progress","tidyr","tibble",
          "dplyr","foreach","doParallel") 
for (pkg in pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}
# Parameters. ----
n_obs  <- 10000
burn   <- 200
simu_n <- 1000
p_max  <- 20
sample_n <- round(seq(2, 5000, length.out = 1000)) |> unique()
idx      <- sapply(sample_n, 
                   FUN = \(n) unique(round(seq(1, n_obs, length.out = n))))
ranks    <- expand.grid(3:p_max, sapply(3:p_max, \(x) max(3:x)), sample_n) |>
  setNames(c("obs_dimension","sample_n"))
# Parallel backend. ----
cl <- parallel::makeCluster(max(1, parallel::detectCores() - 1))
doParallel::registerDoParallel(cl)
# Compute (parallel over p). ----
ranks_sim <- foreach::foreach(
  sim = 1:simu_n,
  .combine  = "rbind"
) %dopar% {
  for (p in 2:p_max) {
    # Simulate VAR(1) for this p
    A <- matrix(0, p, p)
    A[upper.tri(A, diag = FALSE)] <- rnorm((p * (p - 1) / 2))
    diag(A)               <- runif(p, 0.05, 0.95) * (-1L)*(round(runif(p)))
    true_simu <- MTS::VARMAsim(nobs = n_obs, arlags = 1, phi = A, 
                               sigma = diag(p), skip = burn)
    acc <- vector("list", length = length(sample_n))
    k <- 0L
    for (j in seq_along(sample_n)) {
      obs_series <- true_simu$series[idx[[j]],  1:p, drop = FALSE]
      est_dim    <- sum(eigen(cov(obs_series))$values > .Machine$double.eps)
      k          <- k + 1L
      acc[[k]]   <- tibble::tibble(
        obs_dimension  = p,
        sample_n       = sample_n[j],
        est_dimension  = est_dim,
        simulation     = sim
      ) }}
  dplyr::bind_rows(acc)}
parallel::stopCluster(cl)
# Save the result. ----
if(F) saveRDS(ranks_sim, file = "./Datas/intro_simulation.rds")