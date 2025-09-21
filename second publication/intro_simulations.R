# Simulation data of estimated vs. true dimension.
source("Libraries.R")
# Parameters. ----
n_obs  <- 10000
burn   <- 200
simu_n <- 100
p_max  <- 40
sample_n <- floor(seq(2, 5000, length.out = 1000)) |> unique()
idx      <- lapply(sample_n, 
                   FUN = \(n) unique(floor(seq(1, n_obs, length.out = n))))
# Parallel backend. ----
cl <- parallel::makeCluster(max(1, parallel::detectCores() - 1))
doSNOW::registerDoSNOW(cl)
# Compute (parallel over sim). ----
# Progress bar.
pb <- progress_bar$new(
  format = "sim = :simulation [:bar] :elapsed | eta: :eta",
  total = simu_n, width = 60)
progress <- function(n){pb$tick(tokens = list(simulation = n))} 

ranks_sim <- foreach::foreach(
  sim = 1:simu_n,
  .combine      = "rbind",
  .options.snow = list(progress = progress)
) %dopar% {
  k <- 0L
  acc <- vector("list", length = (p_max - 1L) * length(sample_n))
  for (p in 2:p_max) {
    # Simulate VAR(1) for this p
    A <- matrix(0, p, p)
    A[upper.tri(A, diag = FALSE)] <- rnorm((p * (p - 1) / 2))
    diag(A)               <- runif( #Random sign.
      p, 0.05, 0.95) * (-1L)*(sample(c(1,-1),p,replace=T))
    true_simu <- MTS::VARMAsim(nobs = n_obs, arlags = 1, phi = A, 
                               sigma = diag(p), skip = burn)
    for (j in seq_along(sample_n)) {
      obs_series <- true_simu$series[idx[[j]],  1:p, drop = FALSE]
      est_dim    <- sum(
        eigen(cov(obs_series),only.values = T)$values > 1e-8)
      k          <- k + 1L
      acc[[k]]   <- tibble::tibble(
        true_dimension = p,
        sample_n       = sample_n[j],
        est_dimension  = est_dim,
        simulation     = sim
      ) }}
  dplyr::bind_rows(acc)}
parallel::stopCluster(cl)
# Save the result. ----
saveRDS(ranks_sim, file = "./Datas/intro_simulation.rds")