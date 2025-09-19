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
p_max  <- 40
sample_n <- round(seq(2, 5000, length.out = 1000)) |> unique()
ranks <- expand.grid(3:p_max, sapply(3:p_max, \(x) max(3:x)), sample_n) |>
  setNames(c("true_dimension","obs_dimension","sample_n"))
# Parallel backend. ----
cl <- parallel::makeCluster(max(1, parallel::detectCores() - 1))
doParallel::registerDoParallel(cl)
# Compute (parallel over p). ----
ranks_est <- foreach::foreach(
  p = 3:p_max,
  .combine  = "rbind"
) %dopar% {
  # Simulate VAR(1) for this p
  A <- matrix(0, p, p)
  
  # Could add another case, when A is not a DAG.
  
  A[upper.tri(A, diag = TRUE)] <- runif((p * (p - 1) / 2) + p, -0.99, 0.99)
  true_simu <- MTS::VARMAsim(
    nobs = n_obs, arlags = 1,
    phi = A, sigma = diag(p), skip = burn
  )
  idx_list <- lapply(sample_n, function(n) unique(ceiling(seq(1, n_obs, length.out = n))))
  acc <- vector("list", length = (p - 2L) * length(sample_n))
  k <- 0L
  for (p_obs in 3:p) {
    for (j in seq_along(sample_n)) {
      idx <- idx_list[[j]]
      obs_series <- true_simu$series[idx,  1:p_obs, drop = FALSE]
      est_dim    <- sum(eigen(cov(obs_series))$values > .Machine$double.eps)
      k          <- k + 1L
      acc[[k]]   <- tibble::tibble(
        true_dimension = p,
        obs_dimension  = p_obs,
        sample_n       = sample_n[j],
        est_dimension  = est_dim
      ) }}
  dplyr::bind_rows(acc)}
parallel::stopCluster(cl)
# Merge with your grid (fix the join keys). ----
ranks <- ranks |>
  dplyr::select(-dplyr::any_of("est_dimension")) |>
  dplyr::left_join(ranks_est,
                   by = c("true_dimension","obs_dimension","sample_n"))
# Compute statistics. ----
ranks <- ranks %>%
  mutate( dim_diff = obs_dimension - est_dimension )
plotdf <- ranks %>% group_by(obs_dimension) %>%
  summarise(
    mean_diff = mean(dim_diff, na.rm = T),
    obs_dim   = first(obs_dimension)
  )
plot(plotdf$mean_diff, plotdf$obs_dim)
