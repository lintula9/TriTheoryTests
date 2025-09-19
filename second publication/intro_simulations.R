# Install packages. ----
pkgs <- c("MTS", "ggplot2", "gridExtra", "rstan", "progress", "tidyr", "tibble", "dplyr")
for (pkg in pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg)}  else {
      library(pkg, character.only = T)}}
# Parameters.
n_obs <- 10000       # kept observations
burn  <- 200         # burn-in dropped by skip
p_max <- 40
frequencies <- seq(n_obs / 50, n_obs / 2, length.out = 1000)
ranks <- expand.grid(3:p_max, frequencies) |> 
         setNames(c("true_dimension", "sampling_freq"))
# Pb
pb <- progress::progress_bar$new(total = p_max-2 )
for( p in 3:p_max) {
  # Random VAR time series.
  A <- matrix(0,p,p)
  A[upper.tri(A, diag = T)] <- runif((p*(p-1) / 2) + p, -0.99, 0.99)
  true_simu <- MTS::VARMAsim(
    nobs  = n_obs, arlags = 1, 
    phi   = A, 
    sigma = diag(p), skip=burn)
  pb2 <- progress::progress_bar$new(total = length(frequencies))
  for( frq in frequencies ){
    # Pick samples.
    obs_series <- true_simu$series[ seq(1, n_obs, by = frq), ]
    # Save rank of the covariance.
    result     <- cov(obs_series) |> eigen() |> 
                  (\(x) x$values > 1e-06)()  |> sum()
    ranks[ ranks$true_dimension == p & 
           ranks$sampling_freq  == frq, 
           c("est_dimension") ] <- result
    pb2$tick()}
  pb2$terminate()
  pb$tick()}
pb$terminate() ; gc()
# Plot the result.
ggplot(data = ranks,
       aes(x     = (n_obs / sampling_freq) / n_obs, 
           y     = est_dimension,
           color = true_dimension -  est_dimension)) + 
  geom_point() +
  theme_bw() + 
  ylab("Estimated dimension") + 
  xlab("Proportion of equidistant observations") +
  scale_y_continuous(breaks = 1:p_max) +
  scale_color_gradientn(
    colours = hcl.colors(2),
    values  = scales::rescale(c(min(ranks$true_dimension),
                                max(ranks$true_dimension)))
  )
