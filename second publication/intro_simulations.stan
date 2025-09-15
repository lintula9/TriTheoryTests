// Model X(t) as a latent process with gaussian process noise terms.
data {
  int<lower=2>       N;            // number of time points
  int<lower=1>       p;            // dimension
  matrix[N,p]        Y_case3;      // observed series, scaled to mean 0 and SD 1.
}
parameters {
  // State-space params
  vector[p]    c;                          // intercept
  matrix[p, p] A;                          // coefficient matrix
  vector[p]               X_start;         // initial latent state X[1]
  array[N-1] vector[p]    z;               // standard-normal shocks for t = 2..N
  vector<lower=0>[p]      sigma;           // std devs
  cholesky_factor_corr[p] L_Omega;         // Cholesky of correlation
  // Likelihood param(s)
  vector<lower=0>[p]      residual_sd;     // Residual SD.
}
transformed parameters {
  matrix[p, p] L_Sigma = diag_pre_multiply(sigma, L_Omega);
  array[N]   vector[p] X;                  // latent state
  // Latent process
  X[1] = X_start;
  for (t in 2:N)
    X[t] = c + A * X[t - 1] + L_Sigma * z[t - 1];
}
model {
  c                  ~ normal(0,0.5);
  to_vector(A)       ~ normal(0,0.1);
  for(t in 1:(N - 1))
    to_vector(z[t])  ~ normal(0,  1);
  to_vector(X_start) ~ normal(0,  1);
  sigma              ~ exponential(1);
  // Favor diagonal innovation covariance.
  L_Omega            ~ lkj_corr_cholesky(2); 
  // Likelihood for the observations
  residual_sd       ~ exponential(1);
  for (k in 1:p)
    Y_case3[, p ~ normal(X[, p], residual_sd[p]);
}
