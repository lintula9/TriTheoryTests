// Model X(t) as a latent process with gaussian process noise terms.
data {
  int<lower=2>       N;          // number of time points
  int<lower=1>       p;          // dimension
  matrix[N,p]        Y;    // observed series, scaled to mean 0 and SD 1.
}
parameters {
  // VAR process params
  vector[p]               c;               // intercept
  matrix[p,p]             A;               // coefficient matrix
  // Serial noise term params
  vector[p]               eps_ar;      // noise (residual) autocorrelations
  // Residual distribution params
  vector<lower=0>[p]      sigma;           // std devs
  cholesky_factor_corr[p] L_Omega;         // Cholesky of correlation
}
transformed parameters {
  // Process prediction
  matrix[N-1,p]          Y_hat;
  for(n in 2:N)
    Y_hat[n-1, 1:p]    = (c + A * (Y[n-1, 1:p])')';
  // Likelihood covariance
  matrix[p, p] L_Sigma = diag_pre_multiply(sigma, L_Omega);
  // noise term (autocorrelated residual)
  matrix[N,p]          eps;
  eps[1,   1:p]      = rep_row_vector(0,p);     // Initialize eps.
  eps[2:N, 1:p]      = Y[2:N, 1:p] - Y_hat;     // Define eps as resid.
  matrix[p, p] A_eps = diag_matrix(eps_ar); // Diagonal AR for eps.
}
model {
  // VAR process
  c                  ~ normal(0,1);
  to_vector(A)       ~ normal(0,0.5);
  sigma              ~ exponential(1);
  L_Omega            ~ lkj_corr_cholesky(1); 
  
  // noise (residual) likelihood
  eps_ar         ~ normal(0,0.5);
  for(n in 3:N) {
    vector[p] mu = A_eps * (eps[n-1,1:p])';
    target      += multi_normal_cholesky_lpdf(eps[n,1:p] | mu,L_Sigma);}
  
  }
generated quantities{
  // Residual covariance
  matrix[p,p] Sigma   = L_Sigma * L_Sigma';
  
}