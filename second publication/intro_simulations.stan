// Model X(t) as a latent process with gaussian process noise terms.
data {
  int<lower=2>       N;           // number of time points
  int<lower=1>       p;           // dimension
  matrix[N,p]        Y_case3;     // observed series, scaled to mean 0 and SD 1.
}
parameters {
  // State-space params
  vector[p]               c;               // intercept
  matrix[p,p]             A;               // coefficient matrix
  vector<lower=0>[p]      sigma;           // std devs
  cholesky_factor_corr[p] L_Omega;         // Cholesky of correlation
}
transformed parameters {

  matrix[p, p] L_Sigma = diag_pre_multiply(sigma, L_Omega);
  matrix[p, p] Sigma   = L_Sigma * L_Sigma';
}
model {
  c                  ~ normal(0,1);
  to_vector(A)       ~ normal(0,1);
  sigma              ~ exponential(1);
  L_Omega            ~ lkj_corr_cholesky(1); 
  // Likelihood
  for (n in 2:N)
    Y_case3[n, 1:p]  ~ multi_normal(c + A * (Y_case3[n-1, 1:p])', Sigma);
}
