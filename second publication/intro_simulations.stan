// Model X(t) as a latent process with gaussian process noise terms.

data {
  int<lower=2>  N;            // number of time points
  int<lower=1>  p;            // dimension
  matrix[N, p]  Y_case1;      // observed series
}

functions {
  // Squared Exponential (RBF) kernel
  matrix cov_se(vector x, real alpha, real ell) {
    int N = num_elements(x);
    matrix[N, N] K;
    for (i in 1:N) {
      K[i, i] = square(alpha);
      for (j in (i + 1):N) {
        real r = (x[i] - x[j]) / ell;
        real v = square(alpha) * exp(-0.5 * r * r);
        K[i, j] = v;
        K[j, i] = v;
      }
    }
    return K;
  }

parameters {
  // VAR(1) params
  vector[p]    c;             // intercept
  matrix[p, p] A;             // coefficient matrix

  // State-space params
  vector[p]               X_star;          // initial latent state X[1]
  array[N-1] vector[p]    z;               // standard-normal shocks for t = 2..N

  // Innovation covariance (Cholesky parameterization)
  vector<lower=0>[p]      sigma;           // std devs
  cholesky_factor_corr[p] L_Omega;         // Cholesky of correlation
  
  // Noise.
  matrix[N,p]             epsilon;         // Noise matrix.
  
}

transformed parameters {
  matrix[p, p] L_Sigma = diag_pre_multiply(sigma, L_Omega);

  array[N-1] vector[p] epsilon;         // ε_t for transitions 2..N
  array[N]   vector[p] X;               // latent states

  // Scale the non-centered shocks
  for (t in 1:(N - 1))
    epsilon[t] = L_Sigma * z[t];        // ε_t = L_Sigma z_t

  // Latent process
  X[1] = X_star;
  for (t in 2:N)
    X[t] = c + A * X[t - 1] + epsilon[t - 1];
}