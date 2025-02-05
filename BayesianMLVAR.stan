data {
  int<lower=1> S;               // Number of subjects
  int<lower=1> K;               // Number of variables (here, K = 3)
  int<lower=1> N;               // Total number of observations across all subjects
  int<lower=1> subject[N];      // Subject ID for each observation (sorted by subject and time)
  vector[K] X[N];               // Observed data: time-series data for each time point
  int<lower=0,upper=1> missing_mask[N, K]; // Mask: 1 if X[i,k] is missing, 0 otherwise
  int<lower=1> start[S];        // start[s]: starting index for subject s
  int<lower=1> end[S];          // end[s]: ending index for subject s
  int<lower=0> M;  // number of missing entries
  // This was not needed: int<lower=1, upper=N> missing_idx[M];
  int<lower=1, upper=K> missing_var[M];
}
parameters {
  // 1. VAR(1) model parameters
  vector[K] c;              // Intercept vector
  matrix[K, K] A;           // VAR(1) coefficient matrix
  
  // 2. Covariance matrix for residuals
  cholesky_factor_corr[K] L_Omega;  // Cholesky factor of the correlation matrix
  vector<lower=0>[K] sigma;         // Standard deviations
  
  // 3. Random intercepts for each subject
  matrix[S, K] subject_intercept;  // subject-specific intercept deviations
  
  // 4. Hyperparameters for the random intercepts
  vector<lower=0>[K] tau_subj;     // Standard deviation for subject intercepts
  
  // 5. Missing values vector, as a parameter.
  vector<lower=1>[M]missing_vals;
  }
transformed parameters {
  
  // 1. Create the covariance matrix.
  matrix[K, K] Omega;--------
  Omega = multiply_lower_tri_self_transpose(L_Omega) .* (sigma * sigma');

  // 2. Handle missing data by constructing a complete data matrix X_full.
  matrix[K, N] X_full;
  {
    int m_counter = 1;        // Counter for indexing into missing_vals
  for (i in 1:N) {
    for (k in 1:K) {
      if (missing_mask[i, k] == 1) {
        // Use imputed value when data is missing
        X_full[k, i] = missing_vals[ m_counter ];
        m_counter += 1; //Move to the next missing value.
      } else {
        // Use observed data when available
        X_full[k, i] = X[i, k];
      }
    }
  }
  }
}
model {
  
  // 1. Priors for parameters
  c ~ normal(2.5, 1);                           // Prior for global intercept
  to_vector(A) ~ normal(0, 0.5);                 // Prior for VAR coefficients
  L_Omega ~ lkj_corr_cholesky(2);              // Prior for correlation matrix
  sigma ~ exponential(1);                      // Prior for residual standard deviations
  tau_subj ~ exponential(1);                   // Prior for subject intercept SDs
  
  // 2. Subject-specific intercepts
  for (s in 1:S) {
    subject_intercept[s] ~ normal(0, tau_subj);
  }

  
  // 3. Model the missing values:
  for(i in 1:M) {
    if(missing_var[ i ]  == 1 ) {
      missing_vals[ i ] ~ normal(2.5, 0.5);
    } if (missing_var[ i ]  == 2 ) {
      missing_vals[ i ] ~ normal(2.5, 0.5);
    } else {
      missing_vals[ i ] ~ normal(2.5, 0.5);
    }
  }
  
  // 4. VAR(1) model loop over subjects and time points
  for (s in 1:S) {
    // Loop over time for each subject
    for (t in start[s]:end[s]) {
      if (t == start[s]) {
        // Initial condition for each subject
        // You might want a separate model for initial state; here we use the intercept directly.
        X_full[, t] ~ multi_normal(c + subject_intercept[s]', Omega);
      } else {
        // VAR(1) relation: X_t ~ MVN(c + subject-specific adjustments + A * X_{t-1}, Omega)
        X_full[, t] ~ multi_normal(c + subject_intercept[s]' + A * X_full[, t - 1], Omega);
      }
    }
  }
}