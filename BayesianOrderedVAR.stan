data {
  int<lower=1> S;               // Number of subjects
  int<lower=1> K;               // Number of variables (here, K = 3)
  int<lower=1> N;               // Total number of observations across all subjects
  int<lower=1> subject[N];      // Subject ID for each observation (sorted by subject and time)
  int<lower=1,upper=5> X[N,K];               // Observed data: time-series data for each time point
  int<lower=0,upper=1> missing_mask[N, K]; // Mask: 1 if X[i,k] is missing, 0 otherwise
  int<lower=1> start[S];        // start[s]: starting index for subject s
  int<lower=1> end[S];          // end[s]: ending index for subject s
  // Number of missing entries not needed in ordered case:
  // int<lower=0> M;  // number of missing entries
  // This was not needed: int<lower=1, upper=N> missing_idx[M];
  // This is not needed: int<lower=1, upper=K> missing_var[M];
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
  // NOT NEEDED in the ordered case. X* 'deals' with missingness.
  
  // 6. UNIQUE to ordered: X* as the latent standard multinormal variables
  // Note the K times N structure, with time intervals for each subject
  // on the second dimension.
  matrix[K,N] X_star;
  
  // 7. UNIQUE to ordered: X* (time and subject invariant) thresholds:
  ordered[4] cutpoints[K];
  }
  
  
transformed parameters {
  
  // 1. Create the covariance matrix.
  matrix[K, K] Omega;
  Omega = multiply_lower_tri_self_transpose(L_Omega) .* (sigma * sigma');

  // 2. Handle missing data by constructing a complete data matrix X_full.
  // NOT NEEDED in ordered case. We only use the X*, which is missing anyhow.
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

  // 3. Model the missing values: NOT NEEDED since X* is missing anyhow!
  // 3. UNIQUE to ordered: Declare priors for the cutoffs.
  for (k in 1:K) {
    cutpoints[k] ~ normal(0, 5);                // or some other weakly informative prior
    }

  // 4. VAR(1) model loop over subjects and time points
  for (s in 1:S) {
    // Loop over time for each subject
    for (t in start[s]:end[s]) {
      // UNIQUE TO ORDERED: X_full is modeled as a multivariate normal, with 
      // autoregressive structure.
            if(t==start[s]){
          X_star[, t] ~ multi_normal(c + subject_intercept[s]', Omega);
    } else {
          X_star[, t] ~ multi_normal(c + subject_intercept[s]' + A * X_star[, t - 1], Omega);
    }
      for (k in 1:K){
      if(missing_mask[t,k] == 0){
      // Note, that X is N times K, whereas X full is K times N, 
      // hence the time index is placed differently.
        target += ordered_probit_lpmf( X[t,k] | X_star[k,t], cutpoints[k])}
        }
    }
  }
  }