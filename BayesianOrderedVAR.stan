data {
  int<lower=1> S;               // Number of subjects
  int<lower=1> K;               // Number of variables (here, K = 3)
  int<lower=1> N;               // Total number of observations across all subjects
  int<lower=1> subject[N];      // Subject ID for each observation (sorted by subject and time)
  int<lower=1,upper=5> X[N,K];               // Observed data: time-series data for each time point
  int<lower=0,upper=1> missing_mask[N, K]; // Mask: 1 if X[i,k] is missing, 0 otherwise
  int<lower=1> start[S];        // start[s]: starting index for subject s
  int<lower=1> end[S];          // end[s]: ending index for subject s
  int<lower=1> cutpoint_count; // Note that we assume (know) that category count is the same for all vars.
  real cutpoint_prior_locations[cutpoint_count];
  vector[K] mu0;              // Mean for initial states
  cov_matrix[K] Sigma0;       // Covariance for initial states
}
parameters {
  // 1. VAR(1) model parameters
  matrix[K, K] A;           // VAR(1) coefficient matrix
  
  // 2. Covariance matrix for residuals
  cholesky_factor_corr[K] L_Omega;  // Cholesky factor of the correlation matrix
  vector<lower=0>[K] sigma;         // Standard deviations
  
  // 3. Random intercepts for each subject
  matrix[S, K] subject_intercept_raw;  // subject-specific intercept deviations
  vector<lower=0>[K] tau_subj;     // Standard deviation for subject intercepts
  
  // 4. X_star as the latent, true, symptom structure.
  matrix[K,N] X_star;
  vector[K] X_star_zero[S];   // Initial state for each subject at t=0

  
  // 5. UNIQUE to ordered: X* (time and subject invariant) cutoffs:
  ordered[cutpoint_count] cutpoints[K];
  }
  
  
transformed parameters {
  // 2. Create the residual covariance matrix.
  matrix[K, K] Omega;
  Omega = multiply_lower_tri_self_transpose(L_Omega) .* (sigma * sigma');
  
  // 3. Subject intercepts. Better posterior geometry with non-centralized priors.
  vector[K] subject_intercept[S];
  for(s in 1:S){
    for(k in 1:K){
  subject_intercept[s, k] = subject_intercept_raw[s,k] * tau_subj[k];
  }}
  
  // 99. c = 0 (written in for compatibility of all syntax, though redundant).
  vector[K] c;
  c = rep_vector(0, K);
}
model {
  
  // 1. Priors for VAR model parameters
  // NOT needed, set to 0. c ~ normal(0, 1);                           // Prior for global intercept
  to_vector(A) ~ normal(0, 0.2);                 // Prior for VAR coefficients
  L_Omega ~ lkj_corr_cholesky(2);              // Prior for correlation matrix
  sigma ~ exponential(1);                      // Prior for residual standard deviations
  
  // 3. Priors for subject-specific intercepts
  tau_subj ~ normal(0,1);                   // Prior for subject intercept SDs
  for (s in 1:S) {
    subject_intercept_raw[s] ~ normal(0, 0.1);
  }

  // 5. Priors for the cutoffs.
  for (k in 1:K) {
    for(j in 1:cutpoint_count)
    cutpoints[k] ~ normal(cutpoint_prior_locations[j], 0.5); // or some other weakly informative prior
    }
  
  // 6. VAR(1) model loop over subjects and time points
  for (s in 1:S) {
    // Initiate X_star at t = 0.
    X_star_zero[s] ~ multi_normal(mu0, Sigma0); 
    // Loop over time for each subject
    for (t in start[s]:end[s]) {
      // UNIQUE TO ORDERED: X_full is modeled as a multivariate normal, with 
      // autoregressive structure.
            if(t==start[s]){
          X_star[, t] ~ multi_normal(c + subject_intercept[s] + A * X_star_zero[s], Omega);
    } else {
          X_star[, t] ~ multi_normal(c + subject_intercept[s] + A * X_star[, t - 1], Omega);
    }
      for (k in 1:K){
      if(missing_mask[t,k] == 0){
      // Note, that X is N times K, whereas X_star is K times N, 
      // hence the time index is placed differently.
        target += ordered_probit_lpmf( X[t,k] | X_star[k,t], cutpoints[k]);
        }
        }
    }
  }
  }
generated quantities {
  // We'll store the pointwise log-likelihood for each observation 
  // as a matrix of size N x K or as a vector if observations are univariate.
  matrix[N, K] log_lik;
  
  for (n in 1:N) {
    for (k in 1:K) {
      if (missing_mask[n, k] == 0) {
        // Compute log-probability only for non-missing observations.
        log_lik[n, k] = ordered_probit_lpmf(X[n,k] | X_star[k,n], cutpoints[k]);
      } else {
        log_lik[n, k] = 0; // or NA, but zero won't affect WAIC/LOO calculations
      }
    }
  }
  
  // Here we could implement the quick and dirty 'close' indistinguishable VAR(1).
  

}
