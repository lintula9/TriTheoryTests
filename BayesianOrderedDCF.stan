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
  int<lower=1> cutpoint_count;
  real cutpoint_prior_locations[cutpoint_count];
}
parameters {
  // 1. D-CF model parameters
  // Will be set to 0. real c;              // Intercept vector.
  real psi;                 // CF autoregression coefficient.
  vector<lower=0>[K] Lambda; // Factor loadings are assumed positive, to better identify the model.
  
  // 2. Covariance matrix for residuals: NOT NEEDED in D-CF model.
  // This will be produced in the transformed parameters block as (1-psi^2) * Lambda * Lambda'.
  // cholesky_factor_corr[K] L_Omega;  // Cholesky factor of the correlation matrix
  // vector<lower=0>[K] sigma;         // Standard deviations
  
  // 3. Random intercepts for each subject
  vector[S] subject_intercept_raw;  // subject-specific intercept deviations
  real<lower=0> tau_subj;     // Standard deviation for subject intercept.
  
  // 5. Missing values vector, as a parameter.
  // NOT NEEDED in the ordered case. X* 'deals' with missingness.
  
  // 6. UNIQUE to ordered: X* as the latent standard multinormal variables
  // Note the K times N structure, with time intervals for each subject
  // on the second dimension.
  // matrix[K,N] X_star; We wont use X_star, to make the notation more distinct.
  vector[N] eta;
  
  // 7. UNIQUE to ordered: X* (time and subject invariant) thresholds:
  ordered[cutpoint_count] cutpoints[K];
  }
  
transformed parameters {
  // Better posterior geometry with non-centralized priors.
  vector[S] subject_intercept;
  for(s in 1:S){
  subject_intercept[s] = subject_intercept_raw[s] * tau_subj;
  }
  
  // 2. Handle missing data by constructing a complete data matrix X_full.
  // NOT NEEDED in ordered case. We only use the X*, which is missing anyhow.
  
  // c = 0 (written in for compatibility of all syntax, though redundant).
  real c;
  c = 0;
  }
model {
  // 1. Priors for parameters
  // Due to identifiability issues, c will be set to 0.
  // Subject-specific drifts are allowed in subject_intercept_raw
  // c ~ normal(0, 1);                           // Prior for global intercept
  psi ~ normal(0, 0.5);                 // Prior for DCF autoregression coefficient.
  Lambda ~ normal(0, 5);                // Half-normal, weak, prior.
  // sigma ~ exponential(1);  NOT NEEDED for DCF // Prior for residual standard deviations. Not needed as Lambda, psi determine Omega.
  // 2. Subject-specific intercepts
  tau_subj ~ normal(0,1); // Prior for subject intercept SDs
  for (s in 1:S) {
    subject_intercept_raw[s] ~ normal(0, 0.1); // Raw intercept, scaled by tau_subj.
  }
  // 3. Model the missing values: NOT NEEDED for DCF since X* is missing anyhow!
  // 3. UNIQUE to ordered: Declare priors for the cutoffs.
  for (k in 1:K) {
    for(tau in 1:cutpoint_count)
    cutpoints[k] ~ normal(cutpoint_prior_locations[tau], 1);
    }
  // 4. VAR(1) model loop over subjects and time points
  for (s in 1:S) {
    // Loop over time for each subject
    for (t in start[s]:end[s]) {
      // UNIQUE TO ORDERED: X_full is modeled as a multivariate normal, with 
      // autoregressive structure.
            if(t==start[s]){
          eta[t] ~ normal(c + subject_intercept[s], 1);
    } else {
          eta[t] ~ normal(c + subject_intercept[s] + psi*eta[t - 1], 1);
    }
      for (k in 1:K){
        // Skip missing values - though they are modeled in the latent eta.
      if(missing_mask[t,k] == 0){
      // Note, that we use now for the log likelihood.
        target += ordered_probit_lpmf( X[t,k] | Lambda[k]*eta[t], cutpoints[k]);
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
        log_lik[n, k] = ordered_probit_lpmf( X[n,k] | Lambda[k]*eta[n], cutpoints[k]);
      } else {
        log_lik[n, k] = 0; // or NA, but zero won't affect WAIC/LOO calculations
      }
    }
  }
  
  matrix[K,K] A; // Map to A matrix (VAR(1) model presentation)
  A = psi * Lambda * Lambda'; // We can input here some estimated A matrix to compute the closest A.
  matrix[K,K] Z; // Map to Z matrix (VAR(1) model presentation)
  Z = (1-psi^2) * Lambda * Lambda';
}
