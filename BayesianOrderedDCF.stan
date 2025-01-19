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
  vector[K] c;              // Intercept vector.
  real psi;                 // CF autoregression coefficient.
  real<lower=0> Lambda_first;
  vector[K-1] Lambda_rest;         // Factor loadings.
  
  // 2. Covariance matrix for residuals: NOT NEEDED in D-CF model.
  // This will be produced in the transformed parameters block as (1-psi^2) * Lambda * Lambda'.
  // cholesky_factor_corr[K] L_Omega;  // Cholesky factor of the correlation matrix
  // vector<lower=0>[K] sigma;         // Standard deviations
  
  // 3. Random intercepts for each subject
  matrix[S, K] subject_intercept;  // subject-specific intercept deviations
  
  // 4. Hyperparameters for the random intercepts
  vector<lower=0>[K] tau_subj;     // Standard deviation for subject intercepts
  
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
  
  // 1. Create the covariance matrix.
  matrix[K, K] Omega;
  Omega = (1-psi^2) * Lambda * Lambda';

  // 2. Handle missing data by constructing a complete data matrix X_full.
  // NOT NEEDED in ordered case. We only use the X*, which is missing anyhow.
  
  // 3. To solve the sign indeterminancy of factor loadings, we construct them separately here.
  vector[K] Lambda;
  Lambda[1] = Lambda_first;
  for (k in 2:K)
    Lambda[k] = Lambda_rest[k-1];
}
model {
  
  // 1. Priors for parameters
  c ~ normal(0, 1);                           // Prior for global intercept
  psi ~ normal(0, 0.5);                 // Prior for DCF autoregression coefficient.
  Lambda_first ~ normal(0,0.5);              // Prior for correlation matrix
  Lambda_rest ~ normal(0,0.5);              // Prior for correlation matrix
  // sigma ~ exponential(1);                      // Prior for residual standard deviations. Not needed as Lambda, psi determine Omega.
  tau_subj ~ exponential(1);                   // Prior for subject intercept SDs
  
  // 2. Subject-specific intercepts
  for (s in 1:S) {
    subject_intercept[s] ~ normal(0, tau_subj);
  }

  // 3. Model the missing values: NOT NEEDED since X* is missing anyhow!
  // 3. UNIQUE to ordered: Declare priors for the cutoffs.
  for (k in 1:K) {
    for(tau in 1:cutpoint_count)
    cutpoints[k] ~ normal(cutpoint_prior_locations[tau], 5);
    }

  // 4. VAR(1) model loop over subjects and time points
  for (s in 1:S) {
    // Loop over time for each subject
    for (t in start[s]:end[s]) {
      // UNIQUE TO ORDERED: X_full is modeled as a multivariate normal, with 
      // autoregressive structure.
            if(t==start[s]){
          eta[t] ~ multi_normal(c + subject_intercept[s]', Omega);
    } else {
          eta[t] ~ multi_normal(c + subject_intercept[s]' + psi*eta[t - 1], Omega);
    }
      for (k in 1:K){
        // Skip missing values - though they are modeled in the latent eta.
      if(missing_mask[t,k] == 0){
      // Note, that we use now for the log likelihood.
        target += ordered_probit_lpmf( X[t,k] | eta[t], cutpoints[k]);
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
}
