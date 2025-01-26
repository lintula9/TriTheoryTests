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
  int<lower=1> beep[N];         // Time of day, for each observation. (4 = evening, 1 = morning.)
  int<lower=1> nbeeps;          // Count of different times of day.
}

parameters {
  // 1. VAR(1) model parameters
  matrix[K, K] A;           // VAR(1) coefficient matrix
  matrix[K, K] B;        // Evening to morning transition adjustment matrix.
  
  // 2. Correlation matrix for residuals
  cholesky_factor_corr[K] L_corr;  // Cholesky factor of the correlation matrix
  vector<lower=0>[K] X_star_innovation_sd; // Scale for the (independent) innovations.

  // 3. Random intercepts for each subject
  matrix[S, K] subject_intercept_raw;  // subject-specific intercept deviations
  vector<lower=0>[K] subject_intercept_sd;     // Standard deviation for subject intercepts
  
  // 4. X_star as latent symptoms.
  matrix[K,N] X_star_innovation; // The non-deterministic part of X_star.
  matrix[K,S] X_star_zero;   // Initial state for each subject at t=0.
  matrix<lower=0>[K,S] subject_innovation_sd; // Subject specific scale for X_star innovations.

  // 5. Probit measurement model cutoffs:
  ordered[cutpoint_count] cutpoints[K];
  
  // 8. Time of day parameter:
  vector[K] time_of_day_intercept[nbeeps];
  }
  
transformed parameters {

  // 3. Subject intercepts. Better posterior geometry with non-centralized priors.
  vector[K] subject_intercept[S];
  for(s in 1:S){
    for(k in 1:K){
  subject_intercept[s, k] = subject_intercept_raw[s,k] * subject_intercept_sd[k];
  }}
  
  // 99. c = 0 (written in for compatibility of all syntax, though redundant).
  vector[K] c;
  c = rep_vector(0, K);
  
  // 4. Latent symptoms X_star as a transformed parameter: 
  // Improved posterior geometry with non-centralized priors
  matrix[K,N] X_star;
  for (s in 1:S) {
    matrix[K,K] L;
    L = diag_pre_multiply(subject_innovation_sd[,s], L_corr);
    // Loop over time for each subject
    for (t in start[s]:end[s]) {
            if(t==start[s]){     // diag_pre_multiply creates the full cholesky factor here.
          X_star[,t] = c + subject_intercept[s] 
                         + A * X_star_zero[,s] 
                         + time_of_day_intercept[beep[t]];
    } else {
      // Scale of X_star is problematic for identifiability.
          X_star[,t] = c + subject_intercept[s] 
                         + A * X_star[,t-1] 
                         + L * ( X_star_innovation[,t] .* X_star_innovation_sd ) 
                         + time_of_day_intercept[beep[t]]
                         + (beep[t] == 4) * B * X_star[,t-1];
    }}}
}

model {
  
  // 1. Priors for VAR model parameters
  to_vector(A) ~ normal(0, 0.5);               // Prior for VAR coefficients
  to_vector(B) ~ normal(0, 0.5);               // Prior for evening-to-morning adjustment

  L_corr ~ lkj_corr_cholesky(1);              // Prior for correlation matrix
  to_vector(X_star_innovation) ~ normal(0,0.5); 
  // After mixing with L: X -> LX, we obtain multinormal variables.
  // Scale of X_star_innovations is difficult to identify. We'll use a vague prior.
  to_vector(X_star_innovation_sd) ~ normal(0.707,1);
  
  // 3. Priors for subject-specific intercepts
  subject_intercept_sd ~ normal(1,0.1);                   
  // Prior for subject intercept SDs, strong prior for identifiability.
  for (s in 1:S) {
    subject_intercept_raw[s] ~ normal(0, 0.1);
  }
  
  // 4. Latent symptoms at t = 0 and subject specific scale.
  for (s in 1:S) {
    X_star_zero[,s] ~ multi_normal(mu0, Sigma0); 
    }
  to_vector(subject_innovation_sd) ~ normal(1,0.1); 
  // Tight prior, for identifiability. 

  // 5. Priors for the cutoffs.
  for (k in 1:K) {
    for(j in 1:cutpoint_count){
    cutpoints[k][j] ~ normal(cutpoint_prior_locations[j], 0.1); 
    // Empirical CDFs of data justify strong prior, assuming X_star is standard multinormal.
    // The problem is that cutpoints correlate with scaling of X_star.
    }}
  
  // 8. Time of day prior:
  for(time in 1:nbeeps){
    time_of_day_intercept[time] ~ normal(0,0.5); }
  
  // 6. Model target likelihood and variables. OLD:

  for (k in 1:K){
    for (n in 1:N) {
      if(missing_mask[n,k] == 0){
       // Note, that X is N times K, whereas X_star is K times N.
       target += ordered_probit_lpmf( X[n,k] | X_star[k,n], cutpoints[k]);
        }}}
  // Possibly change something that does not require missing_mask;
  // e.g., vectorize whole data and create indicator variables for n (row), k(column.
  // Then we have:
  // for (m in 1:M) {
  //  target += ordered_probit_lpmf(X[m] | X_star[observed_k[m], observed_n[m]], cutpoints[observed_k[m]]);
  //  }

  }

generated quantities {

  // 1. Create the residual covariance matrix ( for interpretation)
  matrix[K, K] Omega;
  Omega = multiply_lower_tri_self_transpose(L_corr);
  // Here we could implement the quick and dirty 'close' indistinguishable VAR(1).
  
}
