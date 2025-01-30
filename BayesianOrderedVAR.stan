data {
  // General integers.
  int<lower=1> S;               // Number of subjects
  int<lower=1> K;               // Number of variables (here, K = 3)
  int<lower=1> N;               // Total number of observations across all subjects
  array[N] int<lower=1> subject;      // Subject ID for each observation (sorted by subject and time)
  int<lower=1> N_obs;        // Number of observations

  // Matrix form data.
  array[N,K] int<lower=1,upper=4> X;   // Observed data: matrix form. (Missing values contained!)
  array[N,K] int<lower=0,upper=1> missing_mask; // Mask: 1 if X[i,k] is missing, 0 otherwise. Matrix form.

  // Subject indexing.
  array[S] int<lower=1>  start;        // start[s]: starting index for subject s
  array[S] int<lower=1>  end;          // end[s]: ending index for subject s
  
  // Cutpoints for probit measurement model.
  int<lower=1> cutpoint_count; // Note that we assume (know) that category count is the same for all vars.
  array[cutpoint_count] real cutpoint_prior_locations;
  
  // Initial, t = 0, unobserved state.
  array[K] real mu0;              // Mean for initial states
  cov_matrix[K] Sigma0;       // Covariance for initial states
  
  // Time of day variables.
  array[N] int<lower=1> beep;         // Time of day, for each observation. (4 = evening, 1 = morning.)
  int<lower=1> nbeeps;          // Count of different times of day.
}

parameters {
  // 1. VAR(1) model parameters
  matrix[K, K] A;           // VAR(1) coefficient matrix

  // 2. Correlation matrix for innovations.
  cholesky_factor_corr[K] L_corr;  // Cholesky factor of the correlation matrix

  // 3. Random intercepts for each subject
  matrix[S, K] subject_intercept_raw;  // subject-specific intercept deviations
  array[S] real<lower=0> subject_intercept_sd;     // Standard deviation for subject intercepts

  // 4. X_star as latent symptoms.
  // Note: the scale of the latent X* symtpoms is not identified, thus it is user defined.
  matrix[K,N] X_star_innovation; // The non-deterministic part of X_star.
  matrix[K,S] X_star_zero;   // Initial state for each subject at t=0.

  // 5. Probit measurement model cutoffs:
  array[K] ordered[cutpoint_count] cutpoints;
  
  // 8. Time of day parameter:
  array[nbeeps] vector[K] time_of_day_intercept;
  }
  
transformed parameters {

  // 3. Subject intercepts. Better posterior geometry with non-centralized priors.
  matrix[K,S] subject_intercept;
  for(s in 1:S){
    for(k in 1:K){
  subject_intercept[k,s] = subject_intercept_raw[s,k] * subject_intercept_sd[s];
  }}
  
  // 99. c = 0 (written in for compatibility of all syntax, though redundant).
  vector[K] c;
  c = rep_vector(0, K);
  
  // 4. Latent symptoms X_star as a transformed parameter: 
  // Improved posterior geometry with non-centralized priors
  matrix[K,N] X_star;
  for (s in 1:S) {
    // Loop over time for each subject
    for (t in start[s]:end[s]) {
            if(t==start[s]){    
          X_star[,t] = c + A * X_star_zero[,s] 
                         + time_of_day_intercept[beep[t]];
    } else {
          X_star[,t] = c + subject_intercept[,s] 
                         + A * X_star[,t-1] 
                         + L_corr * X_star_innovation[,t]  
                         + time_of_day_intercept[beep[t]];
    }}}
    
}

model {
  
  // 1. Priors for VAR model parameters
  to_vector(A) ~ normal(0,0.2);               // Prior for VAR coefficients
  // Since it is untypical for VAR (especially for stationary VAR) to have coefficients
  // with absolute values much larger than 0.5, we will set the prior tight. This will still
  // allow for large coefficients to appear, but, also stabilizes the scales better as well.

  L_corr ~ lkj_corr_cholesky(1);              // Prior for correlation matrix
  to_vector(X_star_innovation) ~ std_normal(); 
  // After mixing with L: X -> LX, we obtain multinormal, possibly correlated, variables.

  // 3. Priors for subject-specific intercepts
  subject_intercept_sd ~ normal(1,0.1); // Hierarchical prior for the sd                   
  // Strong prior for identifiability.
  for (s in 1:S) {
    subject_intercept_raw[s] ~ normal(0, 0.1);
  }
  
  // 4. Latent symptoms at t = 0 and subject specific scale.
  for (s in 1:S) {
    X_star_zero[,s] ~ multi_normal(to_vector(mu0), Sigma0); 
    }
  
  // 5. Priors for the cutoffs.
  for (k in 1:K) {
    for(j in 1:cutpoint_count){
    cutpoints[k][j] ~ normal(cutpoint_prior_locations[j], 0.1); 
    // Empirical CDFs of data justify strong prior, assuming X_star is close to standard multinormal.
    }}
  
  // 8. Time of day prior:
  for(time in 1:nbeeps){
    time_of_day_intercept[time] ~ std_normal(); }
  
  // 6. Model target likelihood.
  for (k in 1:K){
    for (n in 1:N) {
      if(missing_mask[n,k] == 0){
       // Note, that X is N times K, whereas X_star is K times N.
       target += ordered_probit_lpmf( X[n,k] | X_star[k,n], cutpoints[k]);
        }}}

  }

generated quantities {

  // 1. Create the residual covariance matrix ( for interpretation)
  matrix[K, K] Omega;
  Omega = multiply_lower_tri_self_transpose(L_corr);
  
  // 2. Extract log densities
  array[N,K] log_lik;
  for (k in 1:K){
    for (n in 1:N) {
      if(missing_mask[n,k] == 0){
       // Note, that X is N times K, whereas X_star is K times N.
       log_lik[n,k] = ordered_probit_lpmf( X[n,k] | X_star[k,n], cutpoints[k]);
        } else {
      log_lik[n, k] = 0;
    }}}
  
  }
