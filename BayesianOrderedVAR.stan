data {
  // General integers.
  int<lower=1> S;               // Number of subjects
  int<lower=1> K;               // Number of variables (here, K = 3)
  int<lower=1> N;               // Total number of observations across all subjects
  int<lower=1> subject[N];      // Subject ID for each observation (sorted by subject and time)
  int<lower=1> N_obs;        // Number of observations
  int<lower=1> obsval_slice_start[K]; // Each variabe gets its own slice indicator, which tells where in the X_observed the variable starts.
  int<lower=1> obsval_slice_end[K];   // And ends
  
  // Matrix form data: Commented out 26.1.
  // int<lower=1,upper=5> X[N,K];   // Observed data: matrix form. (Missing values contained!)
  // int<lower=0,upper=1> missing_mask[N, K]; // Mask: 1 if X[i,k] is missing, 0 otherwise. Matrix form.
  
  // Vector form data:
  // We can ignore this. int<lower=1> X_vector[N*K];           // Observed data, in vector form. (Missing values contained!)
  // And this. int<lower=1, upper=K> X_k[N*K];       // Observed data variable indicator for vector form.
  int<lower=1> observed_indices[N_obs]; // Observed data indicator. (Leave out missing data in target likelihood.)
  // The above indices are used to obtain only the data points of the X* parameters, which correspond to an observed variable.
  int<lower=1> X_observed[N_obs];
  
  // Subject indexing.
  int<lower=1> start[S];        // start[s]: starting index for subject s
  int<lower=1> end[S];          // end[s]: ending index for subject s
  
  // Cutpoints for probit measurement model.
  int<lower=1> cutpoint_count; // Note that we assume (know) that category count is the same for all vars.
  real cutpoint_prior_locations[cutpoint_count];
  
  // Initial, t = 0, unobserved state.
  vector[K] mu0;              // Mean for initial states
  cov_matrix[K] Sigma0;       // Covariance for initial states
  int<lower=1> beep[N];         // Time of day, for each observation. (4 = evening, 1 = morning.)
  int<lower=1> nbeeps;          // Count of different times of day.
}

parameters {
  // 1. VAR(1) model parameters
  matrix[K, K] A;           // VAR(1) coefficient matrix
  matrix[K, K] B;        // Evening to morning transition adjustment matrix.
  
  // Note: the scale of the latent X* symtpoms is not identified, thus it is user defined.
  
  // 2. Correlation matrix for innovations.
  cholesky_factor_corr[K] L_corr;  // Cholesky factor of the correlation matrix

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
            if(t==start[s]){    
          X_star[,t] = c + subject_intercept[s] 
                         + A * X_star_zero[,s] 
                         + time_of_day_intercept[beep[t]];
    } else {
          X_star[,t] = c + subject_intercept[s] 
                         + A * X_star[,t-1] 
                         + L * ( X_star_innovation[,t] ) 
                         + time_of_day_intercept[beep[t]]
                         + (beep[t] == 4) * B * X_star[,t-1];
    }}}
    
   // Vectorized parameterization:
    vector[N_obs] X_star_observed_vector;
    X_star_observed_vector = to_vector(X_star)[observed_indices]; 
    // Unrolls column-wise, only take those that correspond to observed variables.
}

model {
  
  // 1. Priors for VAR model parameters
  to_vector(A) ~ normal(0, 0.5);               // Prior for VAR coefficients
  to_vector(B) ~ normal(0, 0.5);               // Prior for evening-to-morning adjustment

  L_corr ~ lkj_corr_cholesky(1);              // Prior for correlation matrix
  to_vector(X_star_innovation) ~ normal(0,0.5); 
  // After mixing with L: X -> LX, we obtain multinormal variables.

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

  // Matrix form target likelihood specification: (commented out 26.1.)
  //for (k in 1:K){
  //  for (n in 1:N) {
  //    if(missing_mask[n,k] == 0){
  //     // Note, that X is N times K, whereas X_star is K times N.
  //     target += ordered_probit_lpmf( X[n,k] | X_star[k,n], cutpoints[k]);
  //      }}}
  
  // Vector form target likelihood: (This can be faster as it does not require missing_mask)
  
  // LOOP OVER THE SHORTER 1:K!!!
   for (k in 1:K) {
    target += ordered_probit_lpmf(X_observed[obsval_slice_start[k]:obsval_slice_end[k]] | X_star_observed_vector[obsval_slice_start[k]:obsval_slice_end[k]], cutpoints[k]);
    }

  }

generated quantities {

  // 1. Create the residual covariance matrix ( for interpretation)
  matrix[K, K] Omega;
  Omega = multiply_lower_tri_self_transpose(L_corr);
  
}
