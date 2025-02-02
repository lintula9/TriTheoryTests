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
  array[S] int<lower=1> start;        // start[s]: starting index for subject s
  array[S] int<lower=1> end;          // end[s]: ending index for subject s
  
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
  matrix[K, K] A;           // VAR(1) coefficient matrix.
  real psi; // The CF autoregression coefficient.

  // 2. Declare Lambda (factor loadings) explicitly. 
  // COMMENT: constrain to 1 to identify scale: real<lower=0> Lambda_first;             // First element, constrained to be positive for sign identification.
  vector[K-1] Lambda_rest;                // Remaining elements with no sign constraint
  array[K-1] real<lower=0> L_rest_diag;       // Must be non-negative, to ensure PSD of covariance.
  vector[(K*(K-1)%/%2)] L_rest_offdiag; // Unrestricted.

  // 3. Random intercepts for each subject
  matrix[S, K] subject_intercept_raw;  // subject-specific intercept deviations
  array[S] real<lower=0> subject_intercept_sd;     // Standard deviation for subject intercepts

  // 4. X_star as latent symptoms.
  // Note: the scale of the latent X* symtpoms is not identified, thus it is user defined.
  matrix[K,N] X_star_innovation; // The non-deterministic part of X_star.
  matrix[K,S] X_star_zero;   // Initial state for each subject at t=0.

  // 5. Probit measurement model cutoffs:
  array[K] ordered[cutpoint_count] cutpoints;
  
  // 8. Time of day parameter: We set the first time of day effect as the reference
  vector[nbeeps] ref_time_of_day_effect;
  array[nbeeps] vector[K-1] specific_time_of_day_effect; // and the other are relative to that.
  // In this way, we can model their discrepancies as well (or deviation from the reference, at least).
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
  
  // 1. Create the Cholesky factor.
  vector[K] Lambda; // Declare Lambda as the whole vector.
  // COMMENT: constrain to one to identify scale: Lambda[1] = Lambda_first;
  Lambda[1] = 1;
  Lambda[2:K] = Lambda_rest;
  
  matrix[K,K] L_;
  L_ = rep_matrix(0.0, K, K); // L_ will be the cholesky factor (of the covariance).
  L_[,1] = Lambda;
  for( i in 2:K) L_[i,i] = L_rest_diag[i-1];
  for( i in 2:(K)) {
    int counter = 2 + (i - 2) * K - ((i - 1) * i) %/% 2;
    L_[(i+1):K, i] = L_rest_offdiag[counter:(counter + (K - i) - 1)];
  }
      // Compute L_corr
  matrix[K,K] L_corr;
  L_corr = diag_matrix(1.0 ./ sqrt(diagonal(L_*L_'))) * L_;
  
  // 1. Create A.
  matrix[K,K] A_effective;
  A_effective = add_diag(A, psi);
  
  // 8. Create time of day effect. The reference is set to be the first.
  array[nbeeps] vector[K] time_of_day_effect;
  for(i in 1:nbeeps){
    time_of_day_effect[i][1] = ref_time_of_day_effect[i];
    for(k in 2:(K))
    time_of_day_effect[i][k] = ref_time_of_day_effect[i] + specific_time_of_day_effect[i][k-1];
  }
  
  // 4. Latent symptoms X_star as a transformed parameter: 
  matrix[K,N] X_star;
  for (s in 1:S) {
    // Loop over time for each subject
    for (t in start[s]:end[s]) {
            if(t==start[s]){    
          X_star[,t] = c + A_effective * X_star_zero[,s] 
                         + time_of_day_effect[beep[t]];
    } else {
          X_star[,t] = c + subject_intercept[,s] 
                         + A_effective * X_star[,t-1] 
                         + L_corr * X_star_innovation[,t]  // Transforms our standard normal (white noise) innovations to correlated innovations.
                         + time_of_day_effect[beep[t]];
    }}}
    
}

model {
  
  // 1. Priors for VAR model parameters
  // Brief note:
  // A diagonal A, with structure A = psi * I, satisfies indistinguishability conditions.
  // Thus, _any_ VAR(1) satisfying indistinguishability conditions is indistinguishable from
  // psi*I, Lambda * Lambda' * constant - type of parameterization. Also note, in our case we assume mvnormality, 
  // thus the condition is stronger.
  // Using this knowledge, we present the prior as follows:
  psi ~ normal(0,0.5);              // CF autoregression prior.
  to_vector(A) ~ normal(0,0.01);    // As per our prior belief, we set a strict prior for A parameters, outside of the diagonal psi.
  // To finalize, the computation: A_effective = A + psi * I is conducted in transformed parameters.
  // Effectively, using non-centralized parameterizatoin, psi will denote the location of the diagonal, and deviations from
  // psi * I are punished by our prior belief heavily. Strong prior should also identify our model.
  
  // Prior for innovation covariance.
  // COMMENT: Set to 1 for identification. Lambda_first ~ normal(0.5,0.5);         // Prior, assuming that the factor loadings are non-negligible.        
  Lambda_rest ~ normal(0.5,0.5);          
  L_rest_diag ~ normal(0,0.01);           // Tight prior - our prior is that the true model is CF model.
  L_rest_offdiag ~ normal(0,0.01); 
  
  // Prior for the white noise innovatoins, afterwards transformed to correlated innovations.
  to_vector(X_star_innovation) ~ std_normal();

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
  
  // 8. Time of day parameter:
  // Note, that the effect of time of day should be equivalent, up to lambda. I.e., the effect is informally
  // time_of_day -> CF -> Lambda_i * X_i. This will be reflected on our priors as well.
  ref_time_of_day_effect ~ std_normal();
  for(i in 1:nbeeps){
  specific_time_of_day_effect[i] ~ normal(0,0.01); // Strong prior, as our belief is that they are all the same.
  }
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
  array[N,K] real log_lik;
  for (k in 1:K){
    for (n in 1:N) {
      if(missing_mask[n,k] == 0){
       // Note, that X is N times K, whereas X_star is K times N.
       log_lik[n,k] = ordered_probit_lpmf( X[n,k] | X_star[k,n], cutpoints[k]);
        } else {
      log_lik[n, k] = 0;
    }}}
  
  }

