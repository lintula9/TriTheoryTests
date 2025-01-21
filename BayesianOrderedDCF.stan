data {
  int<lower=1> S;               // Number of subjects
  int<lower=1> K;               // Number of variables (here, K = 3)
  int<lower=1> N;               // Total number of observations across all subjects
  int<lower=1> subject[N];      // Subject ID for each observation (sorted by subject and time)
  int<lower=1,upper=5> X[N,K];               // Observed data: time-series data for each time point
  int<lower=0,upper=1> missing_mask[N, K]; // Mask: 1 if X[i,k] is missing, 0 otherwise
  int<lower=1> start[S];        // start[s]: starting index for subject s
  int<lower=1> end[S];          // end[s]: ending index for subject s
  int<lower=1> cutpoint_count;
  real cutpoint_prior_locations[cutpoint_count];
}

parameters {
  // 1. D-CF model parameters
  real psi;                 // CF autoregression coefficient.
  vector<lower=0>[K] Lambda; // Factor loadings are assumed positive, to better identify the model.
  real<lower=0> subject_innovation_scale[S]; // Innovation variance is assumed to be time invariant.
  vector[N] eta_innovation;  // Innovations (for non-centralized parameterization).
  
  // 2. Initial state for eta, eta*:
  real eta_star[S];
  
  // 3. Random intercepts for each subject
  // The subject specific intercepts are defined for the CF.
  vector[S] subject_intercept_raw;  // subject-specific intercept deviations
  real<lower=0> subject_intercept_sd;     // Standard deviation for subject intercept.

  // 7. UNIQUE to ordered: X* (time and subject invariant) thresholds:
  ordered[cutpoint_count] cutpoints[K];
  }
  
transformed parameters {
  
  // 3. Improved posterior geometry with non-centralized priors.
  vector[S] subject_intercept;
  subject_intercept = subject_intercept_raw * subject_intercept_sd;
  
  // 99. c = 0 (written in for compatibility of all syntax, though redundant).
  real c;
  c = 0;
  
  // 1. Eta as a transformed parameter. Improved posterior geometry with non-centralized priors
  vector[N] eta;
  for (s in 1:S) {
    // Loop over time for each subject
    for (t in start[s]:end[s]) {
            if(t==start[s]){ 
          eta[t] = c + subject_intercept[s] + psi * eta_star[s] + eta_innovation[t]*subject_innovation_scale[s];
    } else {
          eta[t] = c + subject_intercept[s] + psi * eta[t-1] + eta_innovation[t]*subject_innovation_scale[s];
    }

        }
        } 
        
model {
  
  // 1. Priors for parameters
  psi ~    normal(0, 0.1);                 // Prior for DCF autoregression coefficient.
  Lambda ~ normal(0, 0.5);                // (half)-normal prior.
  eta_innovation ~ normal(0,0.1);     // Innovation prior.
  // 1a Subject specific innovation scale.
  subject_innovation_scale ~ normal(0,0.1);
  
  // 2. Initial state prior.
  to_vector(eta_star) ~ normal(0,0.5);

  // 3. Subject-specific intercepts
  subject_intercept_sd ~ normal(0,0.1); // Prior for subject intercept SDs
  to_vector(subject_intercept_raw) ~ normal(0, 0.1); // Raw intercept, scaled by SD.
  
  // 7. Declare priors for the cutoffs.
  for(k in 1:K){
    for(j in 1:cutpoint_count){
      cutpoints[k,j] ~ normal(cutpoint_prior_locations[j],0.2);
    }}
  
  
  // 4. Target likelihood loop over observations and variables.
  for (t in 1:N) {
    for (k in 1:K){
    // Skip missing values - though they are modeled in the latent eta.
      if(missing_mask[t,k] == 0){
        target += ordered_probit_lpmf( X[t,k] | Lambda[k]*eta[t], cutpoints[k]);
        }}}
        
    }
  
generated quantities {
  // We'll store the pointwise log-likelihood for each observation 
  // as a matrix of size N x K.
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
  // Possibly introduce a quick and dirty map to a VAR(1) model from D-CF(1)
  }
