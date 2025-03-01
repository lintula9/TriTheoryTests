
# Helper functions -------------
# Dependencies:
if(!require(fastmatrix)) install.packages("fastmatrix"); library(fastmatrix)
if(!require(expm)) install.packages("expm"); library(expm)

# Helper function, for the (cross-)covariance
var_ccov <- function(coefmat,innocov,Delta) {

  if(any(Re(eigen(coefmat)$values) > 1)) {print(simpleError("A is not stationary. Aborting."))
    stop()}
  
  # Yule-Walker equation solution for stationary A.
  return( covmat = (coefmat %^% Delta) %*% 
            matrix(Matrix::solve(diag(1,ncol=ncol(coefmat)^2,nrow=nrow(coefmat)^2) - 
                                   Matrix::kronecker(coefmat,coefmat)) %*% 
                     fastmatrix::vec(innocov), 
                   ncol = ncol(coefmat), nrow = nrow(coefmat)))
}

# Helper function to form M times M block matrix, of the (cross-)covariances.
var_ccov_Mblock <- function(A, Z, M) {
  K <- nrow(A)      # dimension of each block (A is K x K)
  big_dim <- M * K  # dimension of final matrix
  S_big <- matrix(0, nrow = big_dim, ncol = big_dim)
  
  for (i in seq_len(M)) {
    for (j in seq_len(M)) {
      
      delta <- abs(i - j)
      
      block_ij <- var_ccov(A, Z, Delta = delta)
      
      if (i > j) {
        block_ij <- t(block_ij)
      }
      
      row_idx <- ((i-1)*K + 1) : (i*K)
      col_idx <- ((j-1)*K + 1) : (j*K)
      
      S_big[row_idx, col_idx] <- block_ij
    }
  }
  
  return(S_big)
}

                # Main function --------------

# Closest VAR(1) model - indistinguishable from a dynamic common factor model. 

# Run the whole script, then civ_find takes the coefficient matrix,
# innovation covariance as arguments. If you want to compute RMSEA, you also
# need to give sample size N and specify cov.difference = T. In this case, we need
# diagonal errors to be allowed - otherwise the maximum likelihood criterion is
# not solvable (because the implied covariance will be (in practice nearly) low rank).


civ_find <- function(A, Z, 
                     n.iter = 2000, tol = 1e-6, 
                     W = NULL,
                     random.init    = T,
                     cov.difference = T,
                     time_points = 2,
                     N = NULL,
                     error_ratio = 0.5) {
  
  if(time_points < 2) simpleError("time_points must be 2 or larger.")

  K <- ncol(A)
  if(!random.init){
    # Initial values are taken based on A.
    Lambda <- Re(eigen(A)$vectors[,1])
    psi    <- Re(eigen(A)$values[1])} else {
    Lambda <- rnorm(K,sd=sd(vec(A)))
    psi    <- min(rexp(1), 0.9)    }
  
  omega0 <- rexp(K) 

  if(!cov.difference){
    
  message("Minimizing distance in terms of parameters (not covariance), as cov.difference was set to False.")
    
    loss_function <- function(pars, A = A, Z = Z) {
      z <- vech(Z)
      s <- c(vec(A),z)
      tilde_z = vech( (1-pars[K+1]^2) * (pars[1:K] %*% t(pars[1:K])) )
      C = pars[K+1]*pars[1:K] %*% t(pars[1:K]) * sum(pars[1:K]^2)^-1
      B = (A - C) %*% (diag(1, nrow = K, ncol = K) - pars[1:K] %*% t(pars[1:K]) * sum(pars[1:K]^2)^-1)
      tilde_s = c(vec( C + B ), tilde_z)
      if(is.null(W)) W = diag(1, nrow = length(s), ncol = length(s))
      M = t(s-tilde_s) %*% W %*% (s-tilde_s)
      return( as.numeric(M) )
    }
    
    gradient_function <- function(pars, A = A, Z = Z) {
      numDeriv::grad(func = loss_function, x = pars, A = A, Z = Z)
    }
    
    params = c(Lambda,psi)
    names(params) <- c(paste0("Lambda", 1:K), "psi")
    optim_result <- optim(par = params, 
                          fn = loss_function, 
                          gr = gradient_function, 
                          method = "BFGS",
                          control = list(abstol = tol, maxit = n.iter),
                          A = A, Z = Z) 
    
    optimized_params <- optim_result$par
    Lambda_opt <- optimized_params[1:K]
    psi_opt    <- optimized_params[K + 1]
    
    C_result <- psi_opt * Lambda_opt %*% t(Lambda_opt) * (sum(Lambda_opt^2))^-1 
    B_result <- (A - C_result) %*% (diag(1,ncol = ncol(A),nrow = nrow(A)) - 
                                      Lambda_opt %*% t(Lambda_opt) * sum(Lambda_opt^2)^-1)
    A_result <- C_result + B_result
    Z_result <- (1 - psi_opt^2) * Lambda_opt %*% t(Lambda_opt)
    
    gc(verbose = F)
    
    # Return results
    return(list(
      "Loadings" = Lambda_opt, 
      "psi"      = psi_opt,
      "A"        = A_result,
      "Z"        = Z_result,
      "Optim_Result" = optim_result
    ))}
  
  
  
  if(cov.difference) {
    # Population 2 times 2 block covariance matrix. (2K times 2K.)
    S  <- var_ccov_Mblock(A,Z,time_points)
    
    if(any(eigen(S)$values  < 1e-7)) {
      message("The VAR(1) predicted covariance has near zero or negative eigenvals. Stopping.\n
              This might be a sign of low rank structure in the predicted covariance. \n
              Check eigen(var_ccov(A,Z,0))$values, for example, for low rank structure.")
      stop()
    }
    
    diag(S) <- diag(S) + (error_ratio * diag(S))
    
    # Maximum likelihood criterion
    F_criterion <- function(theta, S, K, return_cov = F){
      # 1) Parse the parameter vector:
      psi    <- theta[1]                   # scalar
      Lambda <- theta[2:(K+1)]             # length K
      omega0 <- theta[(K+2):(2*K+1)]       # length K
      
      A_indistinguishable = psi * diag(1,ncol(A),nrow(A))        # A_indistinguishable
      Z_indistinguishable = (1-psi^2) * tcrossprod(Lambda)

      # 2) Compute the implied covariance.
      S_implied <- var_ccov_Mblock(A_indistinguishable, Z_indistinguishable, time_points) + diag(rep(omega0, times = time_points), time_points*K, time_points*K)
      
      if(return_cov) {return(S_implied)}
      
      if(any(eigen(S_implied)$values < 1e-6)) {warning("Non positive-definite implied covariance.")}

      # 3) Compute the discrepancy.
      return(F_ML = log(det(S_implied)) + 
                    sum(diag( solve(S_implied) %*% S )) - 
                    log(det(S)) - 
                    2*K) 
    # Note: 2*K, because we have 2K variables.
    }
    
    # Maximum Likelihood statistic, with constraints:
    # Error variances and serial covariances positive.
    # The CF is implicitly assumed to have unit variance.
    # DCF is constrained stationary.
    lower_bounds  <- c(-1, rep(-Inf, times = K), rep(0,   times = K))
    upper_bounds  <- c( 1, rep( Inf, times = K), rep(Inf, times = K))
    F_ML  = optim(
      par = c(psi, Lambda, omega0 ), 
      fn  = F_criterion,
      S   = S,
      K   = K, 
      method = "L-BFGS-B", 
      lower  = lower_bounds,
      upper  = upper_bounds,
      control = list(maxit = n.iter))
    
    # Return T statistic, which is distributed ~ Chisq_DF(F_ML - DF)
    statistic = (N-1) * F_ML$value

    # Compute degree of freedom
    DF    = length(fastmatrix::vech(Z)) + length(A) - length(F_ML$par)
    
    # Compute RMSEA
    RMSEA = sqrt( max(0, statistic - DF) / (DF*(N-1)) )
    
    # Return implied covariance
    S_implied <- F_criterion(theta = F_ML$par, S = S, K = K, return_cov = T)
    
    # Warnings.
    if(F_ML$value < 0) warning("Negative discrepancy observed. Likely the implied covariance was non-invertible.
                               \n Results are not interpretable")
    
    # Return results
    return(list(
      "ML_optimizer_result"           = F_ML,
      T_statistic                     = statistic,
      RMSEA                           = RMSEA,
      CF_ar_coef                      = F_ML$par[1],
      Factor_loadings                 = F_ML$par[2:(K+1)],
      within_time_point_errorvars     = F_ML$par[(K+2):(2*K+1)],
      DF                              = DF,
      "Implied_covariance"            = S_implied,
      "VAR_original_covariance"       = S
    ))
    
  }

}



## Not run: test.

if(F){
# Additional function to check how much additional time points increase RMSEA.

RMSEA_check <- function(A,Z,N,max_time_points){
  RMSEAs <- c()
  for( i in 2:max_time_points){
    RMSEAs[i] <- civ_find(A, Z, N, cov.difference = T, 
                          random.init = T, time_points = i)$RMSEA}
  return(RMSEAs)
}}



# Numerical examples -------------------------------

if(F){
  
  A_2 <- matrix(c(
    0.5, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.4, 0.1, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.4, 0.1, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.5, 0.1, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.3, 0.1, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.4, 0.1,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.2
  ), nrow = 7, byrow = TRUE)
  
  
  # Identity on the diagonal; 0.2 for first off-diagonals
  Z_2 <- diag(7)
  Z_2[cbind(1:6, 2:7)] <- 0.2
  Z_2[cbind(2:7, 1:6)] <- 0.2
  
  # The fit is bad.
  result <- civ_find(A_2, Z_2, N = 4200, cov.difference = T, random.init = T, time_points = 4)

  # Proportion explained is not high, though factor congruency (for the first factor) remains decently large.
  civ_parallel(A_2, Z_2)
  
  }




# CIV parallel --------------

civ_parallel <- function(A,Z,time_points = 10) {
  
  if(any(abs(eigen(A)$values) > 1)) simpleError("Non-stationary A, aborting.")
  
  comp        <- lapply(0:time_points,         function(d){return(eigen(var_ccov(A,Z,d)))})
  all_sum     <- sapply(1:ncol(A), function(j){
                        sum(sapply(1:(time_points+1), function(i) comp[[i]]$values[j]))})
  total_sum   <- sum(sapply(1:(time_points+1), function(i) sum(comp[[i]]$values)))
  cosine_i    <- sapply(1:ncol(A), function(j) {
    sapply(1:time_points, function(i) abs(sum(comp[[i]]$vectors[,j]*comp[[i+1]]$vectors[,j])) )}, simplify = "matrix")
  
  prop_explained             = all_sum / total_sum
  if(any(Im(prop_explained) != 0)) warning("Imaginary proportion explained found, 
  likely due to the covariances being (close) to low rank. Use lower time_points.")

  return(
    list(
      prop_explained                = prop_explained,
      min_factor_congruency         = apply(cosine_i, min, MARGIN = 2),
      all_factor_congruencies       = cosine_i
      )
    
  ) }
