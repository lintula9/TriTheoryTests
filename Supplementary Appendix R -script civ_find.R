# Run the whole script, then civ_find takes the coefficient matrix,
# innovation covariance as arguments. If you want to compute RMSEA, you also
# need to give sample size N, and error_ratio. Error_ratio is by default 0.5
# indicating that half of an observation is error (e.g., 0.7 factor loading relates
# to 0.49 error proportion in a one factor model).


# Dependencies:
if(!require(fastmatrix)) install.packages("fastmatrix"); library(fastmatrix)
if(!require(expm)) install.packages("expm"); library(expm)


# Helper function, cor the (cross-)covariance
var_ccov <- function(coefmat,innocov,Delta) {
  
  if(any(Re(eigen(coefmat)$values) > 1)) print(simpleError("A is not stationary. Aborting."))
  
  return( covmat = (coefmat %^% Delta) %*% 
            matrix(Matrix::solve(diag(1,ncol=ncol(coefmat)^2,nrow=nrow(coefmat)^2) - 
                                                            Matrix::kronecker(coefmat,coefmat)) %*% 
                                                  fastmatrix::vec(innocov), 
                                            ncol = ncol(coefmat), nrow = nrow(coefmat)))
   }

civ_find <- function(A, Z, n.iter = 2000, tol = 1e-6, 
                     W = NULL,
                     random.init = F,
                     cov.difference = F,
                     N = NULL,
                     error_ratio = 0.5) {

  K <- ncol(A)
  if(!random.init){
    # Initial values are taken based on A.
    Lambda <- Re(eigen(A)$vectors[,1])
    psi    <- Re(eigen(A)$values[1])} else {
    Lambda <- rnorm(K,sd=sd(vec(A)))
    psi    <- rnorm(1)
    }
  
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
  psi_opt <- optimized_params[K + 1]
  
  C_result <- psi_opt * Lambda_opt %*% t(Lambda_opt) * (sum(Lambda_opt^2))^-1 
  B_result <- (A - C_result) %*% (diag(1,ncol = ncol(A),nrow = nrow(A)) - 
                                    Lambda_opt %*% t(Lambda_opt) * sum(Lambda_opt^2)^-1)
  A_result <- C_result + B_result
  Z_result <- (1 - psi_opt^2) * Lambda_opt %*% t(Lambda_opt)
  
  gc(verbose = F)

  if(cov.difference) {
    # Population 2 times 2 block covariance matrix. (2K times 2K.)
    S                                  <- matrix(0,ncol=2*K,nrow=2*K)
    S[1:K                ,1:K]         <- var_ccov(A,Z,Delta=0)
    S[(K+1):(2*K)        ,(K+1):(2*K)] <- S[1:K,1:K]
    S[1:K                ,(K+1):(2*K)] <- var_ccov(A,Z,Delta=1)
    S[(K+1):(2*K)        ,1:K]         <- t(S[1:K,(K+1):(2*K)])
    
    # CF implied 2 times 2 -||-.
    S_implied                          <- matrix(0,ncol=2*K,nrow=2*K)
    S_implied[1:K        ,1:K]         <- var_ccov(A_result,Z_result,Delta=0)
    S_implied[(K+1):(2*K),(K+1):(2*K)] <- S_implied[1:K,1:K]
    S_implied[1:K        ,(K+1):(2*K)] <- var_ccov(A_result,Z_result,Delta=1)
    S_implied[(K+1):(2*K),1:K]         <- t(S_implied[1:K,(K+1):(2*K)])  
    
    # Ensure positive definitiveness.
    error = diag(S) * (if(exists("error_ratio")) error_ratio else 1e-6)
    S_implied <- S_implied + diag(max(1e-6,error), ncol = 2*K, nrow = 2*K)
    S         <- S         + diag(max(1e-6,error), ncol = 2*K, nrow = 2*K)
    
    # Maximum Likelihood statistic, with small perturbation to ensure invertiblity.
    F_ML = log(det(S_implied)) + sum(diag( solve(S_implied) %*% S )) - log(det(S)) - 2*K
    statistic = (N-1) * F_ML
    
    # Compute RMSEA
    DF = ((2*K*(2*K+1))/2) - (K + 1)
    RMSEA = sqrt( max(0, statistic - DF) / (DF*(N-1)) )
    
  }
  
  # Return results
  return(list(
      "Loadings" = Lambda_opt, 
      "psi"      = psi_opt,
      "A"        = A_result,
      "Z"        = Z_result,
    "Optim_Result" = optim_result,
    # If the cov.difference is computed.
        "RMSEA"  = if(exists("RMSEA")) RMSEA else NULL,
        "Chisq.statistic"  = if(exists("statistic")) statistic else NULL,
  ))
  }

