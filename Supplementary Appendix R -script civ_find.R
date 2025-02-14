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
                     N = NULL) {

  K <- ncol(A)
  if(!random.init){
    # Initial values are taken based on A.
    Lambda <- Re(eigen(A)$vectors[,1])
    psi    <- Re(eigen(A)$values[1])
    omega0 <- omega1 <- runif(K, 0.5, 1)} else {
    Lambda <- rnorm(K,sd=sd(vec(A)))
    psi    <- rnorm(1)
    omega0 <- omega1 <- runif(K, 0.5, 1)
    omega0 <- omega0 + runif(K, 0.5, 1)
    }
  
  if(!cov.difference){
  
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
      "Optim_Result" = optim_result,
    ))}
  
  
  
  if(cov.difference) {
    # Population 2 times 2 block covariance matrix. (2K times 2K.)
    S                                  <- matrix(0,ncol=2*K,nrow=2*K)
    S[1:K                ,1:K]         <- var_ccov(A,Z,Delta=0)
    S[(K+1):(2*K)        ,(K+1):(2*K)] <- S[1:K,1:K]
    S[1:K                ,(K+1):(2*K)] <- var_ccov(A,Z,Delta=1)
    S[(K+1):(2*K)        ,1:K]         <- t(S[1:K,(K+1):(2*K)])
    
    if(any(eigen(S)$values  < 1e-7)) {
      message("The predicted covariance has near zero or negative eigenvals. Stopping.\n
              This might be a sign of low rank structure in the predicted covariance. \n
              Check eigen(var_ccov(A,Z,0))$values, for example, for low rank structure.")
      stop()
    }
    # Maximum likelihood estimate
    F_criterion <- function(theta, S, K){
      # 1) Parse the parameter vector:
      psi    <- theta[1]                   # scalar
      Lambda <- theta[2:(K+1)]             # length K
      omega0 <- theta[(K+2):(2*K+1)]       # length K
      omega1 <- theta[(2*K+2):(3*K+1)]     # length K
      
      # 2) Compute the implied covariance.
      S_implied                          <- matrix(0,ncol=2*K,nrow=2*K) 
      S_implied[1:K        ,1:K]         <- tcrossprod(Lambda)          + diag(omega0) #Within-Cov
      S_implied[(K+1):(2*K),(K+1):(2*K)] <- S_implied[1:K,1:K]          
      S_implied[1:K        ,(K+1):(2*K)] <- psi * tcrossprod(Lambda)    + diag(omega1) #Cross-Cov(1)
      S_implied[(K+1):(2*K),1:K]         <- t(S_implied[1:K,(K+1):(2*K)])

      # 3) Compute the discrepancy.
      return(F_ML=log(det(S_implied)) + sum(diag( solve(S_implied) %*% S )) - log(det(S)) - 2*K) #Note: 2*K, because we have 2K variables now.
    }
    
    # Maximum Likelihood statistic, with small perturbation to ensure invertiblity.
    lower_bounds  <- c(-Inf, rep(-Inf, times = K), rep(1e-6, times = 2*K))
    F_ML = optim(
      par = c(psi, Lambda, omega0, omega1 ), 
      fn  = F_criterion,
      S = S,
      K = K, 
      method = "L-BFGS-B", 
      lower = lower_bounds )
    
    # Return T statistic, which is distributed ~ Chisq_DF(F_ML - DF)
    statistic = (N-1) * F_ML$result$value

    # Compute degrees of freedom, RMSEA
    DF = ((2*K*(2*K+1))/2) - (3*K + 1)
    RMSEA = sqrt( max(0, statistic - DF) / (DF*(N-1)) )
    
    # Return results
    return(list(
      F_ML,
      RMSEA = RMSEA
    ))
    
  }

}



## Not run: test.

if(F){
  
  K=7
  A <- matrix(runif(K^2,-0.5,0.5), ncol = 7, nrow = 7)
  Z <- cov(t(A %*% t(matrix(rep(rnorm(1000), each = K), ncol = 7))) %*% t(A))
  civ_find(A,Z,cov.difference = T)
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
}



