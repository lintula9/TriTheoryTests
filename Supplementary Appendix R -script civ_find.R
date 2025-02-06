# Dependencies:
if(!require(fastmatrix)) install.packages("fastmatrix"); library(fastmatrix)

civ_find <- function(A, Z, n.iter = 2000, tol = 1e-6, 
                     W = NULL,
                     random.init = F) {

  K <- ncol(A)
  if(!random.init){
    # Initial values are taken based on A.
    Lambda <- Re(eigen(A)$vectors[,1])
    psi    <- Re(eigen(A)$values[1])} else {
    Lambda <- rnorm(K,sd=sd(vec(A)))
    psi <- rnorm(1)
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

  
  # Return results
  return(list(
    "Loadings" = Lambda_opt, 
    "psi" = psi_opt,
    "A" = A_result,
    "B" = B_result,
    "C" = C_result,
    "Z" = Z_result,
    "Optim_Result" = optim_result
  ))
}

