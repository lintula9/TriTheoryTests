# Supplementary Appendix R -script: Function for finding the closest stationary VAR(1) that is indistinguishable from a CF model.
# Running (sourcing) the whole script provides the function for the computation.

# Dependencies:
dependencies <- c("torch", "Matrix", "numDeriv")
for( i in dependencies ) {
  
  if( !require(i, character.only = T) ) install.packages( i )
  library(i, character.only = T)
}
# Helper:
var_ccov <- function(A,Z,Delta) {
  
  return( (A %^% Delta) %*% matrix(Matrix::solve(diag(1,ncol=ncol(A)^2,nrow=nrow(A)^2) - Matrix::kronecker(A,A)) %*% 
                                     fastmatrix::vec(Z), ncol = ncol(A), nrow = nrow(A)))
}

# Covariance squared difference
se_stat <- function(A,Z,A_,Z_,T_=10,nodiag = F) {
  selector = matrix(T,ncol(A),nrow(A))
  if(nodiag == T) selector <- lower.tri(A) + upper.tri(A)
  se = sum(sapply(  1:T_, 
    function(k) sum( ( var_ccov(A,Z,Delta=k)[selector] - var_ccov(A_,Z_,Delta=k)[selector] )^2 ))) +
      sum( ( tril(var_ccov(A,Z,Delta=0) - var_ccov(A_,Z_,Delta=0) )[selector ] )^2 )
  return(se)
  }

# Option 1, based on Torch (does not work?)
civ_find <- function( A, Z, n.iter = 1000L, tol = 1e-6 ) {
  
  K <- ncol(A)
  A <- torch_tensor(A)
  Z <- torch_tensor(Z)
  Lambda <- torch_randn(size = c(K, 1), requires_grad = TRUE)
  psi    <- torch_randn(size = c(1),     requires_grad = TRUE)
  
  loss_function <- function(Lambda = Lambda, psi = psi, A = A, Z = Z) {
    
    P_Lambda = Lambda$matmul(Lambda$transpose(1,2)) * 
      Lambda$square()$sum()$reciprocal()
    A_tilde = psi * P_Lambda
    Z_tilde = (1 - psi$square())*Lambda$matmul(Lambda$transpose(1,2))
    Objective_1 = A$matmul(P_Lambda) - psi * P_Lambda 
    Objective_2 = torch_tril(Z - (1-psi$square()) * Lambda$matmul(Lambda$transpose(1,2)))
    return( Objective_1$pow(2)$sum() + Objective_2$pow(2)$sum()  )}
  
  optimizer <- optim_adam(params = list(Lambda, psi), lr = 1e-2)
  
  num_iterations <- n.iter
  
 
  prev_loss <- Inf
  for (i in seq_len(num_iterations)) {
    
    optimizer$zero_grad()  
    loss <- loss_function(Lambda, psi, A, Z)
    loss$backward()
    optimizer$step()
    
    # Check for convergence
    loss_value <- loss$item()
    if (abs(prev_loss - loss_value) < tol) {
      cat(sprintf("Converged at iteration %d | loss = %.4f\n", i, loss_value))
      break
    }
    prev_loss <- loss_value
    
    if (i %% 50 == 0) {
      cat(sprintf("Iter %d | loss = %.4f\n", i, loss$item()))
    }
    
  }
  
  A_result = psi * Lambda$matmul(Lambda$transpose(1,2)) * Lambda$square()$sum()$reciprocal() 
  A_result = A_result +  (A - A_result)$matmul((torch_diag(rep(1, times = ncol(A))) - Lambda$matmul(Lambda$transpose(1,2)) * 
                                                  Lambda$square()$sum()$reciprocal())) 
  Z_result = (1 - psi$square())*Lambda$matmul(Lambda$transpose(1,2))
  
  gc(verbose = F)
  
  return(list("Loadings" = Lambda, "psi" = psi,
              "A_result" = A_result,
              "Z_result" = Z_result,
              "Optimizer" = optimizer)
         )
  
    }

# Option 2, simpler, seems faster and more reliable.

civ_find2 <- function(A, Z, n.iter = 2000, tol = 1e-6, W = NULL) {
  
  K <- ncol(A)
  # Initial values are taken based on A.
  Lambda <- Re(eigen(A)$vectors[,1])
  psi    <- Re(eigen(A)$values[1])
  
  loss_function <- function(pars, A = A, Z = Z) {
    # Loss written differently than in above:
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
    "A_result" = A_result,
    "B" = B_result,
    "C" = C_result,
    "Z_result" = Z_result,
    "Optim_Result" = optim_result
  ))
}



# Not run: Test
if(F) {
  
  A = matrix(c(0.5,0.2,0.3,
               0.1,0.1,0.1,
               0.2,0.3,0.2), ncol = 3)
  
  Z = matrix(c(1.0,0.6,0.0,
               0.6,1.0,0.0,
               0.0,0.0,0.2), ncol = 3)
  
  test_res2 <- civ_find2(A, Z)
  A - test_res2$A_result
  
  civ_find(A, Z)
  
  }

# Option 3, minimize the discrepancy in cross-covariance.
#Redunandt, see below.

civ_find3 <- function(A, Z, n.iter = 2000, tol = 1e-6, W = NULL, T_ = 10) {
  K <- ncol(A)
  # Initial values are taken based on A.
  Lambda <- Re(eigen(A)$vectors[,1])
  psi    <- Re(eigen(A)$values[1])
  
  loss_function <- function(pars, A = A, Z = Z) {
    # Loss based on cross-covariance.
    
    Z_ =  (1-pars[K+1]^2) * (pars[1:K] %*% t(pars[1:K])) 
    A_ = pars[K+1]*pars[1:K] %*% t(pars[1:K]) * sum(pars[1:K]^2)^-1
    
    M = sum(sapply(  1:T_, 
                 function(k) sum( ( var_ccov(A,Z,Delta=k) - var_ccov(A_,Z_,Delta=k) )^2 ))) +
      sum( (fastmatrix::vech( var_ccov(A,Z,Delta=0) - var_ccov(A_,Z_,Delta=0) ))^2 )
    return( M )
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
    "A_result" = A_result,
    "B" = B_result,
    "C" = C_result,
    "Z_result" = Z_result,
    "Optim_Result" = optim_result
  ))
}

# Option 4, minimize the discrepancy in cross-covariance disregarding diagonal cross-covariance.
# I.e., we can assume there being idiosyncratic serial correlatedness to the symptoms, not captured by the CF.
# Thus the diagonal are always a perfect match - as long as the predicted diagonal cross-covariance is less than what the
# original VAR predicts.

civ_find4 <- function(A, Z, n.iter = 2000, tol = 1e-6, W = NULL, T_ = 10, nodiag = T) {
  K <- ncol(A)
  # Initial values are taken based on A.
  Lambda <- Re(eigen(A)$vectors[,1])
  psi    <- Re(eigen(A)$values[1])
  
  loss_function <- function(pars, A = A, Z = Z) {
    # Loss based on cross-covariance. Omitting the diagonals of cross-covariances (if the discrepancy is positive).
    
    Z_ =  (1-pars[K+1]^2) * (pars[1:K] %*% t(pars[1:K])) 
    A_ = pars[K+1]*pars[1:K] %*% t(pars[1:K]) * sum(pars[1:K]^2)^-1
    
    # Sum the non-diagonal entries
    nodiag <- lower.tri(A) + upper.tri(A)
    M = se_stat(A = A,Z = Z,A_ = A_,Z_ = Z_,T_ = T_, nodiag = nodiag)
    return( M )
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
    "A_result" = A_result,
    "B" = B_result,
    "C" = C_result,
    "Z_result" = Z_result,
    "Optim_Result" = optim_result
  ))
}

if(F){
  
  civ_find3(A, Z) # Result is likely to be pretty close to civ_find2.
  testres_temp <- civ_find4(A, Z) # Result is likely very different.
  
  se_stat(A,Z,testres_temp$A_result,testres_temp$Z_result,nodiag = T)
  
}











