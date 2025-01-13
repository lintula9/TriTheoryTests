# Supplementary Appendix R -script: Function for finding the closest stationary VAR(1) that is indistinguishable from a CF model.
# Running (sourcing) the whole script provides the function for the computation.

# Dependencies:
dependencies <- c("torch", "Matrix", "numDeriv")
for( i in dependencies ) {
  
  if( !require(i, character.only = T) ) install.packages( i )
  library(i, character.only = T)
  }

# Option 1, based on Torch (does not work?)
civ_find <- function( A, Z, n.iter = 1000L ) {
  
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
    M = A$matmul(P_Lambda) - psi * P_Lambda + Z - (1-psi$square()) * Lambda$matmul(Lambda$transpose(1,2))
    return( M$pow(2)$sum() )}
  
  optimizer <- optim_adam(params = list(Lambda, psi), lr = 1e-2)
  
  num_iterations <- n.iter
  
  for (i in seq_len(num_iterations)) {
    
    optimizer$zero_grad()  
    loss <- loss_function(Lambda, psi, A, Z)
    loss$backward()
    optimizer$step()
    
    if (i %% 50 == 0) {
      cat(sprintf("Iter %d | loss = %.4f\n", i, loss$item()))
    }
    
  }
  
  A_result = psi * Lambda$matmul(Lambda$transpose(1,2)) * 
    Lambda$square()$sum()$reciprocal()
  Z_result = (1 - psi$square())*Lambda$matmul(Lambda$transpose(1,2))
  
  return(list("Loadings" = Lambda, "psi" = psi,
              "A_result" = A_result,
              "Z_result" = Z_result) )
  
    }





# Option 2, simpler.

civ_find2 <- function(A, Z, n.iter = 500) {
  
  K <- ncol(A)
  # Initial values are taken based on A.
  Lambda <- eigen(A)$vectors[,1]
  psi    <- eigen(A)$values[1]
  
  loss_function <- function(pars, A = A, Z = Z) {
    
    P_Lambda = pars[1:K] %*% t(pars[1:K]) * (sum(pars[1:K]^2))^-1
    M = A %*% P_Lambda - pars[K+1] * P_Lambda + Z - (1-pars[K+1]^2) * pars[1:K] %*% t(pars[1:K])
    
    return( sum(M^2) )
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
        A = A, Z = Z) 
  
  optimized_params <- optim_result$par
  Lambda_opt <- optimized_params[1:K]
  psi_opt <- optimized_params[K + 1]
  
  A_result <- psi_opt * Lambda_opt %*% t(Lambda_opt) * (sum(Lambda_opt^2))^-1
  Z_result <- (1 - psi_opt^2) * Lambda_opt %*% t(Lambda_opt)
  
  # Return results
  return(list(
    "Loadings" = Lambda_opt, 
    "psi" = psi_opt,
    "A_result" = A_result,
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

# Evaluate the amount of misfit somehow...