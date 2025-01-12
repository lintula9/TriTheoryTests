# Supplementary Appendix R -script: Function for finding the closest stationary VAR(1) that is indistinguishable from a CF model.
# Running (sourcing) the whole script provides the function for the computation.

# Dependencies:
dependencies <- c("torch", "Matrix")
for( i in dependencies ) {
  
  if( !require(i, character.only = T) ) install.packages( i )
  library(i, character.only = T)
  }

# Function itself:
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


# Not run: Test
if(F) {
  
  A = matrix(c(0.5,0.2,0.3,
               0.1,0.1,0.1,
               0.2,0.3,0.2), ncol = 3)
  Z = matrix(c(0.2,0.1,0.0,
               0.1,0.2,0.0,
               0.0,0.0,0.2), ncol = 3)
  
  civ_find(A,Z)
  
  
}
