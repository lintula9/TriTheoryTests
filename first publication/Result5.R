# Result 5.
library(fastmatrix); library(Matrix); library(MASS)

  # Stationary 2 dimensional D-CF(1) model onto VAR(1).

  # Define the parameters.

L = matrix(runif(6,min = 0.2, max = 1), nrow = 3) # Time invariant CF loadings randomly from uniform ( 0.2, 1).

A_eta = matrix(c(0.5,0.2,0.2,0.5), ncol = 2) # CF autoregression.
A_Z = matrix(c(0.2 , 0.0 ,0.0 , 0.2), ncol = 2) # CF innovation covariance
Sigma_eta = matrix( solve(diag(1, ncol = 2*2, nrow = 2*2) - kronecker.prod(A_eta) ) %*% vec(A_Z),
        ncol = 2, nrow = 2)

  # Apply the map.
A = L %*% A_eta %*% ginv(L) # Remember to ensure we have produced a stationary VAR(1).
Z = L %*% Sigma_eta %*% t(L) - L %*% A_eta %*% Sigma_eta %*% t(A_eta) %*% t(L) 

  # Compute the (cross-)covariances.

Sigma_VAR = matrix( solve(diag(1, 3*3, 3*3) - kronecker.prod(A)) %*% vec(Z),
                    ncol = 3, nrow = 3)

Sigma_DCF = matrix(matrix( L %*% Sigma_eta %*% t(L) ),
                   ncol = 3, nrow = 3)

  # Check compatibility.

round(Sigma_VAR, 10) == round(Sigma_DCF, 10)

  # Compute cross-covariance up to Deltat = 2.

K_1_VAR = A %*% Sigma_VAR
K_2_VAR = A %*% A %*% Sigma_VAR
K_1_DCF = L %*% ( A_eta %*% Sigma_eta ) %*% t(L)
K_2_DCF = L %*% ( A_eta %*% A_eta %*% Sigma_eta ) %*% t(L)

round(K_1_VAR, 10) == round(K_1_DCF, 10)
round(K_2_VAR, 10) == round(K_2_DCF, 10)

  #... and so forth.

  # Add in B(t), orthogonal to L.

B_0 = rWishart(1, df = 4, Sigma = diag(0.2, nrow = 3))[,,1] # Initiate a matrix B_0 randomly.
B = B_0 %*% ( diag(1, nrow = 3, ncol = 3) - L %*% solve(t(L) %*% L) %*% t(L) )

round(B %*% L, 10)

A_adj = A + B

Sigma_VAR_adj = matrix( solve(diag(1, nrow = 3*3, ncol = 3*3) - kronecker.prod(A_adj)) %*% vec(Z),
                        ncol = 3, nrow = 3)
K_1_VAR_adj = A_adj %*% Sigma_VAR_adj
K_2_VAR_adj = A_adj %*% A_adj %*% Sigma_VAR_adj

round(K_1_VAR_adj, 10) == round(K_1_VAR, 10)
round(K_2_VAR_adj, 10) == round(K_2_VAR, 10)

  #... and so forth.


