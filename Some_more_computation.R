# Install packages if not already installed
install.packages("lavaan")
install.packages("vars")
install.packages("MASS")
install.packages("expm")

# Load the packages
library(lavaan)
library(vars)
library(MASS)
library(fastmatrix)
library(expm)


set.seed(123)

library(MASS)

library(MASS)

# Define the Lambda vector (4x1)
Lambda <- matrix(c(1, 2, 3, 4), ncol = 1)

# Define the autoregressive coefficient delta
delta <- 0.5

# Define Z as a rank-1 matrix with the 1/(1 - delta^2) factor
Z <- (1 - delta^2) * Lambda %*% t(Lambda)

# Compute sigma
Sigma = Lambda %*% t(Lambda)

# Compute A (in any way you want)
A = (delta / (1 - delta^2)) * Z %*% ginv(Sigma)
A = (Lambda %*% t(Lambda) * delta) %*% ginv(Sigma)
eigen(A)$vectors[,1] * (Lambda / eigen(A)$vectors[,1] )

# Compute covariance from A, Z vector
Sigma_VAR = matrix(solve(diag(1,n*n,n*n) - kronecker.prod(A)) %*% vec(Z), ncol = 4)

# Compute cross-covariance
Cross_VAR = A %*% Sigma
Cross_LMI = Sigma * delta

# Compute multiple cross-covariances
Cross_VARs = list()
for( i in 1:4) {Cross_VARs[[i]] = A %^% i %*% Sigma}
Cross_LMIs = list()
for( i in 1:4) {Cross_LMIs[[i]] = (delta^i) * Sigma}
