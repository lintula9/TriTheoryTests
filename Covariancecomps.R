#For VAR(1) and a dynamic factor model, we can obtain distinguishable (cross-) covariance.

# Define the factor loadings
Lambda = c(1,2,3)

# Compute within time point covariance from the factor loadings, assuming variance of C-F is 1
Sigma = Lambda %*% t(Lambda)

# Define covariance for the subsequent C-Fs
psi = 0.5

# Compute A = C + B 
C = psi * Sigma * as.vector(( t(Lambda) %*% Lambda )^( -1 ))
B = matrix(c(-2,-3,0,1,0,-3,0,1,2), ncol = 3)
A = C + B

# Solve for Z
Z = Sigma * (1-psi^2)
I = diag(1,ncol = 3*3,nrow = 3*3)

# Ascertain that the within time point covariance equals for both models:
# VAR(1) within time point covariance
Sigma_VAR = matrix(solve(I - fastmatrix::kronecker.prod(A)) %*% fastmatrix::vec(Z), ncol = 3, nrow = 3)

Sigma == round(Sigma_VAR, 10) #Rounding is necessary due to floating point error

# Ascertain that between time point cross-covariance equals for both models:
#VAR(1) cross-covariance
A%*%Sigma
A%*%A%*%Sigma # And so forth.

#C-F cross-covariance
Sigma * psi 
Sigma * psi^2 # And so forth.

