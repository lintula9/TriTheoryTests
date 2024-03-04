# Covariance numerical calulations
library(fastmatrix)
A = matrix(c(.3,.2,.2,.2,
             .2,.3,.2,.2,
             .2,.2,.3,.2,
             .2,.2,.2,.3), ncol = 4, byrow = F)
Psi = diag(0.2, ncol = ncol(A), nrow = nrow(A))
I = diag(1, ncol = ncol(A)*ncol(A), nrow = nrow(A)*nrow(A))
fastmatrix::kro
Sigma = matrix(solve(I - kronecker.prod(A)) %*% vec(Psi), ncol = 4, nrow = 4)
Sigma_12 = A %*% Sigma

