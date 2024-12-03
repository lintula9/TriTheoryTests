# Result 4.
library(Matrix); library(fastmatrix)

Lambda = c(1,2,3)
psi = 1
phi = 0.5

Sigma = Lambda %*% t(Lambda)
A = phi *  Lambda %*% solve( t(Lambda) %*% Lambda ) %*% t(Lambda)
Z = (psi - phi^2) * Sigma

matrix(solve(diag(1, nrow = 9, ncol = 9) - kronecker.prod(A)) %*% vec(Z), ncol = 3)
Lambda %*% t(Lambda)
