Lambda = c(1,2,3)
Sigma = Lambda %*% t(Lambda)
psi = 0.5
C = psi * Sigma * as.vector(( t(Lambda) %*% Lambda )^( -1 ))
B = matrix(c(-2,-3,0,1,0,-3,0,1,2), ncol = 3)
A = C + B
Z = Sigma * (1-psi^2)

Sigma_VAR = matrix(solve(I - kronecker.prod(A)) %*% vec(Z), ncol = 4, nrow = 4)
