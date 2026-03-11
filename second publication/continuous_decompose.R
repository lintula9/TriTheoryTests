#' Decompose Stationary Covariance of a Continuous-Time VAR (OU) Model
#' @param A Matrix: The coefficient (drift) matrix. Must be Hurwitz (stable).
#' @param Q Matrix: The diffusion matrix (Sigma * Sigma').
continuous_decompose <- function(A, Q) {
  
  # 1. Stability Check: Real parts of eigenvalues must be negative
  ev_A <- eigen(A)$values
  if (any(Re(ev_A) >= 0)) {
    stop("The system is not stationary. Real parts of eigenvalues of A must be negative.")
  }
  
  n <- nrow(A)
  I <- diag(n)
  
  # 2. Solve Lyapunov Equation: A %*% Sigma + Sigma %*% t(A) + Q = 0
  # Using the Kronecker sum identity for the vectorization:
  # (I %x% A + A %x% I) %*% vec(Sigma) = -vec(Q)
  L <- (I %x% A) + (A %x% I)
  vec_sigma <- solve(L, -as.vector(Q))
  stationary_cov <- matrix(vec_sigma, nrow = n)
  
  # Ensure symmetry (correcting for minor numerical rounding)
  stationary_cov <- (stationary_cov + t(stationary_cov)) / 2
  
  # 3. Eigendecomposition
  e_decomp <- eigen(stationary_cov)
  vals <- e_decomp$values
  vecs <- e_decomp$vectors
  
  # 4. Compute Proportion Explained
  # Proportion explained by the first principal component
  prop_pc1 <- vals[1] / sum(vals)
  
  # Return a structured list
  return(list(
    stationary_covariance = stationary_cov,
    eigenvalues = vals,
    eigenvectors = vecs,
    prop_explained_pc1 = prop_pc1
  ))
}