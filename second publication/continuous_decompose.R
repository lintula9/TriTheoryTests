#' Decompose Stationary Covariance of a VAR Model
#' @param A Matrix: The coefficient (drift) matrix for continuous, or AR(1) matrix for discrete.
#' @param Q Matrix: The diffusion matrix / innovation covariance matrix.
#' @param type Character: "continuous" or "discrete". Defaults to "continuous".
decompose_stationary_cov <- function(A, Q, type = c("continuous", "discrete")) {
  
  # Match the type argument to ensure valid input
  type <- match.arg(type)
  
  ev_A <- eigen(A)$values
  n <- nrow(A)
  
  if (type == "continuous") {
    # 1. Stability Check: Real parts of eigenvalues must be negative
    if (any(Re(ev_A) >= 0)) {
      stop("The continuous-time system is not stationary. Real parts of eigenvalues of A must be negative.")
    }
    
    # 2. Solve Continuous Lyapunov Equation: A %*% Sigma + Sigma %*% t(A) + Q = 0
    I <- diag(n)
    L <- (I %x% A) + (A %x% I)
    vec_sigma <- solve(L, -as.vector(Q))
    
  } else if (type == "discrete") {
    # 1. Stability Check: Modulus of eigenvalues must be less than 1
    # Mod() computes the absolute value / modulus of complex numbers in R
    if (any(Mod(ev_A) >= 1)) {
      stop("The discrete-time system is not stationary. Eigenvalues of A must be strictly inside the unit circle.")
    }
    
    # 2. Solve Discrete Lyapunov Equation: Sigma - A %*% Sigma %*% t(A) = Q
    # Note: We need an identity matrix of size n^2 for the vectorization
    I_n2 <- diag(n^2)
    L <- I_n2 - (A %x% A)
    vec_sigma <- solve(L, as.vector(Q))
  }
  
  stationary_cov <- matrix(vec_sigma, nrow = n)
  
  # Ensure symmetry (correcting for minor numerical rounding)
  stationary_cov <- (stationary_cov + t(stationary_cov)) / 2
  
  # 3. Eigendecomposition
  e_decomp <- eigen(stationary_cov)
  vals <- e_decomp$values
  vecs <- e_decomp$vectors
  
  # 4. Compute Proportion Explained
  prop_pc1 <- vals[1] / sum(vals)
  
  # Return a structured list
  return(list(
    stationary_covariance = stationary_cov,
    eigenvalues = vals,
    eigenvectors = vecs,
    prop_explained_pc1 = prop_pc1,
    type = type
  ))
}
