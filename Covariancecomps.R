# Settings ----
require(igraph) # igraph dependency, for graph visualizations
require(expm) # exmp dependency, for matrix exponent computations
require(nloptr) # nloptr depndency, for more general optimization procedures (if needed)

# A numeric example for VAR(1) indistinguishable from a CF model. ----

set.seed(1337)

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

round(Sigma, 10) == round(Sigma_VAR, 10) #Rounding to ensure floating point errors do not influence the result.

# Ascertain that between time point cross-covariance equals for both models:
#VAR(1) cross-covariance
round(A%*%Sigma, 10) == round(Sigma * psi,10) 
round(A%*%A%*%Sigma,10)  == round(Sigma * psi^2 ,10)
round(A%*%A%*%A%*%Sigma,10) == round(Sigma * psi^3,10)
  # And so forth ..

# Plot A to demonstrate how misleadingly complex the Network might look like.


# Create an adjacency matrix for non-zero entries in A to represent edges
adj_matrix <- (A != 0) * 1

# Convert the adjacency matrix to a graph
g <- graph_from_adjacency_matrix(adj_matrix, mode = "directed", weighted = TRUE)

# Add weights (the values in A) as edge attributes
E(g)$weight <- A[which(A != 0)]

# Plot the graph with weights and node labels
plot(g, edge.label = round(E(g)$weight, 2),
     edge.arrow.size = 0.5, 
     vertex.label = c("Var1", "Var2", "Var3"),
     vertex.size = 30,
     layout = layout_in_circle,
     main = "Graphical Model Representation of Matrix A",
     edge.curved = 0.3)

# Compute distance from a VAR(1) to the nearest VAR(1) indistinguishable from a dynamic CF model. ----

  # Define the estimated VAR(1) [here, above example is used, with randomness added to it]

A_hat = A + stats::rWishart(1, df = 5, diag(0.1, ncol = 3, nrow = 3))[,,1]
Z_hat = Z + stats::rWishart(1, df = 5, diag(0.1, ncol = 3, nrow = 3))[,,1]

  # Stationary VAR(1) cross-covariance function
ccov_var = function(Delta, # Distance in time, integer
                    coef_mat, # 'A' matrix
                    inno_cov) { # Innovation covariance 'Z'
  K = nrow(coef_mat)
  Sigma_hat = matrix(solve(I - fastmatrix::kronecker.prod(coef_mat)) %*% 
                       fastmatrix::vec(inno_cov),
                        ncol = K, nrow = K)
  return(A_hat %^% Delta %*% Sigma_hat) }


  # Obtain VAR(1) estimated covariance, for arbitrary many time points T.

T_ = 9 # Set manually
res_varhat = sapply(0:T_, function(x) ccov_var(x, 
                                  coef_mat = A_hat, 
                                  inno_cov = Z_hat),
       simplify = "array")

  # Compute dynamic factor model using the resulting cross-covariances, up to T.
  # Identification is by setting CF variance 1 at time point 1.

    # Set initial values, manually
K_ = nrow(res_var[,,1])
params = c(Lambda = rep(0.5,times = K_), 
           delta = 0.5)
cf_var_1 = 1 # Variance for the CF at time point 1; assumed stationary.


      #### BELOW UNFINISHED
      #### NEED TO UNDERSTAND THE GRADIENT


  # Optimize the Frobenius norm of the difference, w.r.t. L, delta:

criterion_function = function(p = params) {
  
  ( p[1:K_] %*% t(p[1:K_]) * p["delta"] ^ Delta ) - res_varhat
  
  sqrt((sapply(0:T_, function(x) ccov_cf(x,
                                   L = p[1:K_], 
                                   d = p["delta"]),
         simplify = "array") - res_varhat)^2) }

  # Gradient function



  # Optimize the  of the criterion function

cf_model = nloptr::nloptr(x0 = params, 
               eval_f = criterion_function,
               opts = list(
                 "algorithm" = "NLOPT_LN_NELDERMEAD",  # Use Nelder-Mead for derivative-free optimization
                 "xtol_rel" = 1e-6,                    # Relative tolerance on parameters
                 "ftol_rel" = 1e-6,                    # Relative tolerance on function value
                 "maxeval" = 1000,                     # Maximum number of function evaluations
                 "print_level" = 1                     # Print optimization progress (0 for no output)
               ))







