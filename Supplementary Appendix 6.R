# Settings and libraries. ----

for (i in  c("qgraph", "expm", "nloptr", "numDeriv")) {
    if( require(i, character.only = T) ) library(i, character.only = T) else {install.packages(i); library(i, character.only = T)}
  }


# Part 1.A numeric example for VAR(1) indistinguishable from a dynamic CF model. ----

set.seed(21) # One can change the seed to obtain some other result.

# Number of dimensions and time points
K_ = 5
T_ = 10


# Pick a random sample of factor loadings from the interval (0, 2].
Lambda = runif(K_, min = 0.0001, max = 2)

psi_tt = 1 # Set variance of the CF, 1 for example.
delta_ = 0.5 # Define autoregression coefficient for the subsequent C-Fs.

  # Compute A = C + B, as in the article.
C = delta_ * Lambda %*% t(Lambda) * as.vector(( t(Lambda) %*% Lambda )^( -1 ))

  
B = stats::rWishart( 1, df = K_ + 1, diag(0.2, nrow = K_, ncol = K_) )[,,1] # Obtain some random matrix B
B_tilde = B %*% ( diag(1, nrow = K_, ncol = K_ ) - (Lambda  %*% (as.vector(( t(Lambda) %*% Lambda )^( -1 )) %*% t(Lambda))) ) # Project it onto the orthogonal complement.
A = C + B_tilde # Create A

# Compute innovation covariance Z
Z = Lambda %*% t(Lambda) * (psi_tt-delta_^2)

# Ascertain that the within time point covariance equals for both models:

if(all( abs(eigen(A)$values) < 1 ) ) { # VAR(1) within time point covariance - only works for stationary A.
Sigma_VAR = matrix(solve(diag(1, ncol = K_^2,nrow = K_^2) - fastmatrix::kronecker.prod(A)) %*% fastmatrix::vec(Z), ncol = K_, nrow = K_)

# Check the result:
round(Lambda %*% t(Lambda) * psi_tt, 10) == round(Sigma_VAR, 10) #Rounding to ensure floating point errors do not influence the result.

}

# Ascertain that between time point cross-covariance equals for both models:
Sigma_ = Lambda %*% t(Lambda) * psi_tt # This must equal Sigma_VAR above for stationary VAR, otherwise we cannot 'know' it.
#VAR(1) cross-covariance equals the dynamic CF cross-covariance (of symptoms)
round(A%*%Sigma_, 10) == round(Sigma_ * delta_,10) 
round(A%*%A%*%Sigma_,10)  == round(Sigma_ * delta_^2 ,10)
round(A%*%A%*%A%*%Sigma_,10) == round(Sigma_ * delta_^3,10)
  # And so forth ..

# Plot A to demonstrate how the Network looks like.

qgraph(A, 
       title = "VAR(1), indistinguishable from a dynamic CF model.",
       title.cex = 1.5,  
       mar = c(4, 4, 6, 4) 
)

dev.off()

# Part 1.B Brief, simple, extension to linearly time-varying ------

  # Define a simple linear A_t, which is linear in time
# Initialize the array to store results
A_t <- array(0, dim = c( K_, K_, T_ + 1))

# Fill in each slice of the array
for (t in 0:T_) {
  A_t[,,t + 1] <- A + ((A / 2) * t)
}

# Note that A_t has the eigenvalue delta_ + delta_ / 2 * t, for the eigenvector Lambda.

A_t_eigen = lapply( 0:T_, FUN = function(t) eigen( as.matrix(A_t[,,t+1, drop = T]) ) )
for(i in A_t_eigen) { print(i$values[1])} # The eigenvalues.
for(i in A_t_eigen) { print(i$vectors[,1] / Lambda)} # The eigenvectors, seen as scalar multiples of Lambda 
                                                     # (hence, Lambda is also an eigenvector.)


# Part 1.C Figure 1 in main text: Plot A_t, 
# linear in time, indistinguishable from a dynamic CF model.

par(mfrow = c(2,2))
max_weight <- max(sapply(1:dim(A_t)[3], function(t) max(abs((A_t[,,t]))))) / 2

qgraph(Z, 
       title = "Indistinguishable VAR(1) model\nInnovation or 'Contemporaneous' covariance",
       title.cex = 1.5,
       mar = c(4, 4, 6, 4)
)

qgraph(A_t[,,2], 
       title = "\nTime point 1",
       title.cex = 1.5,  
       mar = c(4, 4, 6, 4),
       maximum = max_weight, # Consistent scale for edges
       edge.width = 3 # Adjust this value for larger edges
)



qgraph(A_t[,,6], 
       title = "\nTime point 5",
       title.cex = 1.5,  
       mar = c(4, 4, 6, 4),
       maximum = max_weight, # Consistent scale for edges
       edge.width = 3 # Adjust this value for larger edges
)

qgraph(A_t[,,10], 
       title = "\nTime point 9",
       title.cex = 1.5,  
       mar = c(4, 4, 6, 4),
       maximum = max_weight, # Consistent scale for edges
       edge.width = 3 # Adjust this value for larger edges
)

par(mfrow = c(1,1))


# Part 2. Compute distance from a VAR(1) to the nearest VAR(1) indistinguishable from a dynamic CF model. ----


  # Change below to something else.



  # Define the estimated VAR(1) [here, above example is used, with randomness added to it]

A_hat = A + stats::rWishart(1, df = 5, diag(0.05, ncol = 3, nrow = 3))[,,1]
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
K_ = nrow(res_varhat[,,1])
params = c(Lambda = rep(0.5,times = K_), 
           delta = 0.5)
cf_var_1 = 1 # Variance for the CF at time point 1; assumed stationary.


  # Optimize the Frobenius norm of the difference, w.r.t. L, delta:

criterion_function = function(x) {
  return(
    log(
    sum(
      (sapply(0:T_, function(Delta) ( x[1:K_] %*% t(x[1:K_]) * x[K_+1] ^ Delta ),
                                simplify = "array") - res_varhat)^2)
    )) }

  # Gradient, evaluated numerically.

criterion_gradient = function(x) {
  return(numDeriv::grad(criterion_function, x))
  }
  
  
  

  # Optimize the  of the criterion function, and obtain the CF model

cf_model <- nloptr::nloptr(
  x0 = params,
  eval_f = criterion_function,
  eval_grad_f = criterion_gradient,    # Set gradient function
  opts = list(
    "algorithm" = "NLOPT_LD_LBFGS",    # Using L-BFGS to utilize gradient information
    "xtol_rel" = 1e-6,                 
    "ftol_rel" = 1e-6,                 
    "maxeval" = 1000,                  
    "print_level" = 1  )
)


# Convert the optained dynamic CF model to VAR(1).
L_tilde = cf_model$solution[1 : K_]
delta_tilde = cf_model$solution[K_ + 1]

    # We construe A_tilde = C_tilde + B_tilde. First C_tilde
C_tilde = (delta_tilde * L_tilde %*% ( t(L_tilde)%*% L_tilde )^(-1) %*% t(L_tilde))

    # Then B_tilde
B_1 = A_hat - C_tilde 
B_tilde = B_1 %*% (diag(1,nrow=K_,ncol=K_) - L_tilde %*% ( t(L_tilde)%*% L_tilde )^(-1) %*% t(L_tilde))

    # Create A_tilde
A_tilde = C_tilde + B_tilde

    # Compute frobenius norm of the distance between A_hat and A_tilde.
sqrt(sum(diag( (A_hat - A_tilde) %*% t( A_hat - A_tilde ) )))
      # Observe, that the distance between A_hat and C_tilde only is larger. 
sqrt(sum(diag( (A_hat - C_tilde) %*% t( A_hat - C_tilde ) )))

  # Compute Z_tilde
Z_tilde = (1-cf_model$solution[K_+1]) * L_tilde %*% ( t(L_tilde)%*% L_tilde )^(-1) %*% t(L_tilde)
round(eigen(Z_tilde, only.values = T)$values, 10) # Floating point error corrected positive semi-definiteness check.

  # Compute Frobenius norm of the distance between Z_hat and Z_tilde
sqrt(sum(diag( (Z_hat - Z_tilde) %*% t( Z_hat - Z_tilde ) )))

  # Ensure that the result obtained is indeed CF compatible (by inspecting that the cross-covariance is rank 1 symmetric)
Sigma_tilde  = matrix(solve(I - fastmatrix::kronecker.prod(A_tilde)) %*% fastmatrix::vec(Z_tilde), ncol = K_, nrow = K_)
round(eigen(Sigma_tilde)$values, 10)
