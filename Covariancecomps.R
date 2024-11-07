# Settings ----

require(igraph) # igraph dependency, for graph visualizations
require(expm) # exmp dependency, for matrix exponent computations
require(nloptr) # nloptr depndency, for more general optimization procedures (if needed)
require(numDeriv) # for numerical gradient solving


# Part 1.A numeric example for VAR(1) indistinguishable from a CF model. ----

set.seed(1337)

# Define the factor loadings
Lambda = c(1,2,3)

# Compute within time point covariance from the factor loadings, assuming variance of C-F is 1
Sigma = Lambda %*% t(Lambda)

# Define covariance for the subsequent C-Fs
psi = 0.5

# Compute A = C + B 
C = psi * Sigma * as.vector(( t(Lambda) %*% Lambda )^( -1 ))
B = matrix(c(-1,-1.5,0,0.5,0,-1.5,0,0.5,1), ncol = 3)
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








# Part 2. Stationary: Compute distance from a (stationary) VAR(1) to the nearest (stationary) VAR(1) indistinguishable from a dynamic CF model. ---- 

 #### STEP A: ####

  # Define the estimated VAR(1) [here, above example is used, with randomness added to it]

K_ = 3 # set manually
A_hat = A + stats::rWishart(1, df = 5, diag(0.05, ncol = K_, nrow = K_))[,,1]
Z_hat = Z + stats::rWishart(1, df = 5, diag(0.1, ncol = K_, nrow = K_))[,,1]


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



  #### STEP B: ####

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





    #### STEP C: Convert the obtained dynamic CF model to VAR(1). #####

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
Sigma_tilde = matrix(solve(I - fastmatrix::kronecker.prod(A_tilde)) %*% fastmatrix::vec(Z_tilde), ncol = K_, nrow = K_)
round(eigen(Sigma_tilde)$values, 10)


# Part 3. Linearly non-stationary: Compute distance from a VAR(1) to the nearest VAR(1) indistinguishable from a dynamic CF model. ---- 

#### STEP A: ####

# Define the estimated VAR(1) as linearly increasing in time.

K_ = 3 # set manually, above as well
A_0 = A_hat # 'Baseline' constant matrix, set manually
A_1 = A_hat / 5 # time-varying part
A_t = function( t ) {
  A_0 + A_1*t   }

  # Z_hat kept constant over time.
Z_hat = Z + stats::rWishart(1, df = 5, diag(0.1, ncol = 3, nrow = 3))[,,1]


# VAR(1) cross-covariance function
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



#### STEP B: ####

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





#### STEP C: Convert the obtained dynamic CF model to VAR(1). #####

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
Sigma_tilde = matrix(solve(I - fastmatrix::kronecker.prod(A_tilde)) %*% fastmatrix::vec(Z_tilde), ncol = K_, nrow = K_)
round(eigen(Sigma_tilde)$values, 10)




















# ALTERNATIVE: Simply estimate both models form the data and then compute the distance between. 
# Finite scenario, with A(t) coefficient linear functions of time.

K_ = 5 # Set manually
T_ = 4 # Set manually


  # Parameter estimates:

  # COMPLETELY random A_hat matrices
A_hats = cbind(stats::rWishart(1, df = 5, diag(0.05, ncol = K_, nrow = K_))[,,1], 
               stats::rWishart(1, df = 5, diag(0.05, ncol = K_, nrow = K_))[,,1], 
               stats::rWishart(1, df = 5, diag(0.05, ncol = K_, nrow = K_))[,,1]) # K times (K*T)
A_array = array(A_hats, #
               dim = c(K_,K_,T_ - 1)) # K times K times T (array)

Z_hats = cbind(Lambda %*% t(Lambda), # need to be set manually, current placeholder.
               Lambda %*% t(Lambda), 
               Lambda %*% t(Lambda)) # K times K times T

Z_array = array(Z_hats, #
                dim = c(K_,K_,T_)) # K times K times T (array) (UNCLEAR if T_ + 1 array)

Lambda = c(0.6,0.9,1.2, 0.3, 0.6) # K times 1
delta_t = c(0.4,0.8,1.1, NA) # T + 1 times 1, last is non-existent (NA)
psi_tt = c(0.8,1,1.2, 1) # T + 1 times 1; note the indexing.


  # Convert dynamic CF model to VAR(1)

C_tildes = array(data = NA, dim = c(K_,K_,T_ - 1))
B_tildes = array(data = NA, dim = c(K_,K_,T_ - 1))
A_tildes = array(data = NA, dim = c(K_,K_,T_ - 1))
Z_tildes = array(data = NA, dim = c(K_,K_,T_ - 1))
for( t in 1:(T_-1) ) {
  
  # Compute C, Z tilde.
  C_tildes[ , , t ] <- delta_t[ t ] * ( (Lambda %*% t(Lambda)) / as.numeric( t(Lambda) %*% Lambda ) )
  Z_tildes[ , , t ] <- (psi_tt[ t + 1  ] - delta_t[ t ]^2 * psi_tt[ t ] ) * (Lambda %*% t(Lambda))

  # Orthogonal projection for B tilde.
  B_1 = A_array[ , , t] - C_tildes[ , , t ]
  B_tildes[ , , t ] = B_1 %*% ( diag(1,nrow=K_,ncol=K_) - Lambda %*% ( t(Lambda)%*% Lambda )^(-1) %*% t(Lambda) )
  
  # Save obtained A tilde.
  A_tildes[ , , t ] = C_tildes[ , , t ] + B_tildes[ , , t] 
  }



  # Compute Forbenius norms of distances
raw_distances <- projected_distances <- c()
for( t in 1:(T_-1)) {
  
  raw_distances[ t ] = sum(diag((C_tildes[ , , t ] - A_array[ , , t ]) %*% t(C_tildes[ , , t ] - A_array[ , , t ])))
  projected_distances[ t ] = sum(diag((A_tildes[ , , t ] - A_array[ , , t ]) %*% t(A_tildes[ , , t ] - A_array[ , , t ])))
  }








# Schumacher 2023, and the δ, Λ approximation.

  # We pick the coefficients from Schumacher (although they are defined not as VAr, but instead of lagged logistic models(s)).

exp_coefmat <- matrix(c( # This ignores the intercept.
  9.94, 0.89, 0.94, 1.22, 1.98, 1.46, 1.06, 1.51, 1.00,
  1.98, 7.31, 2.06, 1.54, 1.50, 1.25, 2.22, 1.60, 1.20,
  3.17, 1.54, 2.65, 1.12, 0.81, 1.78, 1.85, 0.97, 2.29,
  1.34, 0.67, 1.12, 7.12, 1.15, 1.27, 1.22, 0.95, 1.24,
  2.00, 1.94, 1.05, 1.50, 3.92, 1.58, 1.76, 1.44, 1.46,
  1.82, 1.23, 1.31, 1.27, 2.20, 4.74, 3.23, 1.35, 2.37,
  1.07, 1.85, 1.36, 1.13, 1.58, 6.28, 9.93, 1.11, 1.13,
  2.03, 1.31, 1.35, 1.58, 3.32, 1.20, 1.14, 2.86, 2.45,
  1.34, 1.33, 2.07, 1.27, 2.34, 1.52, 1.79, 1.85, 4.78
), nrow = 9, byrow = TRUE) # NOT to be confused with matrix exponential.

colnames(exp_coefmat) <- c("sleep", "reduced pleasure", "psychomotor problems", "change in appetite", "mood", 
                            "reduced self-worth", "suicidal ideation", "tiredness", "concentration problems")
rownames(exp_coefmat) <- c("sleep(t-1)", "reduced pleasure(t-1)", "psychomotor problems(t-1)", 
                            "change in appetite(t-1)", "mood(t-1)", "reduced self-worth(t-1)", "suicidal ideation(t-1)", 
                            "tiredness(t-1)", "concentration problems(t-1)")

  # δ, Λ approximate
coefmat = log( exp_coefmat )
criterion_function_deltalambda = function( x  ) {
  dcf_approx_0 = x[1] * (x[2:length(x)] %*% t(x[2:length(x)])) * as.numeric( t(x[2:length(x)]) %*% x[2:length(x)] )^(-1)
  dcf_approx_1 = (coefmat - dcf_approx_0) %*% ( diag(1, nrow(coefmat)) - ( (x[2:length(x)] / 
                                                                            (as.numeric( t( x[2:length(x)] ) %*% x[2:length(x)] ))) %*% 
                                                                            t(x[2:length(x)])
                                                                           ) 
    )   
  # Return Frobenius norm, squared.                                              
  return( sum(diag((coefmat - (dcf_approx_0 + dcf_approx_1)) %*% t((coefmat - (dcf_approx_0 + dcf_approx_1))))) )  }

gradient_function_deltalambda = function( x ) {
  numDeriv::grad( func = criterion_function_deltalambda, x = x)
  }

params = c( 1, rep( mean( coefmat ), times = ncol(coefmat) ) )
deltalambda_approximation <- nloptr::nloptr(
  x0 = params,
  eval_f = criterion_function_deltalambda,
  eval_grad_f = gradient_function_deltalambda,    # Set gradient function
  opts = list(
    "algorithm" = "NLOPT_LD_LBFGS",    # Using L-BFGS to utilize gradient information
    "xtol_rel" = 1e-6,                 
    "ftol_rel" = 1e-6,                 
    "maxeval" = 1000,                  
    "print_level" = 1  ) )

  # The goes (should) go to 0, because the delta lambda approximation (without further constraints) spans all K x K matrices.
  # Now, of course, if the coefficient matrix grows in linear time such that L(t) * coefmat, L(t) is scalar, then the above approximation
  # is perfect as well.


