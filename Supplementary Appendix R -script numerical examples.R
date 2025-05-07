# Supplementary Appendix: R -script, Part 1.
# Sakari Lintula, sakari.lintula@helsinki.fi
# This R script includes all numerical computations, up to using the methods, which are in another script.

# Notes:
# (1) Numerical instabilities are inevitable in these computations. Package 'Matrix' alleviate the problem to some extent.
# Rounding is still necessary to obtain correct results in many cases.
# (2) Object names mostly follows main text.
# Many objects are overwritten throughout the script. Ensure that you run sections, in whole, of the script.




# Settings and libraries. -------

for (i in  c("qgraph", "expm", "nloptr", "numDeriv", "Matrix", "fastmatrix", "MASS", "mlVAR", "vars")) {
    if( require(i, character.only = T) ) library(i, character.only = T) else {install.packages(i); library(i, character.only = T)}
  }

  # Main text, simulation example. ------

    # We simulate from an indistinguishable VAR(1) model, and then fit a VAR(1).
As <- list()
for( i in 1:4){
  simdat <- cbind(simulateVAR(pars = diag(0.5, 4, 4),   # Diagonal coefficient matrix satisfies indistinguishability.
              Nt        = 14*4,                          # This is the same as in the data used in the empirical example.
              residuals = tcrossprod(c(1,1.5,2,2.5)),
              mean      = c(0,0,0,0)
              ),
              1)
      # Add small white noise.
  simdat[,1:4] <- simdat[,1:4] + MASS::mvrnorm(n = 14*4, mu = rep(0, times = 4), Sigma = diag(0.01,4,4))
      # Fit a VAR(1)
  estimated_var <- vars::VAR(simdat[,1:4], p = 1)
  
      # Plot the estimated VAR.
  
  coefs <- coef(estimated_var)
  As[[i]] <- do.call(cbind, lapply(coefs, FUN = function(x) x[1:4,"Estimate"]))
}
tiff(filename = "Figure_1.tiff", 
     width    = 17, 
     height   = 19, 
     units    = "cm", 
     res      = 300,
     pointsize = 10)

par(mfrow = c(2,2))
par(mar   = c(4,4,2,0.5))
qgraph(As[[1]], layout = "circle")
qgraph(As[[2]], layout = "circle")
qgraph(As[[3]], layout = "circle")
qgraph(As[[4]], layout = "circle")
dev.off()
par(mfrow = c(1,1))

    # Below was removed from the main text: -----------

  # Analytic results, Result 1 & 2. ####
  ## Main text, results Result 1 -section -----------
  
    # Computations for the 'unadjusted' example:
  Lambda = c(1,2,3)
  psi = 1
  phi = 0.5
  
  Sigma_CF = Lambda %*% t(Lambda)
  A = phi *  Lambda %*% solve( t(Lambda) %*% Lambda ) %*% t(Lambda)
  Z = (psi - phi^2) * Sigma_CF
  
  Sigma_VAR = matrix(solve(diag(1, nrow = 9, ncol = 9) - kronecker.prod(A)) %*% vec(Z), ncol = 3)
  Lambda %*% t(Lambda)
  
    #Check for indistinguishability:
  all( round(Sigma_CF, 10) == round(Sigma_VAR, 10) ) # Covariance.
  all( round(phi*Sigma_CF, 10) == round(A %*% Sigma_VAR, 10) ) # Cross-covariance for Delta t = 1.
  all( round(phi^2*Sigma_CF, 10) == round(A %*% A %*% Sigma_VAR, 10) ) # Cross-covariance for Delta t = 2.
      # ... and so forth.
  
    # Computations for the 'adjusted' example:
  
  A_adj = A + matrix( c(0,3,-2,
                        -2,1,0,
                        1,-0.5,0), 
                      ncol = 3, nrow = 3, byrow = T)
  
    #Check for indistinguishability:
  all( round(Sigma_CF, 10) == round(Sigma_VAR, 10) ) # Covariance.
  all( round(phi*Sigma_CF, 10) == round(A_adj %*% Sigma_VAR, 10) ) # Cross-covariance for Delta t = 1.
  all( round(phi^2*Sigma_CF, 10) == round(A_adj %*% A_adj %*% Sigma_VAR, 10) ) # Cross-covariance for Delta t = 2.
      # ... and so forth.
  
  
  
  ## Extension of results section 4 example: More general R script using randomly sampled parameters for the dynamic CF model. ------
  
  set.seed(21) # For replicability.
  
  ### Create a stationary indistinguishable VAR(1) model.
    # Number of dimensions and time points
  K_ = 5
  T_ = 10
  
    # Pick a random sample of factor loadings from the interval [0.1, 2], for example.
  Lambda = runif(K_, min = 0.1, max = 2)
  
  psi_tt = 1 # Set variance of the CF - 1 for example.
  phi_ = 0.5 # Define autoregression coefficient for the subsequent C-Fs - 0.5 for example.
  
    # Compute A_tilde = C_tilde + B_tilde, using the formulas provided in the article.
  C_tilde = phi_ * Lambda %*% t(Lambda) * as.vector(( t(Lambda) %*% Lambda )^( -1 ))
  
    # Create a random B matrix to add to C, so that indistinguishably still holds.
  B = stats::rWishart( 1, df = K_ + 1, diag(0.2, nrow = K_, ncol = K_) )[,,1] # Obtain some random matrix B
  B_tilde = B %*% ( diag(1L, nrow = K_, ncol = K_ ) - (Lambda  %*% (as.vector(( t(Lambda) %*% Lambda )^( -1 )) %*% t(Lambda))) ) # Project it onto the orthogonal complement.
  
    # For stationary of A, we need to ensure that B_tilde has eigenvalues with absolute value less than one.
  if(any( abs(eigen(B_tilde)$values) > 1)) {
    
    while(any( abs(eigen(B_tilde)$values) > 1)) {
      B_tilde <- B_tilde / 1.1
  
    }
    
  }
  
  
  A = C_tilde + B_tilde # Create A
  
  # Compute innovation covariance Z
  Z = Lambda %*% t(Lambda) * psi_tt*(1-phi_^2)
  
  ### Ascertain that the within time point covariance equals for both models:
  
  Sigma_VAR = matrix(solve(diag(1, ncol = K_^2,nrow = K_^2) - fastmatrix::kronecker.prod(A)) %*% fastmatrix::vec(Z), ncol = K_, nrow = K_)
  Sigma_CF = Lambda %*% t(Lambda) * psi_tt
  
    # Check the result:
  all(round(Sigma_CF, 10) == round(Sigma_VAR, 10)) #Rounding to ensure floating point errors do not influence the result.
  
  
    # Ascertain that between time point cross-covariance equals for both models:
  Sigma_CF = Lambda %*% t(Lambda) * psi_tt # This must equal Sigma_VAR above for stationary VAR. 
    #VAR(1) cross-covariance equals the dynamic CF cross-covariance (of symptoms)
  all(round(A%*%Sigma_VAR, 10) == round(Sigma_CF * phi_,10) )
  all(round(A%*%A%*%Sigma_VAR,10)  == round(Sigma_CF * phi_^2 ,10))
  all(round(A%*%A%*%A%*%Sigma_VAR,10) == round(Sigma_CF * phi_^3,10))
      # And so forth ..
  
    # Plot A to demonstrate how the Network looks like.
  qgraph(A, 
         title = "VAR(1) coefficient matrix, indistinguishable from a CF model.",
         title.cex = 1.5,  
         mar = c(4, 4, 6, 4)
  ); gc()
  
  
  # Figure 1. in main text -----------
  library(Matrix); library(qgraph)
  
  tiff(filename = "Figure_1.tiff", width = 8, height = 6, units = "in", res = 480)
  lambda = c(1,2,3)
  A = lambda %*% t(lambda) * c(t(lambda) %*% lambda)^-1 + 
    matrix(c(0,3,-2,
             -2,1,0,
             1,-0.5,0), ncol = 3, nrow = 3, byrow = T)
  Z = (1-0.5^2) * lambda %*% t(lambda)
  
  # Plot.
  labels <- expression(X[1], X[2], X[3])
  par(mfrow=c(1,2))
  coeflabs = expression(tilde(A), tilde(A)[adj.])
  qgraph(lambda %*% t(lambda) * c(t(lambda) %*% lambda)^-1,
         mar = c(4, 4, 6, 4),
         layout = "circle", 
         labels = labels, # Use the expression labels
         directed = T,
         title = coeflabs[1],
         title.cex = 1.5
  )
  qgraph(A,
         mar = c(4, 4, 6, 4),
         layout = "circle", 
         labels = labels, # Use the expression labels
         title = coeflabs[2],
         title.cex = 1.5
  )
  par(mfrow=c(1,1));dev.off();gc()
  
  
  
  
  
  # Main text, Analytic results -section, Result 2: multidimensional dynamic factor model to VAR(1) --------------
  set.seed(21)
  
  # Stationary 2 dimensional D-CF(1) model onto VAR(1).
    # Define the parameters.
  
  L = matrix(runif(6,min = 0.2, max = 1), nrow = K_) # Time invariant CF loadings randomly from uniform ( 0.2, 1).
  
  p_ = 2 # number of CFs.
  A_eta = matrix(c(0.5,0.2,0.2,0.5), ncol = p_) # CF autoregression.
  A_Z = matrix(c(0.2 , 0.0 ,0.0 , 0.2), ncol = p_) # CF innovation covariance
  Sigma_eta = matrix( solve(diag(1, ncol = p_*p_, nrow = p_*p_) - kronecker.prod(A_eta) ) %*% vec(A_Z),
                      ncol = p_, nrow = p_)
  
    # Apply the map.
  A = L %*% A_eta %*% ginv(L) # Remember to ensure we have produced a stationary VAR(1). This required the CF latent autoregression to be stationary.
  Z = L %*% Sigma_eta %*% t(L) - L %*% A_eta %*% Sigma_eta %*% t(A_eta) %*% t(L) 
  
    # Compute the (cross-)covariances.
  Sigma_VAR = matrix( solve(diag(1, K_^2, K_^2) - kronecker.prod(A)) %*% vec(Z),
                      ncol = K_, nrow = K_)
  
  Sigma_DCF = matrix(matrix( L %*% Sigma_eta %*% t(L) ),
                     ncol = K_, nrow = K_)
  
    # Check compatibility.
  round(Sigma_VAR, 10) == round(Sigma_DCF, 10)
  
    # Compute cross-covariance up to Delta t = 2.
  
  K_1_VAR = A %*% Sigma_VAR
  K_2_VAR = A %*% A %*% Sigma_VAR
  K_1_DCF = L %*% ( A_eta %*% Sigma_eta ) %*% t(L)
  K_2_DCF = L %*% ( A_eta %*% A_eta %*% Sigma_eta ) %*% t(L)
  
  round(K_1_VAR, 10) == round(K_1_DCF, 10)
  round(K_2_VAR, 10) == round(K_2_DCF, 10)
  
      #... and so forth.
  
    # Add in B(t), orthogonal to L.
  
  B_0 = rWishart(1, df = K_+1, Sigma = diag(0.2, nrow = K_))[,,1] # Initiate a matrix B_0 randomly.
  B = B_0 %*% ( diag(1, nrow = K_, ncol = K_) - L %*% solve(t(L) %*% L) %*% t(L) )
  
  round(B %*% L, 10)
  
  A_adj = A + B
  
  Sigma_VAR_adj = matrix( solve(diag(1, nrow = K_^2, ncol = K_^2) - kronecker.prod(A_adj)) %*% vec(Z),
                          ncol = K_, nrow = K_)
  K_1_VAR_adj = A_adj %*% Sigma_VAR_adj
  K_2_VAR_adj = A_adj %*% A_adj %*% Sigma_VAR_adj
  
  all(round(K_1_VAR_adj, 10) == round(K_1_VAR, 10))
  all(round(K_2_VAR_adj, 10) == round(K_2_VAR, 10))
  
    #... and so forth.
  
  
  
  
  # Main text, Numerical and empirical examples -section, numerical example ------
  
    # 07.04.2025: This section was removed from the main text.
  # In this section of the script, we plot the Figure of a linearly time-varying VAR(1), which is indistinguishable from a dynamic
  # factor model.
  
  set.seed(123); gc()
  
  # Number of dimensions and time points
  K_ = 9
  T_ = 10
  
    # intialize the CF model parameters, this overwrites previous parameters.
  Lambda = runif(K_, min = 0.5, max = 1)
  psi_tt = 1 
  theta_ = 1.5 
  phi_ = 0.5
  
    # First, create C matrix again (named differently to separate from above.)
  C_t = theta_ * Lambda %*% t(Lambda) * as.vector(( t(Lambda) %*% Lambda )^( -1 ))
  
  # Create B matrix, to produce some asymmetry for our coefficient matrix.
  B_t = stats::rWishart( 1, df = K_, diag(0.2, nrow = K_, ncol = K_) )[,,1] %*%
    ( diag(1, nrow = K_, ncol = K_ ) - (Lambda  %*% (as.vector(( t(Lambda) %*% Lambda )^( -1 )) %*% t(Lambda))) )
  
    # Also ensure we have the correct Z
  Z_0 = Lambda %*% t(Lambda) * (psi_tt-phi_^2)
  
    # Initialize an array to store results
  A_t <- array(0, dim = c( K_, K_, T_ + 1))
  
    # Create a random square matrix, projected onto the orthogonal complement of Lambda,
    # which will then vary linearly in time.
  time_matrix <- stats::rWishart( 1, df = K_, diag(0.3, nrow = K_, ncol = K_) )[,,1] %*%
    ( diag(1, nrow = K_, ncol = K_ ) - (Lambda  %*% (as.vector(( t(Lambda) %*% Lambda )^( -1 )) %*% t(Lambda))) )
  time_matrix %*% Lambda # Ensure it is practically zero.
  
  # Fill in each slice of the array
  for (t in 0:T_) {
    A_t[,,t + 1] <- C_t + ( (B_t * 0.2)  + ((time_matrix / 20 ) * exp((t - 5) / 5 )))
  }
  
  # Plot.
  max_weight <- max(sapply(1:dim(A_t)[3], function(t) max(abs((A_t[,,t])))))
  
  par(mfrow = c(2, 2), oma = c(0, 0, 4, 0)) # Adjust oma for outer margin to accommodate the title
  
  expr_list <- lapply(1:K_, function(i) bquote(X[.(i)]))
  labels <- do.call(expression, expr_list)
  # Plot, if needed:
  if(F) {qgraph(Z_0, title = "'Contemporaneous' covariance",title.cex = 1.5,mar = c(4, 4, 6, 4),layout = "circle", labels = labels );qgraph(A_t[,,1], title = "Lagged effects\nTime point 1",title.cex = 1.5,  mar = c(4, 4, 6, 4),maximum = max_weight, edge.width = 3,   layout = "circle",labels = labels);qgraph(A_t[,,5], title = "\nTime point 5",title.cex = 1.5,  mar = c(4, 4, 6, 4),maximum = max_weight, edge.width = 3, layout = "circle",labels = labels  );qgraph(A_t[,,9], title = "\nTime point 9",title.cex = 1.5,  mar = c(4, 4, 6, 4),maximum = max_weight, edge.width = 3,layout = "circle",labels = labels);mtext("VAR(1) Network model indistinguishable from a dynamic CF model.", outer = TRUE, cex = 1.5, font = 2);par(mfrow = c(1,1))}
  
  ## Adjust slightly, and plot again. ###
  
  # Adjust the time-varying part from above, slightly, to create model 2.
  time_matrix_adj       <- time_matrix
  time_matrix_adj[ 1, ] <- time_matrix[ 1, ] + runif(n = K_, 9,15) %*% (diag(1,nrow=K_,ncol=K_) - Lambda %*% solve(t(Lambda) %*% Lambda) %*% t(Lambda)) 
  time_matrix_adj[ 2, ] <- time_matrix[ 1, ] - runif(n = K_, 9,15) %*% (diag(1,nrow=K_,ncol=K_) - Lambda %*% solve(t(Lambda) %*% Lambda) %*% t(Lambda))
  time_matrix_adj[ 3, ] <- time_matrix[ 1, ] + runif(n = K_, 9,15) %*% (diag(1,nrow=K_,ncol=K_) - Lambda %*% solve(t(Lambda) %*% Lambda) %*% t(Lambda))
  time_matrix_adj[ 4, ] <- time_matrix[ 1, ] - runif(n = K_, 9,15) %*% (diag(1,nrow=K_,ncol=K_) - Lambda %*% solve(t(Lambda) %*% Lambda) %*% t(Lambda))
  time_matrix_adj[ 5, ] <- time_matrix[ 1, ] - runif(n = K_, 9,15) %*% (diag(1,nrow=K_,ncol=K_) - Lambda %*% solve(t(Lambda) %*% Lambda) %*% t(Lambda))
  
  # Adjust the time-varying A_t
  A_t_adj <- array(0, dim = c( K_, K_, T_ + 1))
  for (t in 0:T_) {
    A_t_adj[,,t + 1] <-  C_t + ( (B_t * 0.2)  + ((time_matrix_adj / 20 ) * exp((t - 5) / 5 ) ) )
  }
  
  # Plot the comparison
  max_weight <- max(sapply(1:5, function(t) max(abs((A_t_adj[,,t])))))
  
  par(mfrow = c(3, 2), oma = c(0, 0, 4, 0)) # Adjust oma for outer margin to accommodate the title
  
  
  tiff(filename = "Figure_2.tiff", 
       width    = 17, 
       height   = 19, 
       units    = "cm", 
       res      = 300,
       pointsize = 10) 
  
  par(mfrow = c(3, 2), oma = c(0, 0, 0, 0)) # Adjust oma for outer margin to accommodate the title
  # Plot the Lagged effects at Time point 1 (A_t_adj)
  qgraph(A_t_adj[,,1], 
         title = "Model 1 and 2, time point 0",
         title.cex = 1.5,  
         mar = c(4, 4, 6, 4),
         maximum = max_weight, # Consistent scale for edges
         edge.width = 3,       # Adjust this value for larger edges
         layout = "circle",
         labels = labels # Custom line types
  )
  
  # Empty plot
  plot.new()
  
  # Plot the Lagged effects at Time point 5 (A_t_adj)
  qgraph(A_t_adj[,,6], 
         title = "Model 1, time point 5",
         title.cex = 1.5,  
         mar = c(4, 4, 6, 4),
         maximum = max_weight, # Consistent scale for edges
         edge.width = 3,       # Adjust this value for larger edges
         layout = "circle",
         labels = labels      # Use the expression labels
  )
  
  # Plot at Time point 5 (A_t)
  qgraph(A_t[,,6], 
         title = "Model 2, time point 5",
         title.cex = 1.5,  
         mar = c(4, 4, 6, 4),
         maximum = max_weight, # Consistent scale for edges
         edge.width = 3,       # Adjust this value for larger edges
         layout = "circle",
         labels = labels      # Use the expression labels
  )
  
  # Plot at Time point 9 (A_t_adj)
  qgraph(A_t_adj[,,10], 
         title = "Model 1, time point 9",
         title.cex = 1.5,  
         mar = c(4, 4, 6, 4),
         maximum = max_weight, # Consistent scale for edges
         edge.width = 3,       # Adjust this value for larger edges
         layout = "circle",
         labels = labels      # Use the expression labels
  )
  
  # Plot at Time point 9 (A_t)
  qgraph(A_t[,,10], 
         title = "Model 2, time point 9",
         title.cex = 1.5,  
         mar = c(4, 4, 6, 4),
         maximum = max_weight, # Consistent scale for edges
         edge.width = 3,       # Adjust this value for larger edges
         layout = "circle",
         labels = labels      # Use the expression labels
  ); dev.off(); gc(); par(mfrow = c(1,1))
  

   ### Ascertain the covariance computations.
    # Initialize with time point 0 covariance.
  Sigma_0 = Lambda %*% t(Lambda) * psi_tt
  
  # Initialize a list to store Sigma and Z matrices
  Sigma <- vector("list", 11)
  Sigma[[1]] <- Sigma_0 # Store Sigma_0
  
    # Compute Sigma_k iteratively for t = 1, 2, ..., 10.
  for (t in 1:11) {
    if (t == 1) {
      Sigma[[t + 1]] <- A_t[,,1] %*% Sigma[[t]] %*% t(A_t[,,1]) + Z_0
    } else {
      # Use the appropriate transition matrix for Sigma_k (k > 1)
      Sigma[[t + 1]] <- A_t[,,t] %*% Sigma[[t]] %*% t(A_t[,,t]) + Z_0
    }
  }
    # Be aware: Numerical instabilities make some eigenvalues inevitably deviate from 0, ESPECIALLY when using base R.
  for( i in Sigma) {print(round(eigen(i)$values, 5))}
  
