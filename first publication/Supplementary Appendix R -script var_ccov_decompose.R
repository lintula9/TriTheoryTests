# Readme: 


# First, run (source) the whole script. 
# var_ccov_decompose takes the coefficient matrix and innovation covariance as arguments.
# You can define time_points used in the computation and ask to compute all congruencies.
# RMSEA approximation function is also provided (experimental).

#####################################################
# See numerical examples at the end of this script! #
#####################################################


# Helper functions -------------
# Dependencies:
if(!require(fastmatrix)) install.packages("fastmatrix"); library(fastmatrix)
if(!require(expm)) install.packages("expm"); library(expm)
if(!require(pbapply)) install.packages("pbapply"); library(pbapply)
if(!require(pbapply)) install.packages("mgcv"); library(mgcv)
if(!require(pbapply)) install.packages("MASS"); library(MASS)




# Helper function, for the (cross-)covariance
var_ccov <- function(coefmat,innocov,Delta) {

  if(any(Re(eigen(coefmat)$values) > 1)) {print(simpleError("A is not stationary. Aborting."))
    stop()}
  
  # Yule-Walker equation solution for stationary A.
  return( covmat = (coefmat %^% Delta) %*% 
            matrix(Matrix::solve(diag(1,ncol=ncol(coefmat)^2,nrow=nrow(coefmat)^2) - 
                                   Matrix::kronecker(coefmat,coefmat)) %*% 
                     fastmatrix::vec(innocov), 
                   ncol = ncol(coefmat), nrow = nrow(coefmat)))
}

# Helper functions to form M times M block matrices, of the (cross-)covariances. I.e., TK times TK covariance matrices.
var_ccov_Mblock <- function(A, Z, M) {
  K <- nrow(A)      # dimension of each block (A is K x K)
  # Placeholder.
  S_big <- matrix(0, nrow = M * K, ncol = M * K)
  
  for (i in seq_len(M)) {
    for (j in seq_len(M)) {
      
      delta <- abs(i - j)
      
      block_ij <- var_ccov(A, Z, Delta = delta)
      
      if (i > j) {
        block_ij <- t(block_ij)
      }
      
      row_idx <- ((i-1)*K + 1) : (i*K)
      col_idx <- ((j-1)*K + 1) : (j*K)
      
      S_big[row_idx, col_idx] <- block_ij
    }
  }
  
  return(S_big)
}

dcf_ccov_Mblock <- function(L, psi, M) {
  K       <- length(L) 
  # Placeholder.
  S_big <- matrix(0, nrow = M * K, ncol = M * K )
  
  for (i in 1:M) {
    for (j in 1:M) {
      
      delta <- abs(i - j)
      block_ij <- if(delta == 0) tcrossprod(L) else psi^delta * tcrossprod(L)
      
      if (i > j) {
        block_ij <- t(block_ij)
      }
      
      row_idx <- ((i-1)*K + 1) : (i*K)
      col_idx <- ((j-1)*K + 1) : (j*K)
      
      S_big[row_idx, col_idx] <- block_ij
    }
  }
  
  return(S_big)
}

                # Main functions --------------

# Closest VAR(1) model - indistinguishable from a dynamic common factor model. 

civ_find <- function(A, Z, 
                     n.iter = 2000, tol = 1e-6, 
                     W = NULL,
                     random.init    = T,
                     cov.difference = T,
                     time_points = 2,
                     N = NULL,
                     error_ratio = 0.5) {
  
  if(time_points < 2) simpleError("time_points must be 2 or larger.")
  K <- ncol(A)
  
  {

    # Initial values for optimizers.
  if(!random.init){
    Lambda <- Re(eigen(A)$vectors[,1])
    psi    <- Re(eigen(A)$values[1])} else {
    Lambda <- rnorm(K,0.5)
    psi    <- runif(1, 0.5, 0.95)    }
  
  omega0   <- rexp(K, rate = 5)
  
  }
  
  if(!cov.difference){
    
  message("Minimizing distance in terms of parameters (not covariance), as cov.difference was set to False.")
    
    loss_function <- function(pars, A = A, Z = Z) {
      z <- vech(Z)
      s <- c(vec(A),z)
      tilde_z = vech( (1-pars[K+1]^2) * (pars[1:K] %*% t(pars[1:K])) )
      C = pars[K+1]*pars[1:K] %*% t(pars[1:K]) * sum(pars[1:K]^2)^-1
      B = (A - C) %*% (diag(1, nrow = K, ncol = K) - pars[1:K] %*% t(pars[1:K]) * sum(pars[1:K]^2)^-1)
      tilde_s = c(vec( C + B ), tilde_z)
      if(is.null(W)) W = diag(1, nrow = length(s), ncol = length(s))
      M = t(s-tilde_s) %*% W %*% (s-tilde_s)
      return( as.numeric(M) )
    }
    
    gradient_function <- function(pars, A = A, Z = Z) {
      numDeriv::grad(func = loss_function, x = pars, A = A, Z = Z)
    }
    
    params = c(Lambda,psi)
    names(params) <- c(paste0("Lambda", 1:K), "psi")
    optim_result <- optim(par = params, 
                          fn = loss_function, 
                          gr = gradient_function, 
                          method = "BFGS",
                          control = list(abstol = tol, maxit = n.iter),
                          A = A, Z = Z) 
    
    optimized_params <- optim_result$par
    Lambda_opt <- optimized_params[1:K]
    psi_opt    <- optimized_params[K + 1]
    
    C_result <- psi_opt * Lambda_opt %*% t(Lambda_opt) * (sum(Lambda_opt^2))^-1 
    B_result <- (A - C_result) %*% (diag(1,ncol = ncol(A),nrow = nrow(A)) - 
                                      Lambda_opt %*% t(Lambda_opt) * sum(Lambda_opt^2)^-1)
    A_result <- C_result + B_result
    Z_result <- (1 - psi_opt^2) * Lambda_opt %*% t(Lambda_opt)
    
    gc(verbose = F)
    
    # Return results
    return(list(
      "Loadings" = Lambda_opt, 
      "psi"      = psi_opt,
      "A"        = A_result,
      "Z"        = Z_result,
      "Optim_Result" = optim_result
    ))}
  
  
  
  
  if(cov.difference) {
    # Population covariance.
    S  <- var_ccov_Mblock(A,Z,time_points)
    
    if(any(eigen(S)$values  < 1e-7)) {
      message("The VAR(1) predicted covariance has near zero or negative eigenvals. Stopping.\n
              This might be a sign of low rank structure in the predicted covariance. \n
              Check eigen(var_ccov(A,Z,0))$values, for example, for low rank structure.")
      stop()
    }
    
    diag(S) <- diag(S) + (error_ratio * diag(S))
    
    # Maximum likelihood criterion
    F_criterion <- function(theta, S, K, return_cov = F){
      # 1) Parse the parameter vector:
      psi    <- theta[1]                   # scalar
      Lambda <- theta[2:(K+1)]             # length K
      omega0 <- theta[(K+2):(2*K+1)]       # length K
      
      A_indistinguishable = psi * tcrossprod(Lambda) / sum(Lambda^2) # A_indistinguishable
      Z_indistinguishable = (1-psi^2) * tcrossprod(Lambda)

      # 2) Compute the implied covariance.
      S_implied <- dcf_ccov_Mblock(Lambda, psi, time_points) + 
        diag(rep(omega0, times = time_points), time_points*K, time_points*K)
      
      if(return_cov) {return(S_implied)}
      
      if(any(eigen(S_implied)$values < 1e-6)) {warning("Non positive-definite implied covariance.")}

      # 3) Compute the discrepancy.
      return(F_ML = log(det(S_implied)) + 
                    sum(diag( solve(S_implied) %*% S )) - 
                    log(det(S)) - 
                    time_points*K) 
    }
    
    # Maximum Likelihood statistic, with constraints:
    # Error variances (perturbations in the Supplementary Material) and serial covariances positive.
    # The CF is implicitly assumed to have unit variance.
    # DCF is constrained stationary.
    lower_bounds  <- c(-1, rep(-Inf, times = K), rep(0,   times = K))
    upper_bounds  <- c( 1, rep( Inf, times = K), rep(Inf, times = K))
    F_ML  = optim(
      par = c(psi, Lambda, omega0 ), 
      fn  = F_criterion,
      S   = S,
      K   = K, 
      method = "L-BFGS-B", 
      lower  = lower_bounds,
      upper  = upper_bounds,
      control = list(maxit = n.iter))
    
    # Return T statistic.
    statistic = (N-1) * F_ML$value

    # Compute degrees of freedom.
    DF    = time_points*K*(time_points*K + 1)/2 - (K + 1)
    
    # Compute RMSEA.
    RMSEA = sqrt( max(0, statistic - DF) / (DF*(N-1)) )
    
    # Return implied covariance.
    S_implied <- F_criterion(theta = F_ML$par, S = S, K = K, return_cov = T)
    
    # Warnings.
    if(F_ML$value < 0) warning("Negative discrepancy observed.\nResults are not interpretable")
    
    # Return results
    return(list(
      "ML_optimizer_result"           = F_ML,
      T_statistic                     = statistic,
      RMSEA                           = RMSEA,
      CF_ar_coef                      = F_ML$par[1],
      Factor_loadings                 = F_ML$par[2:(K+1)],
      within_time_point_errorvars     = F_ML$par[(K+2):(2*K+1)],
      DF                              = DF,
      "Implied_covariance"            = S_implied,
      "VAR_original_covariance"       = S
    ))
    
  }

}

# Experimental: RMSEA approximation for two subsequent time points.
RMSEA_approx <- function(A, Z, N = NULL, error_ratios = seq(0.01, 1, length.out = 10), ...) {

  rmsea_samples <- as.vector(pbapply::pbsapply(error_ratios, 
                                                   FUN = function(x) {
                                                     rmsea <- try(civ_find(A, Z, error_ratio = x, 
                                                                           N = N, cov.difference = T, ...)$RMSEA,
                                                                  silent = T)
                                                     return(if(is.character(rmsea)) NA else rmsea)}))
  
  # Fit a polynomial regression model
  fit <- try(lm(rmsea_samples ~ error_ratios + I(error_ratios^2) + I(error_ratios^3)), silent = TRUE)
  
  pred_error_ratios <- seq(0, max(error_ratios), length.out = 30)
  predicted <- if (inherits(fit, "try-error")) {rep(NA, length(pred_error_ratios))  # Return NA if fitting fails
    } else {predict(fit, newdata = data.frame(error_ratios = pred_error_ratios))}
  
  
  # Prepare result as a list
  result <- list(RMSEA_approximation = predicted[1],
                 rmsea_samples       = rmsea_samples, 
                 predicted           = predicted, 
                 pred_error_ratios   = pred_error_ratios)
  class(result) <- "rmsea_approximation"
  attr(result, "error_ratios") <- error_ratios
  
  return(result)
  }
  # Plot method for the RMSEA approximation.
plot.rmsea_approximation <- function(x, ...) {
  if (!inherits(x, "rmsea_approximation")) {
    stop("Input must be of class 'rmsea_approximation'.")
  }
  
  error_ratios  <- attr(x, "error_ratios")
  rmsea_samples <- x$rmsea_samples
  pred_error_ratios <- x$pred_error_ratios
  predicted     <- x$predicted
  
  # Create the plot
  plot(error_ratios, rmsea_samples,
       xlab = "Error Ratio", ylab = "RMSEA", 
       main = "RMSEA vs. Error Ratio",
       ylim = c(0, max(max(rmsea_samples, na.rm = TRUE), 0.25)),
       xlim = c(0, max(error_ratios)),
       pch = 19)
  
  # Add quadratic regression line
  lines(pred_error_ratios, predicted, lty = 2)
  
  # Add a cross mark to what RMSEA we predict.
  points(x = 0, y = predicted[1], pch = "x", cex = 1.5)
  
  # Add RMSEA approximation text
  text(0.2, max(max(rmsea_samples, na.rm = TRUE), 0.25) - 0.01, 
       bquote(RMSEA %~~% .(round(predicted[1], 3))))
  
  # Add grid for better readability
  grid() 
  }

# Main function of article: var_ccov_decompose --------------

# var_ccov_decompose
var_ccov_decompose <- function(A,Z,time_points  = 10,
                               all_congruencies = F) {
  
  if(any(abs(eigen(A)$values) > 1)) simpleError("Non-stationary A, aborting.")
  if(any(round(eigen(var_ccov(A,Z,0))$values, digits = 10) == 0)) simpleWarning("Near zero eigenvalues detected in predicted the within time point covariance. Generalized inverse used where needed.")
  
  comp        <- lapply(0:time_points,         function(d){return(eigen(var_ccov(A,Z,d)))})
  all_sum     <- sapply(1:ncol(A), function(j){
                        sum(sapply(1:(time_points+1), function(i) comp[[i]]$values[j]))})
  total_sum   <- sum(sapply(1:(time_points+1), function(i) sum(comp[[i]]$values)))
  cosine_i    <- sapply(1:(time_points+1),
    function(j) {
    sapply(1:(time_points+1),
           function(i) abs( sum( comp[[i]]$vectors[,1] * comp[[j]]$vectors[,1] ) ) ) },
    simplify = "matrix")

  # Multiple cosine similarities (congruencies)
  multiple_congruencies = NULL
  if(all_congruencies) {
    multiple_congruencies  <- lapply(0:time_points,
        FUN = \(i) crossprod((var_ccov(A,Z,i)   |> eigen())$vectors,
                             (var_ccov(A,Z,i+1) |> eigen())$vectors) |>
          abs() # When imaginary eigenvectors are analysed, this is the Hermitian angle.
      ) |> setNames(paste0("(",0:time_points,", ",
                           1:(time_points + 1),")"))
    # Angular similarities for a standardized 
    angular_similarities     <- lapply(multiple_congruencies, \(x) 1 - (acos(x)/pi))
  }

  eigenvals     <- sapply(0:time_points, function(t) eigen( var_ccov(A,Z,t) )$values)
  singularvals  <- sapply(0:time_points, function(t) svd(   var_ccov(A,Z,t  ) )$d )
  # Canonical correlations between x_t and x_{t+Delta}: singular values of the whitened cross-covariance
  Sigma0_inv_sqrt <- {
    e    <- eigen((var_ccov(A,Z,0) + t(var_ccov(A,Z,0)))/2, symmetric = TRUE)
    keep <- e$values > 1e-8 * max(e$values)
    e$vectors[, keep, drop = FALSE] %*% diag(1/sqrt(e$values[keep]), sum(keep)) %*% t(e$vectors[, keep, drop = FALSE])
  }
  ccorrcoefs    <- sapply(0:time_points, function(t) svd( Sigma0_inv_sqrt %*% var_ccov(A,Z,t) %*% Sigma0_inv_sqrt )$d )

  colnames(eigenvals)    <- paste("Increment ",0:time_points)
  colnames(singularvals) <- paste("Increment ",0:time_points)
  colnames(ccorrcoefs)   <- paste("Increment ",0:time_points)

  result <- list(
    eigenvals                     = eigenvals,
    singularvals                  = singularvals,
    canonical_correlations        = ccorrcoefs,
    min_factor_congruency         = min(cosine_i),
    all_factor_congruencies       = cosine_i,
    subsequent_pair_congruencies  = mgcv::sdiag(cosine_i,1),
    all_congruencies              = multiple_congruencies)
  
  class(result) <- c("var_ccov_decompose", "list")
  
  return(result)
}

plot.var_ccov_decompose <- function(x, ...) {
  answer <- readline(
  "What do you want to plot?
    1: eigenvalues,
    2: subsequent congruencies,
    3: singular values,
    4: Canonical correlation,
    5: congruency horizon.")
  if(answer == 1)  {
    matplot((t(x$eigenvals)), type = "b", ylab = "Eigenvalue", 
                           xlab = expression(paste("Increment in time ", Delta, "t")),
                           xaxt = "n", 
            ...)
    axis(1, at  = 1: length(x$eigenvals), 
         labels = 0:(length(x$eigenvals)-1) )
    }
  if(answer == 2)  matplot(x$subsequent_pair_congruencies,  ylim = c(0,1), type = "b", ylab = "Congruency coefficient", 
                           xlab = "T, T+1", ...)
  if(answer == 3)  {
    matplot(t(x$singularvals), type = "b", ylab = "Singular value", 
            xlab = expression(paste("Increment in time ", Delta, "t")),
            xaxt = "n",
            ...)
    axis(1, at  = 1:length(x$singularvals), 
         labels = 0:(length(x$singularvals)-1) )
  }
  if(answer == 4)  {
    matplot(t(x$canonical_correlations), type = "b", ylab = "Canonical correlation coefficient", 
            xlab = expression(paste("Increment in time ", Delta, "t")),
            xaxt = "n",
            ...)
    axis(1, at  = 1:length( x$canonical_correlations), 
         labels = 0:(length(x$canonical_correlations)-1) )
  }
  if(answer == 5)  matplot(x$all_factor_congruencies[1,], ylim = c(0,1), type = "b", ylab = "Congruency coefficient", 
                           xlab = expression(paste("T = 0, T = ", Delta)), ...)
  }


# Numerical examples, also used in main text. -------------------------------

if(F){
  
  # VAR(1), distinguishable
    # Create a VAR(1) model, which rotates and scales.
    # Rotation violates indistinguishability conditions w.r.t. a one CF model.
  Rotation <- matrix(c(
   cos(90*pi/180), -sin(90*pi/180), 0.0, 0.0, 0.0, 0.0, 0.0,
    sin(90*pi/180), cos(90*pi/180), 0.0, 0.0, 0.0, 0.0, 0.0,
                0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
                0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
                0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
                0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
                0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0), 
     nrow = 7, byrow = T)
  Scaling <- diag(seq(0.7,0.4,length.out = 7))
  A_2     <- Rotation %*% Scaling
  Z_2     <- diag(7)
  
  # VAR(1), indistinguighable
  lambdas <- tcrossprod(seq(0.1,0.7,length.out = 7))
  A_3     <- (0.5 *     lambdas ) # + matrix(rnorm(7*7, sd = 0.1), ncol = 7) 
  Z_3     <- tcrossprod(lambdas)  # + matrix(rnorm(7*7, sd = 0.1), ncol = 7) 
  
  # Summary:
    # Figure shown in main text:
  tiff(filename = "Figure_2_tentative.tiff", 
       width    = 17, 
       height   = 19, 
       units    = "cm", 
       res      = 300,
       pointsize = 10)
  
  par(mfrow = c(2,2))
  par(mar   = c(4,4,2,0.5))
  if(!requireNamespace("viridisLite")) {install.packages("viridisLite")
    library(viridisLite) } else library(viridisLite)
  
  #A
  parallel_A <- var_ccov_decompose(A_2,Z_2,all_congruencies=T)
  matplot(t(parallel_A$singularvals), type = "n", 
          ylab = "Singular value", 
          main = "Distinguishable cross-covariance",
          font.main = 1,
          col  = cividis(7),
          xlab = expression(paste("Increment in time ", Delta, "t")),
          xaxt = "n"); grid()
  axis(1, labels = paste0(0:10),
       at = 1:11)
  matplot(t(abs(parallel_A$singularvals)), type = "b",
          col  = cividis(7), add = T,
          lty = 1)
  
  #B
  parallel_B <- var_ccov_decompose(A_3, Z_3)
  matplot(t(parallel_B$singularvals), type = "n", 
          ylab = "",
          main = "Perfectly indistinguishable cross-covariance",
          font.main = 1,
          col  = cividis(7),
          xlab = expression(paste("Increment in time ", Delta, "t")),
          xaxt = "n"); grid()
  axis(1, labels = paste0(0:10),
       at = 1:11)
  matplot(t(abs(parallel_B$singularvals)), type = "b", 
          col  = cividis(6),
          add  = T,
          lty = 1)
  
  #C
  matplot( parallel_A$subsequent_pair_congruencies, type = "n",
          ylab = "Congruency coefficient", 
          xlab = expression(paste("Cross-covariance pair")),
          ylim = c(0,1), 
          main = "Unstable factor loadings",
          font.main = 1,
          xaxt = "n"); grid()
  axis(1, labels = paste0("(",0:10,", ",1:11,")"),
       at = 1:11, cex.axis = 0.7)
  matplot( parallel_A$subsequent_pair_congruencies,
        type   = "b", col = cividis(6), add = T
        )
  
  #D
  matplot( parallel_B$subsequent_pair_congruencies , 
           ylab = "",
           xlab = expression(paste("Cross-covariance pair")),
           type = "n",
           ylim = c(0,1), main = "Perfectly stable factor loadings",
           font.main = 1,
           col  = cividis(6),
           xaxt = "n"); grid()
  axis(1, labels = paste0("(",0:10,", ",1:11,")"),
       at = 1:11, cex.axis = 0.7)
  matplot( parallel_B$subsequent_pair_congruencies,
           type = "b",
           add  = T,
           col  = cividis(6))
  
  dev.off(); par(mfrow=c(1,1));gc()

  }

# Numerical examples, supplement. -------------------------------

if (F) {
  # VAR(1), indistinguishable from a 2 dimensional CF model.
  # Two components shrink at different rates.
  Scaling <- diag(c(0.7, rep(0.8, times = 6)))
  A_4     <- Scaling
  # Create the innovation covariance via the eigendecomposition.
  direction_1 <- eigen(A_4)$vectors[,1]
  direction_2 <- eigen(A_4)$vectors[,2]
  directions <- cbind(c(1,rep(0,times=6)),
                      c(0,1,rep(0,times=5)))
  uneven_eigens <- matrix(c(
                2.0, 0.0,
                0.0, 1.0), 
     nrow = 2, byrow = T)
  Z_4 <- directions %*% uneven_eigens %*% t(directions)

  parallel_C <- var_ccov_decompose(A_4,Z_4,
                                   all_congruencies = T)

  # VAR(1), in which two components exchange place as largest.
  angle <- 90
  Rotation <- matrix(c(
    cos(angle*pi/180), -sin(angle*pi/180), 0.0, 0.0, 0.0, 0.0, 0.0,
    sin(angle*pi/180),  cos(angle*pi/180), 0.0, 0.0, 0.0, 0.0, 0.0,
                0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
                0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
                0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
                0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
                0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0), 
     nrow = 7, byrow = T)
  # Scaling.
  Scaling <- diag(c(0.8, rep(0.9, times = 6)))
  A_5     <- Rotation %*% Scaling
  # Same innovation covariance as above.
  unit_eigens <- matrix(c(
                1.0, 0.0,
                0.0, 1.0), 
     nrow = 2, byrow = T)
  Z_5 <- directions %*% unit_eigens %*% t(directions)

  parallel_D <- var_ccov_decompose(A_5,Z_5, all_congruencies = T)

  tiff(filename = "Figure_supplement_review.tiff", 
       width    = 17, 
       height   = 19, 
       units    = "cm", 
       res      = 300,
       pointsize = 10)
  
  par(mfrow = c(2,2))
  par(mar   = c(4,4,2,0.5))
  if(!requireNamespace("viridisLite")) {install.packages("viridisLite")
    library(viridisLite) } else library(viridisLite)
  
  # Panel A
  matplot(t(parallel_C$singularvals), type = "n", 
          ylab = "Singular value", 
          main = "Components shrink at different rates",
          font.main = 1,
          col  = cividis(7),
          xlab = expression(paste("Increment in time ", Delta, "t")),
          xaxt = "n"); grid()
  axis(1, labels = paste0(0:10),
       at = 1:11)
  matplot(t(abs(parallel_C$singularvals)), type = "b",
          col  = cividis(7), add = T,
          lty = 1)

  # Panel B
  matplot( parallel_C$subsequent_pair_congruencies , 
           ylab = "",
           xlab = expression(paste("Cross-covariance pair")),
           type = "n",
           ylim = c(0,1), main = "Largest component direction changes",
           font.main = 1,
           col  = cividis(6),
           xaxt = "n"); grid()
  axis(1, labels = paste0("(",0:10,", ",1:11,")"),
       at = 1:11, cex.axis = 0.7)
  matplot( parallel_C$subsequent_pair_congruencies,
           type = "b",
           add  = T,
           col  = cividis(6))
  
  # Panel C
  matplot(t(parallel_D$singularvals), type = "n", 
          ylab = "Singular value", 
          main = "Two components varying in size",
          font.main = 1,
          col  = cividis(7),
          xlab = expression(paste("Increment in time ", Delta, "t")),
          xaxt = "n"); grid()
  axis(1, labels = paste0(0:10),
       at = 1:11)
  matplot(t(abs(parallel_D$singularvals)), type = "b",
          col  = cividis(7), add = T,
          lty = 1)

  # Panel D
  matplot( parallel_D$subsequent_pair_congruencies , 
           ylab = "",
           xlab = expression(paste("Cross-covariance pair")),
           type = "n",
           ylim = c(0,1), main = "Rotation alters congruency",
           font.main = 1,
           col  = cividis(6),
           xaxt = "n"); grid()
  axis(1, labels = paste0("(",0:10,", ",1:11,")"),
       at = 1:11, cex.axis = 0.7)
  matplot( parallel_D$subsequent_pair_congruencies,
           type = "b",
           add  = T,
           col  = cividis(6))
  
  dev.off();gc()

  # One component that changes direction. ----
  {
  angle      <- 90
  Rotation_1 <- matrix(c(
    cos(angle*pi/180), -sin(angle*pi/180), 0.0,
    sin(angle*pi/180),  cos(angle*pi/180), 0.0,
                                 0.0, 0.0, 1.0),
     nrow = 3, byrow = T)
  Rotation_2 <- matrix(c(
    1.0, 0.0, 0.0,
    0.0, cos(angle*pi/180), -sin(angle*pi/180),
    0.0, sin(angle*pi/180),  cos(angle*pi/180)),
     nrow = 3, byrow = T)
  Rotation_3 <- matrix(c(
      cos(angle*pi/180), 0.0, sin(angle*pi/180),
      0.0, 1.0, 0.0,
     -sin(angle*pi/180), 0.0,  cos(angle*pi/180)),
     nrow = 3, byrow = T)
  scaling <- diag(c(0.4,0.9,0.9), ncol = 3, nrow = 3)
  lambdas <- tcrossprod(c(1,rep(0,times=2)))
  Z_6     <- tcrossprod(lambdas)
  A_6     <- (Rotation_1 %*% Rotation_2 %*% Rotation_3 %*% scaling) 

  parallel_E <- var_ccov_decompose(A_6,Z_6, 
                                   all_congruencies = T)
  tiff(filename = "Figure_supplement_horizonvssubsequent.tiff", 
       width    = 17, 
       height   = 19, 
       units    = "cm", 
       res      = 300,
       pointsize = 10)
  
  par(mfrow = c(2,2))
  par(mar   = c(4,4,2,0.5))
  if(!requireNamespace("viridisLite")) {install.packages("viridisLite")
    library(viridisLite) } else library(viridisLite)
    # Panel A
  matplot(t(parallel_E$singularvals), type = "n", 
          ylab = "Singular value", 
          main = "One component is dominant",
          font.main = 1,
          col  = cividis(7),
          xlab = expression(paste("Increment in time ", Delta, "t")),
          xaxt = "n"); grid()
  axis(1, labels = paste0(0:10),
       at = 1:11)
  matplot(t(abs(parallel_E$singularvals)), type = "b",
          col  = cividis(7), add = T,
          lty = 1)
  # Panel B
  matplot( parallel_E$subsequent_pair_congruencies , 
           ylab = "",
           xlab = expression(paste("Cross-covariance pair")),
           type = "n",
           ylim = c(0,1), main = "Subsequent congruencies",
           font.main = 1,
           col  = cividis(6),
           xaxt = "n"); grid()
  axis(1, labels = paste0("(",0:10,", ",1:11,")"),
       at = 1:11, cex.axis = 0.7)
  matplot( parallel_E$subsequent_pair_congruencies,
           type = "b",
           add  = T,
           col  = cividis(6))
    # Panel C

  matplot( parallel_E$all_factor_congruencies[1,] , 
           ylab = "",
           xlab = expression(paste("Cross-covariance pair")),
           type = "n",
           ylim = c(0,1), main = "Congruency development t = 0",
           font.main = 1,
           col  = cividis(6),
           xaxt = "n"); grid()
  axis(1, labels = paste0("(",0,", ",0:10,")"),
       at = 0:10, cex.axis = 0.7)
  matplot( parallel_E$all_factor_congruencies[1,],
           type = "b",
           add  = T,
           col  = cividis(6))
  dev.off();gc()
  }

  # Weak eigenvalue or singular value separation. ----
  {
  scaling <- diag(c(0.9, 0.9, 0.900000001), ncol = 3, nrow = 3)
  Z_7     <- diag(1,nrow=3,ncol=3)
  A_7     <- scaling

  parallel_F <- var_ccov_decompose(A_7,Z_7, time_points = 100,
                                   all_congruencies = T)
  tiff(filename = "Figure_supplement_weakseparation.tiff", 
       width    = 17, 
       height   = 19, 
       units    = "cm", 
       res      = 300,
       pointsize = 10)
  
  par(mfrow = c(2,2))
  par(mar   = c(4,4,2,0.5))
  if(!requireNamespace("viridisLite")) {install.packages("viridisLite")
    library(viridisLite) } else library(viridisLite)
    # Panel A
  matplot(t(parallel_F$singularvals), type = "n", 
          ylab = "Singular value", 
          main = "Three comparably large components",
          font.main = 1,
          col  = cividis(7),
          xlab = expression(paste("Increment in time ", Delta, "t")),
          xaxt = "n"); grid()
  axis(1, labels = paste0(0:100),
       at = 1:101)
  matplot(t(abs(parallel_F$singularvals)), type = "b",
          col  = cividis(7), add = T,
          lty = 1)
  # Panel B
  matplot( parallel_E$subsequent_pair_congruencies , 
           ylab = "",
           xlab = expression(paste("Cross-covariance pair")),
           type = "n",
           ylim = c(0,1), main = "Subsequent congruencies",
           font.main = 1,
           col  = cividis(6),
           xaxt = "n"); grid()
  axis(1, labels = paste0("(",0:100,", ",1:101,")"),
       at = 1:101, cex.axis = 0.7)
  matplot( parallel_E$subsequent_pair_congruencies,
           type = "b",
           add  = T,
           col  = cividis(6))
    # Panel C

  matplot( parallel_F$all_factor_congruencies[1,] , 
           ylab = "",
           xlab = expression(paste("Cross-covariance pair")),
           type = "n",
           ylim = c(0,1), main = "Congruency horizon from t = 0",
           font.main = 1,
           col  = cividis(6),
           xaxt = "n"); grid()
  axis(1, labels = paste0("(",0,", ",0:100,")"),
       at = 0:100, cex.axis = 0.7)
  matplot( parallel_F$all_factor_congruencies[1,],
           type = "b",
           add  = T,
           col  = cividis(6))
  dev.off();gc()
  }


  }
