# Packages -----
# nolint start
# List of required packages
required_packages <- c(
  "Matrix", "fastmatrix", "expm", "rstan",
  "qgraph", "tidyverse", 
  "ggplot2", "rstantools", "bayesplot", "cmdstanr", "posterior",
  "viridisLite", "dplyr")

# Function to check and install missing packages
for (pkg in required_packages) {
  if(pkg == "cmdstanr" & !requireNamespace(pkg, quietly = T)) install.packages("cmdstanr", repos = c('https://stan-dev.r-universe.dev', getOption("repos"))) 
  else if(!requireNamespace(pkg, quietly = T)) {install.packages(pkg, dependencies = T)}
    library(pkg, character.only = TRUE)
}
# Read data ----
inference_vars_regex_alpha <- c("A_effective","A","psi","Lambda", "Omega","L_", "cutpoints", "time_of_day_effect",
                                "ref_time_of_day_effect", "specific_time_of_day_effect")
nchains = 8
fit_Net    <- as_cmdstan_fit(files = paste0("Datas/3VAR_CF_relaxedpriors-202502031854-",1:nchains,".csv") ); gc()
draws_data <- as_draws_df(fit_Net$draws(variables = c(inference_vars_regex_alpha) ), .nhcains = nchains ); gc(); rm(fit_Net); gc()

# Compute eigen parallel ---
K = 3
varLabs <- c("Relax", "Worry",  
             "Nervous")
varLabs2 <- c("Relax", "Worry",  
              "Nervous", "Tired", "Hungry",
              "Alone", "Angry")
# Extract posterior means as the VAR parameters:
A <- matrix( unlist(colMeans(draws_data[,grep("A", names(draws_data))])),     ncol = K, nrow = K);
Z <- matrix( unlist(colMeans(draws_data[,grep("Omega", names(draws_data))])), ncol = K, nrow = K)

# Obtain complete VAR(1) model samples:
estimated_var_samples <- draws_data[,c( grep("A", names(draws_data)) , grep("Omega", names(draws_data)) )]; As <- grep("A_", names(estimated_var_samples)); Os <- grep("Omega", names(estimated_var_samples)); 
var_samples           <- pbapply::pblapply(1:nrow(draws_data),
                                 FUN = function(i){
                                   A_temp <- matrix( unlist(estimated_var_samples[i,As]), ncol = sqrt(length(As)), nrow = sqrt(length(Os)));
                                   Z_temp <- matrix( unlist(estimated_var_samples[i,Os]), ncol = sqrt(length(As)), nrow = sqrt(length(Os)));
                                   return(list(A = A_temp, Z = Z_temp))
                                 }); gc()
# Source the methods.
source("./first publication/Supplementary Appendix R -script var_ccov_decompose.R")

## Figure 4 in main text ----
  # Compute parallel analysis imitation and RMSEA, for the posterior mean.
result_parallel  <- var_ccov_decompose(A, Z, time_points = 50)

  # Compute credible intervals for eigenvalues, congruencies.
eigen_congurency <- pbapply::pblapply(var_samples, FUN = function(x){
              res  <- try(var_ccov_decompose(x$A,x$Z,time_points=50))
            eigens <- t(abs(res$eigenvals))
      congruencies <- res$subsequent_pair_congruencies
    return(list(eigens = eigens, congruencies = congruencies)) }); gc()

eigen_dat <- data.frame(
  Re(do.call(rbind,
             lapply(eigen_congurency, 
                    FUN = function(x) cbind( x$eigens, 1:51 ) ))))

upper <- as.matrix(eigen_dat %>% group_by(X4) %>% reframe( across(paste0( "X", 1:(length(eigen_dat)-1) ), ~ quantile(.x, c(.975))) ))
lower <- as.matrix(eigen_dat %>% group_by(X4) %>% reframe( across(paste0( "X", 1:(length(eigen_dat)-1) ), ~ quantile(.x, c(.025))) ))

# Figure ----
tiff(filename = "Figure_4_moretimepoints.tiff", 
     width    = 17, 
     height   = 19, 
     units    = "cm", 
     res      = 300,
     pointsize = 10)
par(mfrow     = c(2,2) )
par(mar       = c(4,4,2,0.5) )
colvec <- cividis(ncol(upper)-1) |> adjustcolor(alpha.f = 0.15)
matplot(t(result_parallel$eigenvals),
        type = "n",
        ylab = "Absolute value of eigenvalue", 
        xlab = expression(paste("Increment in time ", Delta, "t")),
        xaxt = "n",
        main = "Cross-covariance eigenvalues",
        font.main = 1); grid()
axis(side = 1, at = 1:(ncol(result_parallel$eigenvals)), labels = 0:(ncol(result_parallel$eigenvals) - 1))
for( i in 2:ncol(upper)) {
  polygon(x = c(upper[,1], rev(lower[,1])), y = c(upper[,i], rev(lower[,i])),
          col = colvec[i-1], border = F)
}
matplot(t(result_parallel$eigenvals), 
        type = "b",
        col  = cividis(ncol(upper)-1), add = T)


cong_dat <- data.frame(Re(do.call(rbind,lapply(eigen_congurency, FUN = function(x) cbind( x$congruencies, 1:50 ) ))))
upper_c  <- as.matrix(cong_dat %>% group_by(X2) %>% 
                        reframe( quantile(X1, 0.975) ))
lower_c  <- as.matrix(cong_dat %>% group_by(X2) %>% 
                        reframe( quantile(X1, 0.025) ))

matplot(result_parallel$subsequent_pair_congruencies, type = "n",
        ylim = c(0,1),
        ylab = "Congruency coefficient", 
        xlab = "Cross-covariance pair",
        xaxt = "n",
        main = "Largest eigenvector congruency",
        font.main = 1); grid()
polygon(x = c(upper_c[,1], rev(upper_c[,1])), y=c(upper_c[,2], rev(lower_c[,2])),
        col = colvec[1])
axis(1, labels = paste0("(", 0:50,", ", 1:51,")"),
     at = 1:51, cex.axis = 0.7 )
matplot(result_parallel$subsequent_pair_congruencies, type = "b",
        col = cividis(6), add = T )
qgraph( A, layout = "circle", 
        labels = varLabs, mar = c(2,2,7,2))
title("Coefficient matrix",
      font.main = 1,
      line     = -1)
qgraph( Z, layout = "circle", mar = c(3,3,7,3))
title("Innovation covariance",
      font.main = 1,
      line     = -1)

dev.off();gc();par(mfrow = c(1,1) )
