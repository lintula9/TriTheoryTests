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
# 7 Variables. ----
## Read data ----
inference_vars_regex_alpha <- c("A_effective","A","psi","Lambda", "Omega","L_", "cutpoints", "time_of_day_effect",
                                "ref_time_of_day_effect", "specific_time_of_day_effect")
nchains = 8
K = 7
# Read data
if(F){
    ## CHECK THAT THIS IS THE CORRECT FILE
  fit_Net_7    <- as_cmdstan_fit(files = paste0("Datas/BVAR_7_variables_01_09-202609021852-",1:8,"-032488.csv") ); gc()
  draws_data_7 <- as_draws_df(fit_Net_7$draws(variables = c(inference_vars_regex_alpha, "lp__") )); gc(); rm(fit_Net_7); gc()
}
# Diagnostics and posterior distribution marignal plots. 
plotnams <- inference_vars_regex_alpha;pdf(file = paste0("Datas/Bayespots_7VAR_",format(Sys.time(), "%Y-%m-%d"), ".pdf"));color_scheme_set("viridis")
for(i in plotnams){print(  mcmc_trace(draws_data_7, regex_pars = i ));print(  mcmc_areas_ridges(draws_data_7[ , grep(i, names(draws_data_7))]) );
}; gc(); dev.off()


# Extract posterior means as the VAR parameters:
A_7 <- matrix( unlist(colMeans(draws_data_7[,grep("A", names(draws_data_7))])),     ncol = K, nrow = K);
Z_7 <- matrix( unlist(colMeans(draws_data_7[,grep("Omega", names(draws_data_7))])), ncol = K, nrow = K)

# Obtain complete VAR(1) model samples:
estimated_var_samples_7 <- draws_data_7[,c( grep("A", names(draws_data_7)) , grep("Omega", names(draws_data_7)) )]; As <- grep("A_", names(estimated_var_samples_7)); Os <- grep("Omega", names(estimated_var_samples_7)); 
var_samples_7           <- pbapply::pblapply(1:nrow(draws_data_7),
                                           FUN = function(i){
                                             A_temp <- matrix( unlist(estimated_var_samples_7[i,As]), ncol = sqrt(length(As)), nrow = sqrt(length(Os)));
                                             Z_temp <- matrix( unlist(estimated_var_samples_7[i,Os]), ncol = sqrt(length(As)), nrow = sqrt(length(Os)));
                                             return(list(A = A_temp, Z = Z_temp))
                                           }); gc()
# Source the methods.
source("./first publication/Supplementary Appendix R -script var_ccov_decompose.R")

# Figure 5 in main text ----
# Compute parallel analysis imitation and RMSEA, for the posterior mean.
result_parallel_7  <- var_ccov_decompose(A_7, Z_7, time_points = 50)

# Compute credible intervals for eigenvalues, congruencies.
eigen_congurency_7 <- pbapply::pblapply(var_samples_7, FUN = function(x){
  res  <- try(var_ccov_decompose(x$A,x$Z, time_points = 50))
  eigens       <- t(abs(res$eigenvals))
  singularvals <- t(res$singularvals)
  congruencies <- res$subsequent_pair_congruencies
  return(list(eigens = eigens, singularvals = singularvals, congruencies = congruencies)) }); gc()

# Compute quantiles for eigenvalues.
eigen_distributions_7 <- pbapply::pblapply(var_samples_7, FUN = function(x){
  res          <- var_ccov_decompose(x$A,x$Z,time_points=50)
  eigens       <- t((res$eigenvals))
  singularvals <- t((res$singularvals))
  return(list(eigens = eigens)) }); gc()

eigen_dat_7    <- data.frame(Re(do.call(rbind,lapply(eigen_congurency_7, 
                                                FUN = function(x) cbind( x$eigens, 1:51 ) ))))
singular_dat_7 <- data.frame(Re(do.call(rbind,lapply(eigen_congurency_7, 
                                                     FUN = function(x) cbind( x$singularvals, 1:51 ) ))))

upper_7 <- as.matrix(singular_dat_7 %>% group_by(X8) %>% reframe( across(paste0( "X", 1:(length(singular_dat_7)-1) ), ~ quantile(.x, c(.975))) ))
lower_7 <- as.matrix(singular_dat_7 %>% group_by(X8) %>% reframe( across(paste0( "X", 1:(length(singular_dat_7)-1) ), ~ quantile(.x, c(.025))) ))

tiff(filename = "Figure_4_tentative_more_timepoints.tiff", 
     width    = 17, 
     height   = 19, 
     units    = "cm", 
     res      = 300,
     pointsize = 10)
par(mfrow     = c(2,2) )
par(mar       = c(4,4,2,0.5) )
colvec <- cividis(ncol(upper_7)-1) |> adjustcolor(alpha.f = 0.15)
matplot(t(result_parallel_7$singularvals),
        col  = colvec,
        type = "n",
        ylab = "Singular value", 
        xlab = expression(paste("Increment in time ", Delta, "t")),
        xaxt = "n",
        main = "Cross-covariance singular values",
        font.main = 1); grid()
axis(side = 1, at = 1:(ncol(result_parallel_7$singularvals)), 
     labels = 0:(ncol(result_parallel_7$singularvals) - 1))

matplot(t(result_parallel_7$singularvals), 
        type = "b",
        col  = cividis(ncol(upper_7)-1), add = T)
for( i in 2:ncol(upper_7)) {
  polygon(x = c(upper_7[,1], rev(lower_7[,1])), y = c(upper_7[,i], rev(lower_7[,i])),
          col = colvec[i-1], 
          border = NA)
}

cong_dat_7 <- data.frame(Re(do.call(rbind,lapply(eigen_congurency_7, FUN = function(x) cbind( x$congruencies, 1:51 ) ))))
upper_c_7  <- as.matrix(cong_dat_7 %>% group_by(X2) %>% 
                        reframe( quantile(X1, 0.975) ))
lower_c_7  <- as.matrix(cong_dat_7 %>% group_by(X2) %>% 
                        reframe( quantile(X1, 0.025) ))
matplot(result_parallel_7$subsequent_pair_congruencies, type = "n",
        ylim = c(0,1),
        ylab = "Congruency coefficient", 
        xlab = "Cross-covariance pair",
        xaxt = "n",
        main = "Largest eigenvector congruency",
        font.main = 1); grid()
polygon(x = c(upper_c_7[,1], rev(upper_c_7[,1])), y=c(upper_c_7[,2], rev(lower_c_7[,2])),
        col = adjustcolor(cividis(1), alpha.f = 0.15), border = NA )
axis(1, labels = paste0("(", 0:50,", ", 1:51,")"),
     at = 1:51, cex.axis = 0.7 )
matplot(result_parallel_7$subsequent_pair_congruencies, type = "b",
        col = cividis(6), add = T )
qgraph( A_7, layout = "circle", 
        labels = varLabs2, mar = c(2,2,7,2) )
title("Coefficient matrix",
      font.main = 1,
      line     = -1)
qgraph( Z_7, layout = "circle", 
        labels = varLabs2, mar = c(3,3,7,3))
title("Innovation covariance",
      font.main = 1,
      line     = -1)
dev.off();gc();par(mfrow = c(1,1) )
