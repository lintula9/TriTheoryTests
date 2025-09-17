# cmdstanr options
pkgs <- c("cmdstanr", "posterior", "bayesplot")
for( i in pkgs ) {
  ifelse(
  !requireNamespace(i), 
  yes = \(i) {install.packages(i); library(i, character.only = T)},
  no  = library(i, character.only = T)
  )}
options("cmdstanr_output_dir" = "./Datas/")
color_scheme_set(scheme = "viridis")
