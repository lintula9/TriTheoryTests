# Libraries
packages <- c("lavaan", "qgraph", "psych", "mlVAR", "knitr", 
              "semTools", "semPlot","ggplot2", "fastmatrix",
              "MTS","ggplot2","gridExtra","rstan","progress","tidyr","tibble",
              "dplyr","foreach","doSNOW",
              "parallel") |> unique()
for ( i in packages ) {
  if ( !requireNamespace( i, quietly = T )) {
    install.packages( i )}
  library( i, character.only = TRUE )}
