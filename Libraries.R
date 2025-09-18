# Libraries

packages <- c("lavaan", "qgraph", "psych", "mlVAR", "knitr", 
              "semTools", "semPlot","ggplot2", "fastmatrix")

for ( i in packages ) {
  if ( !requireNamespace( i, 
                          quietly = T )) {
    install.packages( i )
  }
  library( i, character.only = TRUE )

  }
