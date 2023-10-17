# Libraries

packages <- c("lavaan", "qgraph", "psych", "mlVAR", "knitr")

for ( i in packages ) {
  if ( !requireNamespace( i, 
                          quietly = T )) {
    install.packages( i )
  }
  library( i, character.only = TRUE )

  }
