# Cross-sectional basic setting that violates 1-dimensionality assumption.

source("Libraries.R")

Model1 <- "
X1 ~~ .6*X2 
X1 ~~ .6*X3
X2 ~~ .6*X3
X3 ~~ .6*X4

# 'False' correlations are freely estimated.
X1 ~~ X4 
X2 ~~ X4
"

Simulation1 <- simulateData( model = Model1, sample.nobs = 10000 )

FARes1 <- efa( Simulation1, std.lv = T )

par( mfrow = c( 2, 1 ))
qgraph( FARes1$loadings %*% t( FARes1$loadings ) , diag = F ) 
qgraph( cor( Simulation1 ) , diag = F )
par( mfrow = c( 1 , 1 ) )
