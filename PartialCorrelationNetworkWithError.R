# Partial correlation with measurement error test - GGM:

library(MASS)

Niters <- 100000
X <- mvrnorm( n = Niters, 
         Sigma = matrix(c(2,.7,.7,
                          .7,2,.7,
                          .7,.7,2), ncol = 3),
         mu = c(0,0,0))
TrueSigma <- matrix(c(2,.7,.7,
                      .7,2,.7,
                      .7,.7,2), ncol = 3)
solve( cov2cor(TrueSigma) )

solve(cor( X )) ; 
solve(cor( round(X, digits = 1 ) )) ; solve(cor( round(X, digits = 0 ) )) # Rounding error is bad.

X_withError <- X + ( matrix( runif(n = 3 * Niters, min = -2, 2 ), ncol = 3 ) )
solve(cor( X_withError )) ; 

describe(X)
describe(X_withError)

# How about polychoric correlatins ?

X_rounded <- round( X, digits = 0 )
X_rounded[ X_rounded < 0 ] <- 0
X_rounded[ X_rounded > 5 ] <- 5


X_withErrorRounded <- round( X_withError, digits = 0 )
X_withErrorRounded[ X_withErrorRounded < 0 ] <- 0
X_withErrorRounded[ X_withErrorRounded > 5 ] <- 5

solve( polychoric( X_rounded )$rho ); 
solve( polychoric( X_withErrorRounded )$rho )

describe(X_withErrorRounded)

factorResult <- sem(model = 'F1 =~ X1+X2+X3',
                    data = data.frame(X_withErrorRounded), ordered = T, estimator = "WLSMV" )
( inspect(factorResult, "std")$lambda %*% t(inspect(factorResult, "std")$lambda) ) + inspect(factorResult, "std")$theta
lavInspect( factorResult, what = "implied")
lavInspect( factorResult, what = "information")

polychoric(X_withErrorRounded)
