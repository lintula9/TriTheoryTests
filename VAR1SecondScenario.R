# 2.
# Same as 1., but change:
# Between person homogeneous -> heterogeneous.
# D: Create random VAR(1) process coefficient (drift) matrix:


simstime1 <- matrix(NA , ncol = 3, nrow = 100 ) # Placeholder matrices for simulation data.
simstime2 <- matrix(NA , ncol = 3, nrow = 100 )

for ( i in 1:100 ) { 
  VARparams <- matrix(c(.3,.3,.3, 
                        .3,.3,.3,
                        .3,.3,.3) + 
                        runif( 9 , -.1 , .1 ), # D: Add idiosyncratic randomness, ~ Uni( -.1, .1 ), to VAR coefficient matrix.
                      ncol = 3) 
  
  simui <- simulateVAR( # A
    VARparams,  
    means = 0, 
    lags = 1, 
    Nt = 200, 
    init = c( 0,0,0 ), 
    residuals = 0.1 )
  
  time1 <- sample( 1:100, 1) # B
  time2 <- time1 + 100 # B
  
  simstime1[ i , ] <- as.matrix( simui[ time1, ] )
  simstime2[ i , ] <- as.matrix( simui[ time2, ] )
}

round( cov( simstime1 ) - cov( simstime2 ), 2 ) # C: Result: Difference between covariance matrices is clear.

# Examples of DAG presentations of data generating models: 
par(mfrow = c(2,2))
sapply(1:4, 
       function (x) qgraph(matrix(c(.3,.3,.3, 
                                    .3,.3,.3,
                                    .3,.3,.3) + 
                                    runif( 9 , -.1 , .1 ), # D: Add idiosyncratic randomness, ~ Uni( -.1, .1 ), to VAR coefficient matrix.
                                  ncol = 3) , directed = T, diag = T, edge.labels = T, edge.color = "black", layout = "circle", edge.label.cex = 2))
par(mfrow = c(1,1))


# Examples of different VAR processes:
VARexamples <- array( dim = c( 200, 3, 4 ))
ExampleVARCoefs <- array( dim = c( 3, 3, 4 ))
ExampleVARCoefs[ , , 1] <- matrix( c(.3,.3,.3, 
                                     .3,.3,.3,
                                     .3,.3,.3), ncol = 3 )
ExampleVARCoefs[ , , 2] <- matrix( c(.3,.3,.05, 
                                     .3,.3,.05,
                                     .05,.05,.3), ncol = 3 )
ExampleVARCoefs[ , , 3] <- matrix( c(.05,.3,.3, 
                                     .3,.05,.3,
                                     .3,.3,.05), ncol = 3 )
ExampleVARCoefs[ , , 4] <- matrix( c(.3,.05,.05, 
                                     .3,.3,.05,
                                     .3,.05,.3), ncol = 3 )

for ( i in 1:4 ) {
  
  VARexamples[ , , i ] <- as.matrix( simulateVAR( 
    ExampleVARCoefs[ , , i ],  
    means = 0, 
    lags = 1, 
    Nt = 200, 
    init = c( 0, 0, 0 ), 
    residuals = 0.1 ) ) }


par( mfrow = c( 2, 1 ))
sapply( 1 : 4 , function( i ) {
  
  qgraph(ExampleVARCoefs[ , , i ], edge.labels = T, diag = T, mar = c( 12, 12, 10, 10 ) )
  matplot(VARexamples[ , , i ], type = "l", ylab = "Within person symptom process", xlab = "Time", main = "", family = "serif")
  time1 <- sample( 1 : 100, 1 )
  abline( v = time1 ) 
  abline( v = time1 + 100 ) 
  
})
par( mfrow = c( 1, 1) )





