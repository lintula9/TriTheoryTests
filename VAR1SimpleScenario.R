source("Libraries.R")

# 3 Variable scenarios of longitudinal data as VAR(1) process. 

# 1.
# Assume between person-homogeneous, time-homogeneous, VAR(1) process, closed system (no external shocks).
# A: Simulate 100 observations from an VAR process ->
# B: Choose randomly timepoint1 of observations from VAR processes ->
# B: Create timepoint2 = timepoint1 + 100, equidistant intervals between observations.
# C, Result: Observe covariance of time1 (n = 100) and time 2 (n = 200), and calculate difference.


VARparams <- matrix(c(.3,.3,.3, 
                      .3,.3,.3,
                      .3,.3,.3), ncol = 3) 

simstime1 <- matrix(NA , ncol = 3, nrow = 100 ) # Placeholder matrices for simulation data.
simstime2 <- matrix(NA , ncol = 3, nrow = 100 )

if (FALSE) {
for ( i in 1:100 ) { 
  Sim1 <- simulateVAR( # A
            VARparams,  
            means = 0, 
            lags = 1, 
            Nt = 200, 
            init = c( 0,0,0 ), 
            residuals = 0.1)

  time1 <- sample( 1:100, 1) # B
  time2 <- time1 + 100 # B

  simstime1[ i , ] <- as.matrix(Sim1[ time1, ])
  simstime2[ i , ] <- as.matrix(Sim1[ time2, ])
}

round( cov( simstime1 ) - cov( simstime2 ), 3 ) # C, Result: Difference between covariance matrices is small to non-existent, indicating 'measurement invariance'.

# DAG presentation of data generating model: 
qgraph(VARparams, directed = T, diag = T, edge.labels = T)

# Example figure of sampling:

matplot(Sim1, type = "l", ylab = "Within person symptom process", xlab = "Time", main = "Two samples are observed from one person.", family = "serif")
abline( v = time1 ) # Observation 1
abline( v = time2 ) # Observation 2
}
