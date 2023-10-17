source("Libraries.R")

# 3 Variable scenarios of longitudinal data as VAR(1) process. 

# 1.
# Assume between person-homogeneous, time-homogeneous, VAR(1) process, closed system (no external shocks).
# A: Simulate 100 observations from an VAR process ->
# B: Choose randomly timepoint1 of observations from VAR processes ->
# B: Create timepoint2 = timepoint1 + 100, equidistant intervals between observations.
# C, Result: Observe covariance of time1 (n = 100) and time 2 (n = 200), and calculate difference.

Nruns <- 10000
VARparams <- matrix(c(.3,.3,.3, 
                      .3,.3,.3,
                      .3,.3,.3), ncol = 3) 

simstime1 <- matrix(NA , ncol = 3, nrow = Nruns ) # Placeholder matrices for simulation data.
simstime2 <- matrix(NA , ncol = 3, nrow = Nruns )




if (T) {
  pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                       max = Nruns, # Maximum value of the progress bar
                       style = 3,    # Progress bar style (also available style = 1 and style = 2)
                       width = 50,   # Progress bar width. Defaults to getOption("width")
                       char = "=")   # Character used to create the bar
  
  for (i in 1:Nruns) {
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
    
    setTxtProgressBar(pb, i)
    
  }
    

round( cov( na.omit(simstime1) ) - cov( na.omit(simstime2) ), 3 ) # C, Result: Difference between covariance matrices is small to non-existent, indicating 'measurement invariance'.

if(F){
# DAG presentation of data generating model: 
qgraph(VARparams, directed = T, diag = T, edge.labels = T)

# Example figure of sampling:

matplot(Sim1, type = "l", ylab = "Within person symptom process", xlab = "Time", main = "Two samples are observed from one person.", family = "serif")
abline( v = time1 ) # Observation 1
abline( v = time2 ) # Observation 2
}
}
