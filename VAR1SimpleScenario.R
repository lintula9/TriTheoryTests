source("Libraries.R")

# Assume person-homogeneous, time-homogeneous, VAR(1) process, closed system (no external shocks).
# A: Simulate 100 observations from an VAR process ->
# B: Choose randomly timepoint1 of observations from VAR processes ->
# B: Create timepoint2 = timepoint1 + 100, equidistant intervals between observations.
# C: Observe covariance of time1 (n = 100) and time 2 (n = 200), and calculate difference.

VARparams <- matrix(c(.3,.3,.3, # VAR process is the same for every observation.
                      .3,.3,.3,
                      .3,.3,.3), ncol = 3) 

simstime1 <- matrix(NA , ncol = 3, nrow = 100 ) # Placeholder matrices for simulation data.
simstime2 <- matrix(NA , ncol = 3, nrow = 100 )


for (i in 1:100) {
  simui <- simulateVAR( # A
            VARparams,  
            means = 0, 
            lags = 1, 
            Nt = 200, 
            init = c(0,0,0), 
            residuals = 0.1)

  time1 <- sample(1:100, 1) # B
  time2 <- time1 + 100 # B

  simstime1[ i , ] <- as.matrix(simui[ time1, ])
  simstime2[ i , ] <- as.matrix(simui[ time2, ])
}

round(cov(simstime1) - cov(simstime2), 3) # C: Difference between covariance matrices is small, indicating 'measurement invariance'.

# Example figure:

matplot(simui, type = "l", ylab = "Within person symptom process", xlab = "Time", main = "Two samples are observed from one person.", family = "serif")
abline( v = time1 ) # Observation 1
abline( v = time2 ) # Observation 2

