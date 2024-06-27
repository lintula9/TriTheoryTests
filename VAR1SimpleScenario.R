source("Libraries.R")

# 3 Variable scenarios of longitudinal data as VAR(1) process. 

# 1.
# Assume between person-homogeneous, time-homogeneous, VAR(1) process, closed system (no external shocks).
# A: Simulate 100 observations from an VAR process ->
# B: Choose randomly timepoint1 of observations from VAR processes ->
# B: Create timepoint2 = timepoint1 + 100, equidistant intervals between observations.
# C, Result: Observe covariance of time1 (n = 100) and time 2 (n = 200), and calculate difference.

Nruns <- 10000
VARparams <- matrix(c(.4,.3,.3, 
                      .3,.3,.05,
                      .0,.05,.3), ncol = 3) 

simstime1 <- matrix(NA , ncol = 3, nrow = Nruns ) # Placeholder matrices for simulation data.
simstime2 <- matrix(NA , ncol = 3, nrow = Nruns )




if(F) {
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
    
if(F) {utils::write.csv(file = "Simulation1.csv" , x = cbind(simstime1, simstime2), )} # Save simulation
if(F) {simstime1 <- utils::read.csv("Simulation1.csv")[,2:4]
       simstime2 <- utils::read.csv("Simulation1.csv")[,5:7]} # Read simulation.
  
round( cov( na.omit(simstime1) ) - cov( na.omit(simstime2) ), 3 ) # C, Result: Difference between covariance matrices is small to non-existent, indicating 'measurement invariance'.

if(F){
# DAG presentation of data generating model: 
qgraph(VARparams, directed = T, diag = T, edge.labels = T)

# Example figure of sampling:

matplot(Sim1, type = "l", ylab = "Within person symptom process", xlab = "Time", main = "Two samples are observed from one person.", family = "serif")
abline( v = time1 ) # Observation 1
abline( v = time2 ) # Observation 2

ggplot() + 
  geom_line(aes(x = 1:200, y = Sim1[,1], color = "Symptom 1")) +
  geom_line(aes(x = 1:200, y = Sim1[,2], color = "Symptom 2")) +
  geom_line(aes(x = 1:200, y = Sim1[,3], color = "Symptom 3")) +
  scale_color_manual(values = c("Symptom 1" = "darkorange", 
                                "Symptom 2" = "black", 
                                "Symptom 3" = "blue4"),
                     guide = guide_legend(title = "")) +
  theme_minimal() +
  labs(title = "Three variable Network", x = "Time", y = "Within person symptom process")

  }
}



# Real life data, GAD - 7 example

library(openxlsx)
VARpcor <- read.xlsx("VARpcor.xlsx")
pcor = sapply(VARpcor[,-1], as.numeric)
rownames(pcor) <- VARpcor[,1]

# Smooth pcor with cor smoot.
library(psych)
pcor <- cor.smooth(pcor)

# Model
model <- '
  # Autoregressive paths
  GAD701_2 ~ GAD701_1 + GAD702_1 + GAD703_1 + GAD704_1 + GAD705_1 + GAD706_1 + GAD707_1
  GAD702_2 ~ GAD702_1 + GAD701_1 + GAD703_1 + GAD704_1 + GAD705_1 + GAD706_1 + GAD707_1
  GAD703_2 ~ GAD703_1 + GAD701_1 + GAD702_1 + GAD704_1 + GAD705_1 + GAD706_1 + GAD707_1
  GAD704_2 ~ GAD704_1 + GAD701_1 + GAD702_1 + GAD703_1 + GAD705_1 + GAD706_1 + GAD707_1
  GAD705_2 ~ GAD705_1 + GAD701_1 + GAD702_1 + GAD703_1 + GAD704_1 + GAD706_1 + GAD707_1
  GAD706_2 ~ GAD706_1 + GAD701_1 + GAD702_1 + GAD703_1 + GAD704_1 + GAD705_1 + GAD707_1
  GAD707_2 ~ GAD707_1 + GAD701_1 + GAD702_1 + GAD703_1 + GAD704_1 + GAD705_1 + GAD706_1
  
  GAD701_3 ~ GAD701_2 + GAD702_2 + GAD703_2 + GAD704_2 + GAD705_2 + GAD706_2 + GAD707_2
  GAD702_3 ~ GAD702_2 + GAD701_2 + GAD703_2 + GAD704_2 + GAD705_2 + GAD706_2 + GAD707_2
  GAD703_3 ~ GAD703_2 + GAD701_2 + GAD702_2 + GAD704_2 + GAD705_2 + GAD706_2 + GAD707_2
  GAD704_3 ~ GAD704_2 + GAD701_2 + GAD702_2 + GAD703_2 + GAD705_2 + GAD706_2 + GAD707_2
  GAD705_3 ~ GAD705_2 + GAD701_2 + GAD702_2 + GAD703_2 + GAD704_2 + GAD706_2 + GAD707_2
  GAD706_3 ~ GAD706_2 + GAD701_2 + GAD702_2 + GAD703_2 + GAD704_2 + GAD705_2 + GAD707_2
  GAD707_3 ~ GAD707_2 + GAD701_2 + GAD702_2 + GAD703_2 + GAD704_2 + GAD705_2 + GAD706_2

  GAD701_4 ~ GAD701_3 + GAD702_3 + GAD703_3 + GAD704_3 + GAD705_3 + GAD706_3 + GAD707_3
  GAD702_4 ~ GAD702_3 + GAD701_3 + GAD703_3 + GAD704_3 + GAD705_3 + GAD706_3 + GAD707_3
  GAD703_4 ~ GAD703_3 + GAD701_3 + GAD702_3 + GAD704_3 + GAD705_3 + GAD706_3 + GAD707_3
  GAD704_4 ~ GAD704_3 + GAD701_3 + GAD702_3 + GAD703_3 + GAD705_3 + GAD706_3 + GAD707_3
  GAD705_4 ~ GAD705_3 + GAD701_3 + GAD702_3 + GAD703_3 + GAD704_3 + GAD706_3 + GAD707_3
  GAD706_4 ~ GAD706_3 + GAD701_3 + GAD702_3 + GAD703_3 + GAD704_3 + GAD705_3 + GAD707_3
  GAD707_4 ~ GAD707_3 + GAD701_3 + GAD702_3 + GAD703_3 + GAD704_3 + GAD705_3 + GAD706_3

  GAD701_5 ~ GAD701_4 + GAD702_4 + GAD703_4 + GAD704_4 + GAD705_4 + GAD706_4 + GAD707_4
  GAD702_5 ~ GAD702_4 + GAD701_4 + GAD703_4 + GAD704_4 + GAD705_4 + GAD706_4 + GAD707_4
  GAD703_5 ~ GAD703_4 + GAD701_4 + GAD702_4 + GAD704_4 + GAD705_4 + GAD706_4 + GAD707_4
  GAD704_5 ~ GAD704_4 + GAD701_4 + GAD702_4 + GAD703_4 + GAD705_4 + GAD706_4 + GAD707_4
  GAD705_5 ~ GAD705_4 + GAD701_4 + GAD702_4 + GAD703_4 + GAD704_4 + GAD706_4 + GAD707_4
  GAD706_5 ~ GAD706_4 + GAD701_4 + GAD702_4 + GAD703_4 + GAD704_4 + GAD705_4 + GAD707_4
  GAD707_5 ~ GAD707_4 + GAD701_4 + GAD702_4 + GAD703_4 + GAD704_4 + GAD705_4 + GAD706_4

  GAD701_6 ~ GAD701_5 + GAD702_5 + GAD703_5 + GAD704_5 + GAD705_5 + GAD706_5 + GAD707_5
  GAD702_6 ~ GAD702_5 + GAD701_5 + GAD703_5 + GAD704_5 + GAD705_5 + GAD706_5 + GAD707_5
  GAD703_6 ~ GAD703_5 + GAD701_5 + GAD702_5 + GAD704_5 + GAD705_5 + GAD706_5 + GAD707_5
  GAD704_6 ~ GAD704_5 + GAD701_5 + GAD702_5 + GAD703_5 + GAD705_5 + GAD706_5 + GAD707_5
  GAD705_6 ~ GAD705_5 + GAD701_5 + GAD702_5 + GAD703_5 + GAD704_5 + GAD706_5 + GAD707_5
  GAD706_6 ~ GAD706_5 + GAD701_5 + GAD702_5 + GAD703_5 + GAD704_5 + GAD705_5 + GAD707_5
  GAD707_6 ~ GAD707_5 + GAD701_5 + GAD702_5 + GAD703_5 + GAD704_5 + GAD705_5 + GAD706_5

  GAD701_7 ~ GAD701_6 + GAD702_6 + GAD703_6 + GAD704_6 + GAD705_6 + GAD706_6 + GAD707_6
  GAD702_7 ~ GAD702_6 + GAD701_6 + GAD703_6 + GAD704_6 + GAD705_6 + GAD706_6 + GAD707_6
  GAD703_7 ~ GAD703_6 + GAD701_6 + GAD702_6 + GAD704_6 + GAD705_6 + GAD706_6 + GAD707_6
  GAD704_7 ~ GAD704_6 + GAD701_6 + GAD702_6 + GAD703_6 + GAD705_6 + GAD706_6 + GAD707_6
  GAD705_7 ~ GAD705_6 + GAD701_6 + GAD702_6 + GAD703_6 + GAD704_6 + GAD706_6 + GAD707_6
  GAD706_7 ~ GAD706_6 + GAD701_6 + GAD702_6 + GAD703_6 + GAD704_6 + GAD705_6 + GAD707_6
  GAD707_7 ~ GAD707_6 + GAD701_6 + GAD702_6 + GAD703_6 + GAD704_6 + GAD705_6 + GAD706_6

  GAD701_8 ~ GAD701_7 + GAD702_7 + GAD703_7 + GAD704_7 + GAD705_7 + GAD706_7 + GAD707_7
  GAD702_8 ~ GAD702_7 + GAD701_7 + GAD703_7 + GAD704_7 + GAD705_7 + GAD706_7 + GAD707_7
  GAD703_8 ~ GAD703_7 + GAD701_7 + GAD702_7 + GAD704_7 + GAD705_7 + GAD706_7 + GAD707_7
  GAD704_8 ~ GAD704_7 + GAD701_7 + GAD702_7 + GAD703_7 + GAD705_7 + GAD706_7 + GAD707_7
  GAD705_8 ~ GAD705_7 + GAD701_7 + GAD702_7 + GAD703_7 + GAD704_7 + GAD706_7 + GAD707_7
  GAD706_8 ~ GAD706_7 + GAD701_7 + GAD702_7 + GAD703_7 + GAD704_7 + GAD705_7 + GAD707_7
  GAD707_8 ~ GAD707_7 + GAD701_7 + GAD702_7 + GAD703_7 + GAD704_7 + GAD705_7 + GAD706_7

  GAD701_9 ~ GAD701_8 + GAD702_8 + GAD703_8 + GAD704_8 + GAD705_8 + GAD706_8 + GAD707_8
  GAD702_9 ~ GAD702_8 + GAD701_8 + GAD703_8 + GAD704_8 + GAD705_8 + GAD706_8 + GAD707_8
  GAD703_9 ~ GAD703_8 + GAD701_8 + GAD702_8 + GAD704_8 + GAD705_8 + GAD706_8 + GAD707_8
  GAD704_9 ~ GAD704_8 + GAD701_8 + GAD702_8 + GAD703_8 + GAD705_8 + GAD706_8 + GAD707_8
  GAD705_9 ~ GAD705_8 + GAD701_8 + GAD702_8 + GAD703_8 + GAD704_8 + GAD706_8 + GAD707_8
  GAD706_9 ~ GAD706_8 + GAD701_8 + GAD702_8 + GAD703_8 + GAD704_8 + GAD705_8 + GAD707_8
  GAD707_9 ~ GAD707_8 + GAD701_8 + GAD702_8 + GAD703_8 + GAD704_8 + GAD705_8 + GAD706_8

  GAD701_10 ~ GAD701_9 + GAD702_9 + GAD703_9 + GAD704_9 + GAD705_9 + GAD706_9 + GAD707_9
  GAD702_10 ~ GAD702_9 + GAD701_9 + GAD703_9 + GAD704_9 + GAD705_9 + GAD706_9 + GAD707_9
  GAD703_10 ~ GAD703_9 + GAD701_9 + GAD702_9 + GAD704_9 + GAD705_9 + GAD706_9 + GAD707_9
  GAD704_10 ~ GAD704_9 + GAD701_9 + GAD702_9 + GAD703_9 + GAD705_9 + GAD706_9 + GAD707_9
  GAD705_10 ~ GAD705_9 + GAD701_9 + GAD702_9 + GAD703_9 + GAD704_9 + GAD706_9 + GAD707_9
  GAD706_10 ~ GAD706_9 + GAD701_9 + GAD702_9 + GAD703_9 + GAD704_9 + GAD705_9 + GAD707_9
  GAD707_10 ~ GAD707_9 + GAD701_9 + GAD702_9 + GAD703_9 + GAD704_9 + GAD705_9 + GAD706_9

  GAD701_11 ~ GAD701_10 + GAD702_10 + GAD703_10 + GAD704_10 + GAD705_10 + GAD706_10 + GAD707_10
  GAD702_11 ~ GAD702_10 + GAD701_10 + GAD703_10 + GAD704_10 + GAD705_10 + GAD706_10 + GAD707_10
  GAD703_11 ~ GAD703_10 + GAD701_10 + GAD702_10 + GAD704_10 + GAD705_10 + GAD706_10 + GAD707_10
  GAD704_11 ~ GAD704_10 + GAD701_10 + GAD702_10 + GAD703_10 + GAD705_10 + GAD706_10 + GAD707_10
  GAD705_11 ~ GAD705_10 + GAD701_10 + GAD702_10 + GAD703_10 + GAD704_10 + GAD706_10 + GAD707_10
  GAD706_11 ~ GAD706_10 + GAD701_10 + GAD702_10 + GAD703_10 + GAD704_10 + GAD705_10 + GAD707_10
  GAD707_11 ~ GAD707_10 + GAD701_10 + GAD702_10 + GAD703_10 + GAD704_10 + GAD705_10 + GAD706_10

  GAD701_12 ~ GAD701_11 + GAD702_11 + GAD703_11 + GAD704_11 + GAD705_11 + GAD706_11 + GAD707_11
  GAD702_12 ~ GAD702_11 + GAD701_11 + GAD703_11 + GAD704_11 + GAD705_11 + GAD706_11 + GAD707_11
  GAD703_12 ~ GAD703_11 + GAD701_11 + GAD702_11 + GAD704_11 + GAD705_11 + GAD706_11 + GAD707_11
  GAD704_12 ~ GAD704_11 + GAD701_11 + GAD702_11 + GAD703_11 + GAD705_11 + GAD706_11 + GAD707_11
  GAD705_12 ~ GAD705_11 + GAD701_11 + GAD702_11 + GAD703_11 + GAD704_11 + GAD706_11 + GAD707_11
  GAD706_12 ~ GAD706_11 + GAD701_11 + GAD702_11 + GAD703_11 + GAD704_11 + GAD705_11 + GAD707_11
  GAD707_12 ~ GAD707_11 + GAD701_11 + GAD702_11 + GAD703_11 + GAD704_11 + GAD705_11 + GAD706_11

  GAD701_13 ~ GAD701_12 + GAD702_12 + GAD703_12 + GAD704_12 + GAD705_12 + GAD706_12 + GAD707_12
  GAD702_13 ~ GAD702_12 + GAD701_12 + GAD703_12 + GAD704_12 + GAD705_12 + GAD706_12 + GAD707_12
  GAD703_13 ~ GAD703_12 + GAD701_12 + GAD702_12 + GAD704_12 + GAD705_12 + GAD706_12 + GAD707_12
  GAD704_13 ~ GAD704_12 + GAD701_12 + GAD702_12 + GAD703_12 + GAD705_12 + GAD706_12 + GAD707_12
  GAD705_13 ~ GAD705_12 + GAD701_12 + GAD702_12 + GAD703_12 + GAD704_12 + GAD706_12 + GAD707_12
  GAD706_13 ~ GAD706_12 + GAD701_12 + GAD702_12 + GAD703_12 + GAD704_12 + GAD705_12 + GAD707_12
  GAD707_13 ~ GAD707_12 + GAD701_12 + GAD702_12 + GAD703_12 + GAD704_12 + GAD705_12 + GAD706_12
'

# Fit the model using the polychoric correlation matrix
library(lavaan)
fit <- lavaan::sem(model, sample.cov = pcor, sample.nobs = 2178, estimator = "ML")

# Summarize the results
summary(fit, fit.measures = TRUE, standardized = TRUE)

# Plot the lavaan fit.
library(semPlot)
semplot <- semPaths(fit, "std", "est")
         
                    