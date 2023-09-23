# VAR (1) process simulations, simple scenario:
set.seed(123)
VARparams <- matrix(c(.3,.3,.3,
                      .3,.3,.3,
                      .3,.3,.3), ncol = 3) 


Sim3 <- simulateVAR(VARparams,  
            means = 0, 
            lags = 1, Nt = 100, init = c(0,0,0), residuals = 0.1)
matplot(Sim3, type = "l", ylab = "Symptom level", xlab = "Time", ylim = c(-1.5,1.5))
describe(Sim3)

# Latent variable AR (1) process simple scenario:
Sim4 <- simulateVAR(matrix(.5),
                    means = 0, 
                    lags = 1, Nt = 100, init = 0, residuals = 0.1)
Sim4obs <- Sim4[,1] * .7
plot(x = 1:100, y = Sim4obs, type = "l")
