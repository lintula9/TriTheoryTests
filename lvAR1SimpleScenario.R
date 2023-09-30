# AR(1), latent autoregression.

# With independent error terms (not in VAR models..), covariance structure is simple.
n_lvar <- 200
lvar <- arima.sim(model = list( ar = .9 ), n = 200, rand.gen = function(n, ...) rnorm(n))
Sim3 <- ( lvar %*% t( c( .7, .7, .7) ) ) + matrix(ncol = 3, rnorm(3 * n_lvar, sd = 1))
describe( Sim3 ) ; cov( Sim3 ) ; plot( density( Sim3 ) ) ; cor( Sim3 )
matplot( Sim3, type = "l" )

c( .7, .7, .7) %*% t(c( .7, .7, .7) ) # Latent variable implied covariance matrix.




# With AR(1) process residuals:

Sim4 <- ( lvar %*% t( c( .7, .7, .7) ) ) + matrix(ncol = 3, rep(arima.sim(model = list( ar = .9 ), n = n_lvar), 3), byrow = F)
describe(Sim4) ; cov(Sim4) ; plot(density(Sim4)) ; cor(Sim4)
matplot(Sim4, type = "l")




