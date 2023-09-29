# AR(1), latent autoregression.


# With independent error terms (not in VAR models..), covariance structure is simple.
lvar <- arima.sim(model = list( .3, 0, 0 ), n = 200)
X <- ( lvar %*% t( c( .7, .7, .7) ) ) + matrix(ncol = 3, rnorm(600, sd = .1))
cov(X)
