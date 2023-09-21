# Longitudinal setting with 2 measurement points and assumed AR(1) process with 0.2 cross-item regression coefficients, medium .5 AR(1).
# Total effect AR(1) process.

Model2 <- "
# Timepoint 1.
X1_t1 ~~ .5 * X2_t1
X1_t1 ~~ .5 * X3_t1
X2_t1 ~~ .5 * X3_t1

X1_t2 ~ .5 * X1_t1 + .2 * X2_t1 + .2 * X3_t1
X2_t2 ~ .2 * X1_t1 + .5 * X2_t1 + .2 * X3_t1
X3_t2 ~ .2 * X1_t1 + .2 * X2_t1 + .5 * X3_t1

"


Sim2 <- simulateData( model = Model2, sample.nobs = 10000 )
round(cor(Sim2),2)
# ARCoefMat <- matrix( c(.6,.2,.2,
#                        .2,.6,.2,
#                        .2,.2,.6), ncol = 3)
# Sim2[ , c( "X1_t2", "X2_t2", "X3_t2" )] = (as.matrix(Sim2) %*% ARCoefMat) # + error.

Sim2Lng <- reshape(Sim2, varying = list( c("X1_t1", "X1_t2") , c( "X2_t1", "X2_t2" ), c( "X3_t1", "X3_t2" ) ) , direction = "long")

# Fit 1-dimensional model with invariant loadings:

Model3 <- "
# 1-Dimension:
F1_t1 =~ b1 * X1_t1 + b2 * X2_t1 + b3 * X3_t1
F1_t2 =~ b1 * X1_t2 + b2 * X2_t2 + b3 * X3_t2

# Residuals are allowed to correlate:
X1_t1 ~~ X1_t2
X2_t1 ~~ X2_t2
X3_t1 ~~ X3_t2

# Factors are allowed to correlate regression:
F1_t1 ~~ F1_t2
"

FARes2 <- sem( data = Sim2, model = Model3 )
inspect( FARes2, "est" )


