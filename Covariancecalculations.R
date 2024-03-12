# Covariance numerical calulations
library(fastmatrix)
A = matrix(c(.25,.2,.2,.2,
             .2,.25,.2,.2,
             .2,.2,.25,.2,
             .2,.2,.2,.3), ncol = 4, byrow = F)
Psi = diag( 0.1, 
            ncol = ncol(A), 
            nrow = nrow(A) )
I = diag(1, ncol = ncol(A)*ncol(A), nrow = nrow(A)*nrow(A))

Sigma = matrix(solve(I - kronecker.prod(A)) %*% vec(Psi), ncol = 4, nrow = 4)
Sigma_12 = A %*% Sigma
Sigma_13 = A %*% Sigma_12
Sigma_14 = A %*% Sigma_13
Sigma_15 = A %*% Sigma_14

eigen(Sigma - Psi)

# The 2K x 2K

C2K = rbind(cbind(Sigma,Sigma_12),cbind(Sigma_12,Sigma))
colnames(C2K) <- c(paste("X",1:ncol(A),1, sep = ""),paste("X",1:ncol(A),2, sep = ""))
rownames(C2K) <- c(paste("X",1:ncol(A),1, sep = ""),paste("X",1:ncol(A),2, sep = ""))
LMI_model = "
CF_T1 =~ X11 + X21 + X31 + X41
CF_T2 =~ X12 + X22 + X32 + X42
CF_T2 ~ CF_T1
"
library(lavaan)
LMI = sem(model = LMI_model, sample.cov = C2K, sample.nobs = 10000)
summary(LMI)

# We obtain measurement invariance

# In fact we can observe that the regression coefficient of subsequent latent latent variables is 0.8 in this case.
C3K = do.call(args = list(
            cbind(Sigma,Sigma_12,Sigma_13),
            cbind(Sigma_12,Sigma,Sigma_12),
            cbind(Sigma_13,Sigma_12,Sigma)),
            what = rbind); print(round(C3K, 4))
# this is given by simple division observable from
Sigma_12 * 0.8 == Sigma_13
Sigma_13 * 0.8 == Sigma_14

# Manually fit a S-LMI
lambda = sqrt(0.144444444444444444)



# Fit s-lmi model
library(lavaan);library(semTools)
colnames(C2K) = c((paste("X",1:4,sep="")), paste("X", 5:8,sep=""))
rownames(C2K) = c((paste("X",1:4,sep="")), paste("X", 5:8,sep=""))

SLMI = "
  CF1 =~ X1 + X2 + X3 + X4
  CF2 =~ X5 + X6 + X7 + X8

  CF2 ~ CF1

  X1 ~~ X5
  X2 ~~ X6
  X3 ~~ X7
  X4 ~~ X8"
LMI_result = sem(model = SLMI, sample.cov = C2K, 
                 sample.nobs = 1000, 
                 std.lv = T, 
                 estimator = "ml", 
                 parameterization = "theta", auto.var = F)
lavInspect(LMI_result)


efa(sample.cov = C2K, rotation = "oblimin", sample.nobs = 1000, nfactors = 2, rotation.args = list(orthogonal = F))

longFacNames <- list(CF = c("CF1", "CF2"))
LMI_result= measEq.syntax(configural.model = configural.model, 
                          ID.fac = "std.lv", 
                          longFacNames = longFacNames, long.equal = "residual.covariances",
                          return.fit = T, sample.cov = C2K, sample.nobs = 10000)

