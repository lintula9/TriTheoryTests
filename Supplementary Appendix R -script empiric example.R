# Original code was created by Eiko Fried, Faidra Papanikolaou & Sascha Epskamp. 
# This version of code was (heavily) edited in 29.11.2024 for reanalysis of Fried(2022)
# Mental Health and Social Contact During the COVID-19 Pandemic: 
# An Ecological Momentary Assessment Study.
# For the original code and data, see https://osf.io/kj5zh

# --------------- 1. Loading packages & Data ------------------------------
# List of required packages
required_packages <- c(
  "Matrix", "fastmatrix", "BVAR", "expm", "qgraph", "tidyverse"
)

# Function to check and install missing packages
install_if_missing <- function(packages) {
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg)
    }
    library(pkg, character.only = TRUE)
  }
}

# Install and load packages
install_if_missing(required_packages)

# --------------- 2. Estimation -----------------------------


# load 
load(file.path(getwd(), "/Fried_2022 data/clean_network.RData")); gc()
Data5b <- Data2

# Alpha to detrend:
alpha <- 0.05

# Variables to investigate:
vars <- c("Q1","Q2","Q3","Q4","Q5","Q6","Q7","Q8","Q9",
          "Q10","Q11","Q12","Q13","Q14","Q15","Q16","Q17","Q18")

# Labels:
varLabs <- c("Relax","Irritable","Worry","Nervous","Future","Anhedonia",
             "Tired","Hungry","Alone","Angry","Social_offline","Social_online","Music",
             "Procrastinate","Outdoors","C19_occupied","C19_worry","Home")
names(Data5b)[names(Data5b) %in% vars] <- varLabs



# Select a plausible set of variables, which could be explained by a one CF model:
Data5b <- as_tibble(Data5b %>% select("Relax", "Worry", 
                            "Nervous", 
                            "id", "beep", "day", "conc"))
varLabs <- c("Relax", "Worry",  #Adjustment for below to work. 03.12.2024. Sakari Lintula
            "Nervous")


# Obtain the estimates:
# N = 1 estimation is done using ID 57, as it had a long timeseries.
res_bayes_N1 <- bvar(data = Data5b  %>%
                    filter(id == 57) %>% 
                    select(Relax, Worry, Nervous) %>% 
                    na.omit(),
               lags = 1L
               )

# Whole data, for comparison, but it is clustered hence not interpretable.
res_bayes_all <- bvar(data = Data5b  %>%
                       select(Relax, Worry, Nervous) %>% 
                       na.omit(),
                     lags = 1L
)



# --------------------- 3. Obtain the closest indistinguishable model ------

# Get the results: Between coefficient matrix as well as between innovation covariance.
A = coef(res_bayes_N1, type = "mean")[2:4,]
Z = vcov(res_bayes_N1, type = "mean")


# Find the closest indistinguishable VAR. Note, that there are many.
#Method 1
dcf_var <- civ_find(A,Z, tol = 1e-10, n.iter = 3000)
dcf_var$A_result - A
dcf_var$Z_result - Z
dcf_var$Loadings
#Method 2
dcf_var2 <- civ_find2(A,Z, tol = 1e-10, n.iter = 3000)
dcf_var2$A_result - A
dcf_var2$Z_result - Z
dcf_var2$Loadings

# Plot the result

# Plot.
max_weight <- max(A)

par(mfrow = c(2, 2), oma = c(0, 0, 4, 0)) # Adjust oma for outer margin to accommodate the title
labels = varLabs
qgraph(Z, 
       title = "'Contemporaneous' covariance",
       title.cex = 1.5,
       mar = c(4, 4, 6, 4),
       layout = "circle", 
       labels = labels # Use the expression labels
)
qgraph(A, 
       title = "Lagged effects\nTime point 1",
       title.cex = 1.5,  
       mar = c(4, 4, 6, 4),
       maximum = max_weight, # Consistent scale for edges
       edge.width = 3,       # Adjust this value for larger edges
       layout = "circle",
       labels = labels       # Use the expression labels
)
qgraph(dcf_var$A_result, 
       title = "\nTime point 5",
       title.cex = 1.5,  
       mar = c(4, 4, 6, 4),
       maximum = max_weight, # Consistent scale for edges
       edge.width = 3,       # Adjust this value for larger edges
       layout = "circle",
       labels = labels       # Use the expression labels
)
qgraph(dcf_var$Z_result, 
       title = "\nTime point 9",
       title.cex = 1.5,  
       mar = c(4, 4, 6, 4),
       maximum = max_weight, # Consistent scale for edges
       edge.width = 3,       # Adjust this value for larger edges
       layout = "circle",
       labels = labels       # Use the expression labels
)

mtext("VAR(1) Network model indistinguishable from a dynamic CF model.", outer = TRUE, cex = 1.5, font = 2)

par(mfrow = c(1,1))


# Simplest method of inspecting eigenvectors of A and Z:

which.min(eigen(Z)$vectors - eigen(A)$vectors ) 










# Scratch:
# Solve for the covariance
library(expm); library(Matrix)
Sigma_VAR = matrix(solve(diag(1, ncol = ncol(A)^2, nrow = nrow(A)^2) - fastmatrix::kronecker.prod(A)) %*% fastmatrix::vec(Z),
                   ncol = ncol(A), nrow = nrow(A))
K_VAR = function( Delta ) {
  
  A %^% Delta %*% Sigma_VAR
  }

eigen(Sigma_VAR) # Suggestive of there being one large component. 
# Capture Lambda as the normalized eigenvector of the largest eigenvalue.

L = eigen(Sigma_VAR)$vectors[,1] / c(sqrt(t(eigen(Sigma_VAR)$vectors[,1]) %*% eigen(Sigma_VAR)$vectors[,1]))
eigen(A)$vectors[,1] # The first vector of A is nearly in the same direction.
eigen(A)$vectors[,1] - (-1)*L 

  # Project A onto Lambda
A_0 =  L %*% t(L) %*% A %*% L %*% t(L)

  # Project the difference onto the orthogonal complement of Lambda
B_tilde = (A - A_0) %*% (diag(1, nrow = nrow(A), ncol = ncol(A)) - L %*% t(L))

# Create the VAR(1), indistinguishable from one dimensional D-CF(1) model.
A_tilde = A_0 + B_tilde
Z_tilde = L %*% t(L) %*% Z %*% L %*% t(L)# Project onto Lambda, as it must be proportional to LL^T


par(mfrow=c(2,2))
qgraph(Z, layout = "circle", labels = expression(X[1], X[2], X[3], X[4]))
title("VAR(1) and indistinguishable VAR(1).",outer = T)

qgraph(Z_tilde, layout = "circle", labels = expression(X[1], X[2], X[3], X[4]))
qgraph(A_tilde, layout = "circle", labels = expression(X[1], X[2], X[3], X[4]))
qgraph(A, layout = "circle", labels = expression(X[1], X[2], X[3], X[4]))
par(mfrow=c(1,1))

Sigma_DCF = matrix(solve(diag(1, ncol = ncol(A)^2, nrow = nrow(A)^2) - fastmatrix::kronecker.prod(A_tilde)) %*% fastmatrix::vec(Z_tilde),
                   ncol = ncol(A), nrow = nrow(A))
K_VAR = function( Delta ) {
  
  A_tilde %^% Delta %*% Sigma_DCF
}

