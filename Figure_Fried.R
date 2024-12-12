# Figure Fried.
source("Fried_2022_reanalysis.R")

  # Figure.
par(mfrow=c(2,3))

# Maximum for scaling
A_max = max(A, A_tilde)
Z_max = max(Z, Z_tilde)

# Define edge properties for Z
qgraph(Z, 
       layout = "circle", 
       labels = expression(X[1], X[2], X[3], X[4]),
       maximum = Z_max) # Optionally add edge labels
title("VAR(1) and indistinguishable VAR(1).", outer = TRUE)

# Define edge properties for Z_tilde
qgraph(Z_tilde, 
       layout = "circle", 
       labels = expression(X[1], X[2], X[3], X[4]),
       maximum = Z_max)

# Define edge properties for A
qgraph(abs(Z - Z_tilde), 
       layout = "circle", 
       labels = expression(X[1], X[2], X[3], X[4]),
       maximum = Z_max)


# Define edge properties for A_tilde
qgraph(A_tilde, 
       layout = "circle", 
       labels = expression(X[1], X[2], X[3], X[4]),
       maximum = A_max)

# Define edge properties for A
qgraph(A, 
       layout = "circle", 
       labels = expression(X[1], X[2], X[3], X[4]),
       maximum = A_max)

# Diofferences
qgraph(abs(A_tilde - A), 
       layout = "circle", 
       labels = expression(X[1], X[2], X[3], X[4]),
       maximum = A_max)



par(mfrow=c(1,1))