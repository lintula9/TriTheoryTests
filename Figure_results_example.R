# Figure 1.
library(Matrix); library(qgraph)

lambda = c(1,2,3)
A = lambda %*% t(lambda) * c(t(lambda) %*% lambda)^-1 + 
  matrix(c(0,3,-2,
           -2,1,0,
           1,-0.5,0), ncol = 3, nrow = 3, byrow = T)
Z = (1-0.5^2) * lambda %*% t(lambda)

# Plot.
labels <- expression(X[1], X[2], X[3])
par(mfrow=c(1,2))
coeflabs = expression(tilde(A), tilde(A)[adj.])
qgraph(lambda %*% t(lambda) * c(t(lambda) %*% lambda)^-1,
       mar = c(4, 4, 6, 4),
       layout = "circle", 
       labels = labels, # Use the expression labels
       directed = T,
       title = coeflabs[1],
       title.cex = 1.5
       )
qgraph(A,
       mar = c(4, 4, 6, 4),
       layout = "circle", 
       labels = labels, # Use the expression labels
       title = coeflabs[2],
       title.cex = 1.5
       )
par(mfrow=c(1,1))
dev.off()