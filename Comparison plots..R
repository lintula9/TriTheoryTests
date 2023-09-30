# Comparison plots.
# Models do not always have samve innovation variance.

sapply( dir()[ grep("(.R)(1)", dir()) ], source) # Takes a mom

par(mfrow = c(2,2))
matplot(Sim1, type = "l")
matplot(Sim2, type = "l")
matplot(Sim3, type = "l")
matplot(Sim4, type = "l")
