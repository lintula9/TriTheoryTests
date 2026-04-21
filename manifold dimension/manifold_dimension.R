#' Correlation integral example.
#' We simulate data from a 2-d manifold with non-linear charts to 3-d ambient space.
#' Then, using the correlation dimension, we obtain the ordinary dimension.

library(ggplot2)
library(plotly)

# Generate synthetic manifold data: 2D Swiss roll embedded in 3D
set.seed(42)
n_samples <- 5000
t <- runif(n_samples, 0, 4*pi)
h <- runif(n_samples, 0, 10)
x <- t * cos(t)
y <- h
z <- t * sin(t)
data <- cbind(x, y, z)

# Compute pairwise distances
distances <- as.matrix(dist(data, method='euclidean'))

# Correlation integral for various epsilon values
epsilon_values <- 10^seq(-1, 1, length.out=30)
correlation_sum <- numeric(length(epsilon_values))

for (i in seq_along(epsilon_values)) {
  epsilon <- epsilon_values[i]
  C_eps <- sum(distances < epsilon) / (nrow(data) * (nrow(data) - 1))
  correlation_sum[i] <- C_eps
}

# Estimate correlation dimension using linear regression on log-log plot
valid_idx <- (correlation_sum > 0) & is.finite(log(correlation_sum))
log_epsilon <- log(epsilon_values[valid_idx])
log_corr <- log(correlation_sum[valid_idx])

# Fit line: slope = correlation dimension
fit <- lm(log_corr ~ log_epsilon)
D_c <- coef(fit)[2]

cat(sprintf("Estimated Correlation Dimension: %.3f\n", D_c))
cat("Expected Intrinsic Dimension: ~2.0 (Swiss roll)\n")

# Visualization - Plot 1: Log-log plot
plot1_data <- data.frame(log_epsilon, log_corr)
fitted_line <- data.frame(
  log_epsilon = log_epsilon,
  fitted = predict(fit)
)

plot1 <- ggplot(plot1_data, aes(x=log_epsilon, y=log_corr)) +
  geom_point(alpha=0.6) +
  geom_line(data=fitted_line, aes(y=fitted), color='red', linewidth=1) +
  labs(x='log(ε)', y='log(C(ε))', 
       title='Correlation Dimension Estimation',
       subtitle=sprintf('D_c = %.3f', D_c)) +
  theme_minimal() +
  theme(panel.grid=element_line(color='gray90'))

print(plot1)

# Visualization - Plot 2: 3D scatter of Swiss roll
df_3d <- data.frame(x=x, y=y, z=z, t=t)
plot2 <- plot_ly(df_3d, x=~x, y=~y, z=~z, color=~t, type='scatter3d', mode='markers',
                  marker=list(size=3, opacity=0.6)) %>%
  layout(scene=list(xaxis=list(title='X'),
                    yaxis=list(title='Y'),
                    zaxis=list(title='Z')),
         title='Swiss Roll Manifold (2D embedded in 3D)')

print(plot2)