# =============================================================================
# 3-D Ornstein-Uhlenbeck process with a 1-D approximation ("Lambda")
#
#   dX_t = -Theta (X_t - mu) dt + Sigma dW_t ,   Theta = theta * I
#
# Lambda is the first principal axis of the sampled path: the line through the
# sample mean along the leading eigenvector of the path's covariance, i.e. the
# line minimising total squared distance to the sampled points.
# =============================================================================

set.seed(42)

# ---------------------------------------------------------------------------
# 1. Parameters
# ---------------------------------------------------------------------------
theta <- 1.2                 # isotropic mean-reversion rate
Theta <- theta * diag(3)     # drift matrix = theta * I
mu    <- c(0, 0, 0)          # long-run mean

# Theta is isotropic as requested, so the anisotropy lives in the diffusion:
# with Sigma isotropic too the stationary law would be spherical and PC1 pure
# sampling noise. Set Sigma <- diag(3) to see that degenerate case.
D <- diag(c(1.60, 0.55, 0.30))    # noise scale along three orthogonal axes

# Rotate those axes so the elongation is not aligned with x/y/z
rot <- function(a, b, g) {
  Rz <- matrix(c(cos(a), -sin(a), 0, sin(a), cos(a), 0, 0, 0, 1), 3, 3, byrow = TRUE)
  Ry <- matrix(c(cos(b), 0, sin(b), 0, 1, 0, -sin(b), 0, cos(b)), 3, 3, byrow = TRUE)
  Rx <- matrix(c(1, 0, 0, 0, cos(g), -sin(g), 0, sin(g), cos(g)), 3, 3, byrow = TRUE)
  Rz %*% Ry %*% Rx
}
R0    <- rot(pi / 5, pi / 7, pi / 9)
Sigma <- R0 %*% D                 # diffusion matrix

# ---------------------------------------------------------------------------
# 2. Exact simulation
# ---------------------------------------------------------------------------
# For Theta = theta*I the transition law is closed form, so there is no
# Euler-Maruyama discretisation error:
#
#   X_{t+dt} | X_t  ~  N( mu + e^{-theta dt}(X_t - mu),  V )
#   V = (1 - e^{-2 theta dt}) / (2 theta) * Sigma Sigma'

n  <- 1500        # number of steps
dt <- 0.02        # time step

SS   <- Sigma %*% t(Sigma)
a    <- exp(-theta * dt)
V    <- (1 - exp(-2 * theta * dt)) / (2 * theta) * SS
Lv   <- t(chol(V))                 # lower-triangular factor, V = Lv Lv'
Vinf <- SS / (2 * theta)           # stationary covariance
Linf <- t(chol(Vinf))

X      <- matrix(0, nrow = n, ncol = 3, dimnames = list(NULL, c("x", "y", "z")))
X[1, ] <- mu + Linf %*% rnorm(3)   # start from the stationary distribution
for (i in 2:n) {
  X[i, ] <- mu + a * (X[i - 1, ] - mu) + Lv %*% rnorm(3)
}

# ---------------------------------------------------------------------------
# 3. Lambda: first principal component of the sampled path
# ---------------------------------------------------------------------------
pca <- prcomp(X, center = TRUE, scale. = FALSE)
ctr <- pca$center                  # sample mean
v1  <- pca$rotation[, 1]           # leading eigenvector (unit length)
s   <- pca$x[, 1]                  # scores = (X - ctr) %*% v1

# Sign convention: make v1 point towards +x for reproducible labelling
if (v1[1] < 0) { v1 <- -v1; s <- -s }

pad    <- 0.03 * diff(range(s))
Lambda <- rbind(ctr + (min(s) - pad) * v1,
                ctr + (max(s) + pad) * v1)
colnames(Lambda) <- c("x", "y", "z")

var_expl <- pca$sdev^2 / sum(pca$sdev^2)
cat(sprintf("Lambda direction : (%.3f, %.3f, %.3f)\n", v1[1], v1[2], v1[3]))
cat(sprintf("Through point    : (%.3f, %.3f, %.3f)\n", ctr[1], ctr[2], ctr[3]))
cat(sprintf("Variance explained by PC1/2/3: %.1f%% / %.1f%% / %.1f%%\n",
            100 * var_expl[1], 100 * var_expl[2], 100 * var_expl[3]))

# ---------------------------------------------------------------------------
# 4. Colours
# ---------------------------------------------------------------------------
# grDevices::palette.colors(); first entry of "Okabe-Ito" is black, which is
# reserved for Lambda, so take the first non-black colour for the path.
pal        <- palette.colors(palette = "Okabe-Ito")
is_black   <- toupper(substr(pal, 1, 7)) %in% c("#000000")
path_col   <- unname(pal[!is_black][1])          # "#E69F00" (orange)
lambda_col <- "black"
path_alpha <- 1

# ---------------------------------------------------------------------------
# 5. Plot (plotly)
# ---------------------------------------------------------------------------
library(plotly)

# toRGB() turns the hex into "rgba(230,159,0,0.6)"; plotly.js parses that,
# whereas an 8-digit hex from adjustcolor() is not reliably understood.
path_rgba <- toRGB(path_col, alpha = path_alpha)

Xdf <- as.data.frame(X[1:1000,])
Ldf <- as.data.frame(Lambda)
lab <- "\u1D6B2"                    # the character Lambda
p <- plot_ly() |>
  add_trace(data = Xdf, x = ~x, y = ~y, z = ~z,
            type = "scatter3d", mode = "lines",
            line = list(color = path_rgba, width = 2),
            hoverinfo = "skip") |>
  add_trace(data = Ldf, x = ~x, y = ~y, z = ~z,
            type = "scatter3d", mode = "lines",
            line = list(color = lambda_col, width = 8),
            hoverinfo = "skip") |>
  add_trace(x = Lambda[2, 1], y = Lambda[2, 2], z = Lambda[2, 3],
            type = "scatter3d", mode = "text",
            text = lab, textfont = list(size = 28, color = lambda_col),
            hoverinfo = "skip") |>
  layout(title = "",
         showlegend = FALSE,
         scene = list(xaxis = list(title = "", showticklabels = FALSE,
                                   visible = FALSE),
                      yaxis = list(title = "", showticklabels = FALSE,
                                   visible = FALSE),
                      zaxis = list(title = "", showticklabels = FALSE,
                                   visible = FALSE),
                      aspectmode = "data"))

print(p)
