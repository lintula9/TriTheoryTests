# Black-and-white versions of all figures.
#
# Plots from the workspaces saved by
#   "Supplementary Appendix R -script reproduce figures without sampling.R"
# in ./Datas/workspace_*.RData, so nothing is recomputed and no MCMC output is read.
# Run from the repository root. Output: <figure>_BW.tiff next to the colour figures.
#
# Compared with the colour figures only the colouring changes: everything is black,
# negative edges in the networks are dashed instead of orange, credible bands are grey.

library(qgraph)

ws <- function(name) {                      # load one saved workspace into its own environment
  e <- new.env()
  load(file.path("Datas", paste0("workspace_", name, ".RData")), envir = e)
  e
}

band <- adjustcolor("black", alpha.f = 0.15)

bw_qgraph <- function(mat, labels, vsize = 18, mar = c(7, 7, 14, 7), ...) {
  qgraph(mat,
         layout      = "circle",
         labels      = labels,
         posCol      = "black",
         negCol      = "black",
         negDashed   = TRUE,      # line style, not colour, marks negative edges
         esize       = 12,
         vsize       = vsize,
         edge.width  = 1,
         diag        = TRUE,      # always draw the diagonal (innovation variances)
         label.scale = FALSE,
         label.cex   = 1,
         label.font  = 1,
         mar         = mar,
         ...)
}

# Singular values over time increments, optionally with 95% credible bands.
bw_singular <- function(res, main, ylab = "Singular value", upper = NULL, lower = NULL, ...) {
  n <- ncol(res$singularvals)               # increments 0 .. n-1
  matplot(t(res$singularvals), type = "n", ylab = ylab, main = main, font.main = 1,
          xlab = expression(paste("Increment in time ", Delta, "t")), xaxt = "n"); grid()
  axis(1, at = 1:n, labels = 0:(n - 1))
  if (!is.null(upper))
    for (i in 2:ncol(upper))
      polygon(c(upper[, 1], rev(lower[, 1])), c(upper[, i], rev(lower[, i])), col = band, border = NA)
  matplot(t(abs(res$singularvals)), type = "b", col = "black", add = TRUE, ...)
}

# Congruency of the largest eigenvector between subsequent cross-covariances, optionally with a band.
bw_congruency <- function(res, main, ylab = "Congruency coefficient", upper_c = NULL, lower_c = NULL) {
  m <- length(res$subsequent_pair_congruencies)   # pairs (0,1) .. (m-1,m)
  matplot(res$subsequent_pair_congruencies, type = "n", ylim = c(0, 1), ylab = ylab, main = main,
          font.main = 1, xlab = "Cross-covariance pair", xaxt = "n"); grid()
  axis(1, at = 1:m, labels = paste0("(", 0:(m - 1), ", ", 1:m, ")"), cex.axis = 0.7)
  if (!is.null(upper_c))
    polygon(c(upper_c[, 1], rev(upper_c[, 1])), c(upper_c[, 2], rev(lower_c[, 2])), col = band, border = NA)
  matplot(res$subsequent_pair_congruencies, type = "b", col = "black", add = TRUE)
}

# The four-panel empirical figure: singular values, congruencies, coefficient matrix, innovation covariance.
bw_figure4 <- function(file, res, upper, lower, upper_c, lower_c, A, Z, labels) {
  tiff(file, width = 17, height = 19, units = "cm", res = 300, pointsize = 10)
  par(mfrow = c(2, 2), mar = c(4, 4, 2, 0.5))
  bw_singular(res, "Cross-covariance singular values", upper = upper, lower = lower)
  bw_congruency(res, "Largest eigenvector congruency", upper_c = upper_c, lower_c = lower_c)
  bw_qgraph(A, labels); title("Coefficient matrix",    font.main = 1, line = -1)
  bw_qgraph(Z, labels); title("Innovation covariance", font.main = 1, line = -1)
  dev.off()
}

# Figure 1: VAR(1) fits to four samples from an indistinguishable model. ----
e <- ws("numerical_examples")
tiff("Figure_1_BW.tiff", width = 17, height = 19, units = "cm", res = 300, pointsize = 10)
par(mfrow = c(2, 2), mar = c(4, 4, 2, 0.5))
for (i in 1:4) bw_qgraph(e$As[[i]], expression(X[1], X[2], X[3], X[4]), mar = c(7, 7, 7, 7))
dev.off()

# Example for presentation: time-varying, indistinguishable models. ----
tiff("Example_for_pres_BW.tiff", height = 24, width = 16, res = 120, units = "in", pointsize = 24)
par(mfrow = c(3, 2))
pres <- function(mat, main)
  bw_qgraph(mat, e$labels, vsize = 10, mar = c(4, 4, 6, 4), maximum = e$max_weight, title = main, title.cex = 1.5)
pres(e$A_t_adj[,,1],  "Model 1 and 2, time point 0"); plot.new()
pres(e$A_t_adj[,,6],  "Model 1, time point 5");        pres(e$A_t[,,6],  "Model 2, time point 5")
pres(e$A_t_adj[,,10], "Model 1, time point 9");        pres(e$A_t[,,10], "Model 2, time point 9")
dev.off()

# Figure 2: distinguishable vs. perfectly indistinguishable VAR(1). ----
e <- ws("var_ccov_decompose")
tiff("Figure_2_tentative_BW.tiff", width = 17, height = 19, units = "cm", res = 300, pointsize = 10)
par(mfrow = c(2, 2), mar = c(4, 4, 2, 0.5))
bw_singular(e$parallel_A,   "Distinguishable cross-covariance", lty = 1)
bw_singular(e$parallel_B,   "Perfectly indistinguishable cross-covariance", ylab = "", lty = 1)
bw_congruency(e$parallel_A, "Unstable factor loadings")
bw_congruency(e$parallel_B, "Perfectly stable factor loadings", ylab = "")
dev.off()

# Figure 4: empirical example, 3 and 7 symptoms. ----
e <- ws("empiric_example")
bw_figure4("Figure_4_BW.tiff",           e$result_parallel,   e$upper,   e$lower,   e$upper_c,   e$lower_c,   e$A,   e$Z,   e$varLabs)
bw_figure4("Figure_4_tentative_BW.tiff", e$result_parallel_7, e$upper_7, e$lower_7, e$upper_c_7, e$lower_c_7, e$A_7, e$Z_7, e$varLabs2)

# Figure 4 with 50 time increments. ----
e <- ws("sensitivity_Figures")
bw_figure4("Figure_4_moretimepoints_BW.tiff",           e$result_parallel,   e$upper,   e$lower,   e$upper_c,   e$lower_c,   e$A,   e$Z,   e$varLabs)
bw_figure4("Figure_4_tentative_more_timepoints_BW.tiff", e$result_parallel_7, e$upper_7, e$lower_7, e$upper_c_7, e$lower_c_7, e$A_7, e$Z_7, e$varLabs2)

# Figure 4, sensitivity analysis (alternative Stan model). ----
e <- ws("empiric_example_sensitivity_analysis_more_time_points")
bw_figure4("Figure_4a_sensitivity_BW.tiff", e$result_parallel,   e$upper,   e$lower,   e$upper_c,   e$lower_c,   e$A,   e$Z,   e$varLabs)
bw_figure4("Figure_4_sensitivity_BW.tiff",  e$result_parallel_7, e$upper_7, e$lower_7, e$upper_c_7, e$lower_c_7, e$A_7, e$Z_7, e$varLabs2)
