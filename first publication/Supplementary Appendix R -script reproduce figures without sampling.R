# Reproduce the figures of a supplementary script WITHOUT running the Stan samplers.
#
# Usage (from the repository root):
#   Rscript "first publication/Supplementary Appendix R -script reproduce figures without sampling.R" "<script.R>"
#   Rscript "first publication/Supplementary Appendix R -script reproduce figures without sampling.R" "<script.R>" replot
#
# With "replot" the saved workspace is loaded instead, and only the package loading and the
# tiff(...) ... dev.off() blocks of the script are evaluated: figures are re-drawn in seconds.
#
# The script is evaluated top to bottom, except that
#   - expressions calling $sample(), cmdstan_model() or save_output_files() are skipped,
#   - the bayesplot diagnostic PDF blocks (mcmc_trace / Bayespots) are skipped,
#   - the bodies of if(F){...} blocks that read saved draws (as_cmdstan_fit) or draw a tiff() are run,
#   - the remaining if(F){...} blocks (saveRDS / readRDS) are left untouched.
# At the end the workspace is saved to ./Datas/workspace_<script>.RData, which is what
# "Supplementary Appendix R -script black and white figures.R" plots from.

script <- commandArgs(trailingOnly = TRUE)[1]
replot <- identical(commandArgs(trailingOnly = TRUE)[2], "replot")
short  <- gsub(" ", "_", sub("[.]R$", "", sub("^Supplementary Appendix R -script ", "", basename(script))))
ws     <- file.path("Datas", paste0("workspace_", short, ".RData"))
if (replot) { load(ws, envir = globalenv()); cat("workspace loaded:", ws, "\n") }
in_fig <- FALSE
options(repos = c(CRAN = "https://cloud.r-project.org"), warn = 1)
suppressWarnings(pbapply::pboptions(type = "none"))
pdf(NULL)                      # sink for on-screen plots; tiff() opens its own device
t0 <- Sys.time()
exprs <- parse(file = script, keep.source = FALSE)
skip_re <- "[$]sample[(]|cmdstan_model[(]|save_output_files[(]|mcmc_trace|Bayespots"
is_if_false <- function(e) is.call(e) && identical(e[[1]], as.name("if")) &&
  (identical(e[[2]], as.name("F")) || isFALSE(e[[2]]))
for (i in seq_along(exprs)) {
  e   <- exprs[[i]]
  txt <- paste(deparse(e), collapse = "\n")
  head1 <- substr(strsplit(txt, "\n")[[1]][1], 1, 90)
  if (grepl(skip_re, txt)) { cat(sprintf("[%3d] SKIP  %s\n", i, head1)); next }
  if (is_if_false(e)) {
    body_txt <- paste(deparse(e[[3]]), collapse = "\n")
    if (grepl("as_cmdstan_fit|tiff[(]", body_txt)) { e <- e[[3]]; cat(sprintf("[%3d] UNGATE if(F)\n", i)) }
    else { cat(sprintf("[%3d] SKIP  if(F) save/read block\n", i)); next }
  }
  if (identical(e, quote(dev.off())) && dev.cur() == 1) { cat(sprintf("[%3d] SKIP  dev.off() (no device)\n", i)); next }
  if (replot) {                       # only package loading and figure blocks
    is_fig <- in_fig || grepl("tiff[(]", txt)
    if (!is_fig && !grepl("library[(]|require[(]", txt)) next
    in_fig <- is_fig && !grepl("dev[.]off[(][)]", txt)
  }
  cat(sprintf("[%3d] %s  %s\n", i, format(Sys.time(), "%H:%M:%S"), head1)); flush.console()
  eval(e, envir = globalenv())
}
cat(sprintf("DONE %s in %.1f min\n", basename(script), as.numeric(difftime(Sys.time(), t0, units = "mins"))))

# Keep the post-processed objects so figures can be re-drawn without re-reading the draws.
if (!replot) { save.image(ws); cat("workspace saved:", ws, "\n") }
