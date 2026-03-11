source("./second publication/continuous_decompose.R")
# Tammilehto et al., 2025. ----
# http://journals.sagepub.com/doi/suppl/10.1177/21677026241301057
# The most complex models are used - as they were best fitting.

# Clinical cohort model: Depression iCBT.
Tammilehto_model_4_A <- matrix(
  c(-.324,  .055  , 
     .043, -.244 ), 2, 2, byrow = TRUE)
Tammilehto_model_4_Q <- matrix(
  c(.352^2,         .352*.337*.418, 
    .352*.337*.418, .337^2), 2, 2, byrow = TRUE)

# Clinical cohort models: GAD iCBT.
Tammilehto_model_8_A <- matrix(
  c(-.191,  .034  , 
     .044, -.137 ), 2, 2, byrow = TRUE)
Tammilehto_model_8_Q <- matrix(
  c(.233^2,         .233*.292*.467, 
    .233*.292*.467, .292^2), 2, 2, byrow = TRUE)

# Clinical cohort models: Pooled iCBT.
Tammilehto_model_12_A <- matrix(
  c(-.297,  .041, 
    -.016, -.164 ), 2, 2, byrow = TRUE)
Tammilehto_model_12_Q <- matrix(
  c(.304^2,         .304*.318*.504, 
    .304*.318*.504, .318^2), 2, 2, byrow = TRUE)

# Combine Tammilehto et al., 2025 esimates. ----
tammilehto_2025 <- lapply(
  paste0("Tammilehto_model_", 1:12, "_"),
  FUN = \(model_name) {
    A <- get(paste0(model_name, "A"))
    Q <- get(paste0(model_name, "Q"))
    return(continuous_decompose(A,Q))
      })
