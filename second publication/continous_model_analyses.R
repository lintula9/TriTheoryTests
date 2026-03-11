source("./second publication/continuous_decompose.R")
# Tammilehto et al., 2025. ----
# http://journals.sagepub.com/doi/suppl/10.1177/21677026241301057
# The most complex models are used - as they were best fitting.

# Ask for parameter estimates from Jaakko.

# Clinical cohort model: Depression iCBT.
Tammilehto_2025_depiCBT_A <- matrix(
  c(-.3250,  .0558  , 
     .0431, -.2442 ), 2, 2, byrow = TRUE)
Tammilehto_2025_depiCBT_Q <- matrix(
  c( .2244,  .1881, 
     .1881,  .2664), 2, 2, byrow = TRUE)
continuous_decompose(Tammilehto_2025_depiCBT_A,
                     Tammilehto_2025_depiCBT_Q)

# Clinical cohort models: GAD iCBT.
Tammilehto_2025_gadiCBT_A <- matrix(
  c(-.191,  .034  , 
     .044, -.137 ), 2, 2, byrow = TRUE)
Tammilehto_2025_gadiCBT_Q <- matrix(
  c(.233^2,         .233*.292*.467, 
    .233*.292*.467, .292^2), 2, 2, byrow = TRUE)

# Clinical cohort models: Pooled iCBT.
Tammilehto_2025_pooliCBT_A <- matrix(
  c(-4.2833,  .4641, 
      .1359, -3.4495 ), 2, 2, byrow = TRUE)
Tammilehto_2025_pooliCBT_Q <- matrix(
  c(.2953, .0619,
    .0619, .1612), 2, 2, byrow = TRUE)

## Combine Tammilehto et al., 2025 estimates. ----
tammilehto_2025 <- lapply(
  grep("Tammilehto_2025_(.*)_*", ls(), value = T,) |> 
    gsub(pattern = "(_A|_Q)", replacement = "_") |> 
    unique(),
  FUN = \(model_name) {
    A <- get(paste0(model_name, "A"))
    Q <- get(paste0(model_name, "Q"))
    return(continuous_decompose(A,Q))
      })

# Tammilehto et al., 2026. ----
