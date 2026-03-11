source("./second publication/continuous_decompose.R")
# Tammilehto et al., 2025. ----
# http://journals.sagepub.com/doi/suppl/10.1177/21677026241301057

# Clinical cohort models: Depression iCBT.
## Model 1. ----
Tammilehto_model_1_A <- matrix(
  c(-0.179,  0.081 , 
     0.098, -0.152 ), 2, 2, byrow = TRUE)
Tammilehto_model_1_Q <- matrix(
  c(0.481^2,             0.481*0.463*0.152, 
    0.481*0.463*0.152, 0.463^2), 2, 2, byrow = TRUE)

## Model 2. ----
Tammilehto_model_2_A <- matrix(
  c(-.572,  .286   , 
     .371, -.643 ), 2, 2, byrow = TRUE)
Tammilehto_model_2_Q <- matrix(
  c(.559^2,           .559*.602*.034, 
    .559*.602*.034, .602^2), 2, 2, byrow = TRUE)

# Model 3. ----
Tammilehto_model_3_A <- matrix(
  c(-.246,  .021  , 
     .082, -.260 ), 2, 2, byrow = TRUE)
Tammilehto_model_3_Q <- matrix(
  c(.321^2,         .321*.335*.448,
    .321*.335*.448, .335^2), 2, 2, byrow = TRUE)

# Model 4. ----
Tammilehto_model_4_A <- matrix(
  c(-.324,  .055  , 
     .043, -.244 ), 2, 2, byrow = TRUE)
Tammilehto_model_4_Q <- matrix(
  c(.352^2,         .352*.337*.418, 
    .352*.337*.418, .337^2), 2, 2, byrow = TRUE)

# Clinical cohort models: GAD iCBT.
## Model 5. ----
Tammilehto_model_5_A <- matrix(
  c( -.066,  .032 , 
      .040, -.072 ), 2, 2, byrow = TRUE)
Tammilehto_model_5_Q <- matrix(
  c(.299^2,         .299*.351*.172, 
    .299*.351*.172, .351^2), 2, 2, byrow = TRUE)

## Model 6. ----
Tammilehto_model_6_A <- matrix(
  c(-.349,  .188   , 
     .173, -.190 ), 2, 2, byrow = TRUE)
Tammilehto_model_6_Q <- matrix(
  c(.410^2,         .410*.402*.038, 
    .410*.402*.038, .402^2), 2, 2, byrow = TRUE)

# Model 7. ----
Tammilehto_model_7_A <- matrix(
  c(-.224,  .118  , 
     .048, -.080 ), 2, 2, byrow = TRUE)
Tammilehto_model_7_Q <- matrix(
  c(.250^2,         .250*.254*.440,
    .250*.254*.440, .254^2), 2, 2, byrow = TRUE)

# Model 8. ----
Tammilehto_model_8_A <- matrix(
  c(-.191,  .034  , 
     .044, -.137 ), 2, 2, byrow = TRUE)
Tammilehto_model_8_Q <- matrix(
  c(.233^2,         .233*.292*.467, 
    .233*.292*.467, .292^2), 2, 2, byrow = TRUE)

# Clinical cohort models: Pooled iCBT.
## Model 9. ----
Tammilehto_model_9_A <- matrix(
  c( -.136,  .046 , 
      .062, -.100 ), 2, 2, byrow = TRUE)
Tammilehto_model_9_Q <- matrix(
  c(.424^2,         .424*.423*.178, 
    .424*.423*.178, .423^2), 2, 2, byrow = TRUE)
## Model 10. ----
Tammilehto_model_10_A <- matrix(
  c(-.620,  .300   , 
     .320, -.520 ), 2, 2, byrow = TRUE)
Tammilehto_model_10_Q <- matrix(
  c(.548^2,         .548*.566*.030, 
    .548*.566*.030, .566^2), 2, 2, byrow = TRUE)
# Model 11. ----
Tammilehto_model_11_A <- matrix(
  c(-.294,  .098  , 
     .026, -.109 ), 2, 2, byrow = TRUE)
Tammilehto_model_11_Q <- matrix(
  c(.288^2,         .288*.299*.544,
    .288*.299*.544, .299^2), 2, 2, byrow = TRUE)
# Model 12. ----
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
