source("./second publication/continuous_decompose.R")
source("./Libraries.R")
# Tammilehto et al., 2025. ----
# http://journals.sagepub.com/doi/suppl/10.1177/21677026241301057
# The most complex models are used - as they were best fitting.

# Clinical cohort model: Depression iCBT.
Tammilehto_2025_depiCBT_A <- matrix(
  c(-.3245,  .0552  , 
     .0426, -.2437 ), 2, 2, byrow = TRUE)
Tammilehto_2025_depiCBT_Q <- matrix(
  c( .1245,  .0823, 
     .0823,  .1139), 2, 2, byrow = TRUE)

# Clinical cohort models: GAD iCBT.
Tammilehto_2025_gadiCBT_A <- matrix(
  c(-.1906,  .0336  , 
     .0336, -.1366 ), 2, 2, byrow = TRUE)
Tammilehto_2025_gadiCBT_Q <- matrix(
  c(.1245, .0823,
    .0823, .1139), 2, 2, byrow = TRUE)

# Clinical cohort models: Pooled iCBT.
Tammilehto_2025_pooliCBT_A <- matrix(
  c(-.2961,  .0392, 
     .0166, -.1641 ), 2, 2, byrow = TRUE)
Tammilehto_2025_pooliCBT_Q <- matrix(
  c(.0920, .0764,
    .0764, .1015), 2, 2, byrow = TRUE)

## Combine Tammilehto et al., 2025 estimates. ----
tammilehto_2025 <- lapply(
  grep("Tammilehto_2025_(.*)_*", ls(), value = T,) |> 
    gsub(pattern = "(_A|_Q)", replacement = "_") |> 
    unique(),
  FUN = \(model_name) {
    A <- get(paste0(model_name, "A"))
    Q <- get(paste0(model_name, "Q"))
    return(decompose_stationary_cov(A,Q))
      })

# Tammilehto et al., 2026. ----
# Worry <-> somatic symptoms main model.
Tammilehto_2026_main_A <- matrix(
  c(-.44,  -.22, 
     .29,  -1.14 ), 2, 2, byrow = TRUE)
Tammilehto_2026_main_Q <- matrix(
  c(.3011, .2692,
    .2692, .4366), 2, 2, byrow = TRUE)

## Combine Tammilehto et al., 2026 estimates. ----
tammilehto_2026 <- lapply(
  grep("Tammilehto_2026_(.*)_*", ls(), value = T,) |> 
    gsub(pattern = "(_A|_Q)", replacement = "_") |> 
    unique(),
  FUN = \(model_name) {
    A <- get(paste0(model_name, "A"))
    Q <- get(paste0(model_name, "Q"))
    return(decompose_stationary_cov(A,Q,"continuous"))
  })

# Fried et al., 2022. ----
Fried_2022_detrended_A <- read.csv("./Datas/fried_2022_A.csv") |>
  select(-X) |> as.matrix()
Fried_2022_detrended_Q <- read.csv("./Datas/fried_2022_Q.csv") |>
  select(-X) |> as.matrix()
# Select psychopathology variables.
decompose_stationary_cov(Fried_2022_detrended_A,
                         Fried_2022_detrended_Q,
                         type = "discrete")
# Lintula et al., 2026. ----
A_7 <- readRDS("C:/Users/lintu/Documents/TriTheoryTests/Datas/A_7.RDS")
Z_7 <- readRDS("C:/Users/lintu/Documents/TriTheoryTests/Datas/Z_7.RDS")
decompose_stationary_cov(A_7,
                         Z_7,
                         type = "discrete")
