# Original code was created by Eiko Fried, Faidra Papanikolaou & Sascha Epskamp. 
# This version of code was (heavily) edited in 29.11.2024 for reanalysis of Fried(2022)
# Mental Health and Social Contact During the COVID-19 Pandemic: 
# An Ecological Momentary Assessment Study.
# For the original code and data, see https://osf.io/kj5zh

# --------------- 1. Loading packages & Data ------------------------------
# List of required packages
required_packages <- c(
  "Matrix", "fastmatrix", "BVAR", "brms", "expm", 
  "qgraph", "tidyverse", "ggplot2", "blavaan", "lavaan"
)

# Function to check and install missing packages
install_if_missing <- function(packages) {
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg)
    }
    library(pkg, character.only = TRUE)
  }
}

install_if_missing(required_packages)

# --------------- 2. Estimation -----------------------------


# load 
load(file.path(getwd(), "/Fried_2022 data/clean_network.RData")); gc()
Data5b <- Data2

# Variables to investigate:
vars <- c("Q1","Q2","Q3","Q4","Q5","Q6","Q7","Q8","Q9",
          "Q10","Q11","Q12","Q13","Q14","Q15","Q16","Q17","Q18")

# Labels:
varLabs <- c("Relax","Irritable","Worry","Nervous","Future","Anhedonia",
             "Tired","Hungry","Alone","Angry","Social_offline","Social_online","Music",
             "Procrastinate","Outdoors","C19_occupied","C19_worry","Home")
names(Data5b)[names(Data5b) %in% vars] <- varLabs

varLabs <- c("Relax", "Worry",  
             "Nervous")
varLabs2 <- c("Relax", "Worry",  
              "Nervous", "Tired", "Hungry",
              "Alone", "Angry")


# Lagged variables
Data5b <- Data5b %>%
  group_by(id) %>%
  mutate(across(all_of(varLabs2), 
                ~ lag(.x), 
                .names = "{.col}_lag")) %>%
  ungroup()


# Bayesian Lavaan estimation procedures - so that we can obtain WAIC, DIC, for model comparisons maybe. ------
## Network VAR(1) model: -----
model_lavaan <- "
  level:1
    Relax   ~ a1*Relax_lag + a2*Worry_lag + a3*Nervous_lag
    Worry   ~ b1*Relax_lag + b2*Worry_lag + b3*Nervous_lag
    Nervous ~ c1*Relax_lag + c2*Worry_lag + c3*Nervous_lag

    # Within-level residual covariances if desired:
    Relax   ~~  w1*Worry
    Relax   ~~  w2*Nervous
    Worry   ~~  w3*Nervous

  level: 2
    # Freed means for the random intercepts:
    Relax   ~ 1
    Worry   ~ 1
    Nervous ~ 1

    # Covariances among random intercepts:
    Relax   ~~ Worry
    Relax   ~~ Nervous
    Worry   ~~ Nervous
"

# Then fit using bsem() or bcfa() or blavaan(), e.g.:

fit_lavaan <- sem(
  model          = model_lavaan,
  data           = Data5b,
  cluster        = "id",      # multi-level SEM
  estimator      = "MLR" 
  )
summary(fit_lavaan)

## D-CF(1) model: --------
# NOTE: we could allow there to be autoregression of symptoms, i.e., due to measurement error.
model_lavaan_dcf <- "
  level:1
  
  CF_2 =~ l1*Relax + l2*Worry + l3*Nervous
  CF_1 =~ l1*Relax_lag + l2*Worry_lag + l3*Nervous_lag
  CF_2 ~ psi*CF_1
  
  CF_1 ~~ phi*CF_1
  CF_2 ~~ phi*CF_2

  level: 2
  # Freed means for the random intercepts:
  Relax   ~ 1
  Worry   ~ 1
  Nervous ~ 1

  # Covariances among random intercepts:
  Relax   ~~ Worry
  Relax   ~~ Nervous
  Worry   ~~ Nervous
  "

# Then fit using bsem() or bcfa() or blavaan(), e.g.:

fit_lavaan_dcf <- sem(
  model          = model_lavaan_dcf,
  data           = Data5b,
  cluster        = "id",      # multi-level SEM
  estimator      = "MLR"
  )
summary(fit_lavaan_dcf)
fitmeasures(fit_lavaan_dcf)

## DCF with serially correlated 'measurement errors' ------------

model_lavaan_dcf_error <- "
  level:1
  
  CF_2 =~ l1*Relax + l2*Worry + l3*Nervous
  CF_1 =~ l1*Relax_lag + l2*Worry_lag + l3*Nervous_lag
  
  CF_2 ~ psi*CF_1
  
  CF_1 ~~ phi*CF_1
  CF_2 ~~ phi*CF_2

  # Added serial correlation for 'measurement error'.
  Relax ~~ Relax_lag
  Worry ~~ Worry_lag
  Nervous ~~ Nervous_lag

  level: 2
  # Freed means for the random intercepts:
  Relax   ~ 1
  Worry   ~ 1
  Nervous ~ 1

  # Covariances among random intercepts:
  Relax   ~~ Worry
  Relax   ~~ Nervous
  Worry   ~~ Nervous
  "
# Then fit using bsem() or bcfa() or blavaan(), e.g.:

fit_lavaan_dcf_error <- sem(
  model          = model_lavaan_dcf_error,
  data           = Data5b,
  cluster        = "id",      # multi-level SEM
  estimator      = "MLR"
)
summary(fit_lavaan_dcf_error)
fitmeasures(fit_lavaan_dcf_error)

## compare models. -----

# The var has perfect fit, because the innovation covariance is 'allowed'. This perfectly explains the covariance of symptoms and lagged symptoms...
# Fair?

# Obtain fit indices of the dcf model
lavaan::anova(fit_lavaan_dcf, fit_lavaan, fit_lavaan_dcf_error)
lavaan::fitmeasures(fit_lavaan_dcf_error, fit.measures = c("rmsea.scaled", "cfi.scaled")) # Good fit?

# Lavaan, SEM, for more variables. --------------
varLabs2 <- c("Relax", "Worry", "Nervous", "Tired", "Hungry", "Alone", "Angry")

model_lavaan_more <- "
  level: 1
    # Regress each variable on all lagged predictors
    Relax   ~ a1*Relax_lag   + a2*Worry_lag   + a3*Nervous_lag + a4*Tired_lag + a5*Hungry_lag + a6*Alone_lag + a7*Angry_lag
    Worry   ~ b1*Relax_lag   + b2*Worry_lag   + b3*Nervous_lag + b4*Tired_lag + b5*Hungry_lag + b6*Alone_lag + b7*Angry_lag
    Nervous ~ c1*Relax_lag   + c2*Worry_lag   + c3*Nervous_lag + c4*Tired_lag + c5*Hungry_lag + c6*Alone_lag + c7*Angry_lag
    Tired   ~ d1*Relax_lag   + d2*Worry_lag   + d3*Nervous_lag + d4*Tired_lag + d5*Hungry_lag + d6*Alone_lag + d7*Angry_lag
    Hungry  ~ e1*Relax_lag   + e2*Worry_lag   + e3*Nervous_lag + e4*Tired_lag + e5*Hungry_lag + e6*Alone_lag + e7*Angry_lag
    Alone   ~ f1*Relax_lag   + f2*Worry_lag   + f3*Nervous_lag + f4*Tired_lag + f5*Hungry_lag + f6*Alone_lag + f7*Angry_lag
    Angry   ~ g1*Relax_lag   + g2*Worry_lag   + g3*Nervous_lag + g4*Tired_lag + g5*Hungry_lag + g6*Alone_lag + g7*Angry_lag

    # Optional: specify within-level residual covariances among variables
    # For brevity, only a subset is shown; add others as needed.
    Relax   ~~ w1*Worry + w2*Nervous + w3*Tired + w4*Hungry + w5*Alone + w6*Angry
    Worry   ~~ w7*Nervous + w8*Tired + w9*Hungry + w10*Alone + w11*Angry
    Nervous ~~ w12*Tired + w13*Hungry + w14*Alone + w15*Angry
    Tired   ~~ w16*Hungry + w17*Alone + w18*Angry
    Hungry  ~~ w19*Alone + w20*Angry
    Alone   ~~ w21*Angry

  level: 2
    # Freed means for the random intercepts
    Relax   ~ 1
    Worry   ~ 1
    Nervous ~ 1
    Tired   ~ 1
    Hungry  ~ 1
    Alone   ~ 1
    Angry   ~ 1

    # Covariances among random intercepts for all variables
    Relax   ~~ Relax + Worry + Nervous + Tired + Hungry + Alone + Angry
    Worry   ~~ Worry + Nervous + Tired + Hungry + Alone + Angry
    Nervous ~~ Nervous + Tired + Hungry + Alone + Angry
    Tired   ~~ Tired + Hungry + Alone + Angry
    Hungry  ~~ Hungry + Alone + Angry
    Alone   ~~ Alone + Angry
    Angry   ~~ Angry
"

# Fit the model (adjust control parameters as necessary)
fit_lavaan_more <- sem(
  model   = model_lavaan_more,
  data    = Data5b,
  cluster = "id",
  estimator = "MLR"
)

summary(fit_lavaan_more)

## DCF -------------


model_lavaan_more_dcf <- "
  level: 1
  
  # Define latent factors for current and lagged variables
  CF_2 =~ l1*Relax + l2*Worry + l3*Nervous + l4*Tired + l5*Hungry + l6*Alone + l7*Angry
  CF_1 =~ l1*Relax_lag + l2*Worry_lag + l3*Nervous_lag + l4*Tired_lag + l5*Hungry_lag + l6*Alone_lag + l7*Angry_lag
  
  # Regression of current factor on lagged factor
  CF_2 ~ psi*CF_1
  
  # Variances for the latent factors
  CF_1 ~~ phi*CF_1
  CF_2 ~~ phi*CF_2


  level: 2
    # Freed means for the random intercepts for all variables
    Relax   ~ 1
    Worry   ~ 1
    Nervous ~ 1
    Tired   ~ 1
    Hungry  ~ 1
    Alone   ~ 1
    Angry   ~ 1

    # Covariances among random intercepts for all pairs of variables
    Relax   ~~ Worry + Nervous + Tired + Hungry + Alone + Angry
    Worry   ~~ Nervous + Tired + Hungry + Alone + Angry
    Nervous ~~ Tired + Hungry + Alone + Angry
    Tired   ~~ Hungry + Alone + Angry
    Hungry  ~~ Alone + Angry
    Alone   ~~ Angry
"

# Fit the model using lavaan's SEM function:
fit_lavaan_more_dcf <- sem(
  model   = model_lavaan_more_dcf,
  data    = Data5b,
  cluster = "id",
  estimator = "MLR"
)

summary(fit_lavaan_more_dcf, standardized = TRUE)
fitmeasures(fit_lavaan_more_dcf)


## DCF with error -------

model_lavaan_more_dcf_error <- "
  level: 1
  
  # Define latent factors for current and lagged variables
  CF_2 =~ l1*Relax + l2*Worry + l3*Nervous + l4*Tired + l5*Hungry + l6*Alone + l7*Angry
  CF_1 =~ l1*Relax_lag + l2*Worry_lag + l3*Nervous_lag + l4*Tired_lag + l5*Hungry_lag + l6*Alone_lag + l7*Angry_lag
  
  # Regression of current factor on lagged factor
  CF_2 ~ psi*CF_1
  
  # Variances for the latent factors
  CF_1 ~~ phi*CF_1
  CF_2 ~~ phi*CF_2

  # Serial correlations for measurement error across time for each variable
  Relax   ~~ Relax_lag
  Worry   ~~ Worry_lag
  Nervous ~~ Nervous_lag
  Tired   ~~ Tired_lag
  Hungry  ~~ Hungry_lag
  Alone   ~~ Alone_lag
  Angry   ~~ Angry_lag

  level: 2
    # Freed means for the random intercepts for all variables
    Relax   ~ 1
    Worry   ~ 1
    Nervous ~ 1
    Tired   ~ 1
    Hungry  ~ 1
    Alone   ~ 1
    Angry   ~ 1

    # Covariances among random intercepts for all pairs of variables
    Relax   ~~ Worry + Nervous + Tired + Hungry + Alone + Angry
    Worry   ~~ Nervous + Tired + Hungry + Alone + Angry
    Nervous ~~ Tired + Hungry + Alone + Angry
    Tired   ~~ Hungry + Alone + Angry
    Hungry  ~~ Alone + Angry
    Alone   ~~ Angry
"

# Fit the model using lavaan's SEM function:
fit_lavaan_more_dcf_error <- sem(
  model   = model_lavaan_more_dcf_error,
  data    = Data5b,
  cluster = "id",
  estimator = "MLR"
)

summary(fit_lavaan_more_dcf_error, standardized = TRUE)
fitmeasures(fit_lavaan_more_dcf_error)


