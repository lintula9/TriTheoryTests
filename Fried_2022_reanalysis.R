# Original code was created by Eiko Fried, Faidra Papanikolaou & Sascha Epskamp. 
# This version of code was edited in 29.11.2024 for reanalysis of Fried(2022)
# Mental Health and Social Contact During the COVID-19 Pandemic: 
# An Ecological Momentary Assessment Study.
# For the original code and data, see https://osf.io/kj5zh

# --------------- 1. Loading packages & Data ------------------------------
# List of required packages
required_packages <- c(
  "lme4", "ggplot2", "dplyr", "tidyr", 
  "qgraph", "tidyselect", "mlVAR", 
  "viridis", "summarytools", "lm.beta",
  "Matrix", "fastmatrix"
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

# Install and load packages
install_if_missing(required_packages)


# --------------- 2. Detrending fixed effects -----------------------------

# We do not do multilevel detrending per person

# load 
load(file.path(getwd(), "Fried_2022 data/clean_network.RData")); gc()
Data5b <- Data2

# Alpha to detrend:
alpha <- 0.05

# Variables to investigate:
vars <- c("Q1","Q2","Q3","Q4","Q5","Q6","Q7","Q8","Q9",
          "Q10","Q11","Q12","Q13","Q14","Q15","Q16","Q17","Q18")

# Labels:
varLabs <- c("Relax","Irritable","Worry","Nervous","Future","Anhedonia",
             "Tired","Hungry","Alone","Angry","Social_offline","Social_online","Music",
             "Procrastinate","Outdoors","C19_occupied","C19_worry","Home")
names(Data5b)[names(Data5b) %in% vars] <- varLabs



# Select a plausible set of variables, which could be explained by a one CF model:
Data5b <- as_tibble(Data5b %>% select("Relax", "Worry", 
                            "Nervous", "Future",
                            "id", "beep", "day", "conc"))
varLabs <- c("Relax", "Worry",  #Adjustment for below to work. 03.12.2024. Sakari Lintula
            "Nervous", "Future")

# Data frame with empty values for fitted effects (all):
fitted_all <- expand.grid(
  beep = seq(min(Data5b$beep),max(Data5b$beep)),
  day = seq(min(Data5b$day),max(Data5b$day))
)

# Data frame with empty values for day trends:
fitted_day <- data.frame(
  day = seq(min(Data5b$day),max(Data5b$day))
)

# Data frame with empty values for beeps:
fitted_beep <- data.frame(
  beep = seq(min(Data5b$beep),max(Data5b$beep))
)

# Data frame to store p-values:
p_values <- data.frame(
  var = c("day", "beep")
)

# Also empty data frame list for test statistics:
testStatistics <- list()
coefficients <- list()
stdcoefficients <- list()

# Make the beep variable factor in dataset:
Data5b$beepFactor <- factor(Data5b$beep, levels = 0:3, labels = c("09:00 - 12:00","12:00 - 15:00","15:00 - 18:00","18:00 - 21:00"))
fitted_all$beepFactor <- factor(fitted_all$beep, levels = 0:3, labels = c("09:00 - 12:00","12:00 - 15:00","15:00 - 18:00","18:00 - 21:00"))
fitted_beep$beepFactor <- factor(fitted_beep$beep, levels = 0:3, labels = c("09:00 - 12:00","12:00 - 15:00","15:00 - 18:00","18:00 - 21:00"))

# Make day variable for dates:
Data5b$date <- as.Date("2020-03-15") + Data5b$day
fitted_all$date <- as.Date("2020-03-15") + fitted_all$day
fitted_day$date <- as.Date("2020-03-15") + fitted_day$day

# Add the midpoints as time variable:
Data5b$midTime <- as.character(factor(Data5b$beep, levels = 0:3, labels = c("10:30","13:30","16:30","19:30")))
Data5b$midTime <- as.POSIXct(paste(Data5b$date,Data5b$midTime), format = "%Y-%m-%d %H:%M", tz = "Europe/Amsterdam")

fitted_all$midTime <- as.character(factor(fitted_all$beep, levels = 0:3, labels = c("10:30","13:30","16:30","19:30")))
fitted_all$midTime <- as.POSIXct(paste(fitted_all$date,fitted_all$midTime), format = "%Y-%m-%d %H:%M", tz = "Europe/Amsterdam")

# Data frame to store detrended data:
data_detrended <- Data5b

# Fix curves:
for (v in seq_along(varLabs)){
  formula <- as.formula(paste0(varLabs[v], " ~ 1 + day + factor(beep)"))
  lmRes <- lm(formula, data = Data5b)
  
  # Fixed effects:
  fixed <- coef(lmRes)
  
  # make zero if not significant at alpha:
  p_values[[varLabs[v]]] <- anova(lmRes)[["Pr(>F)"]][1:2]
  if (p_values[,varLabs[v]][1] > alpha){
    fixed[2] <- 0
  }
  if (p_values[,varLabs[v]][2] > alpha){
    fixed[3:5] <- 0
  }
  
  # Add to DFs:
  fitted_all[,varLabs[v]] <- fixed[1] + fixed[2] * fitted_all[["day"]]  +  fixed[3] * (fitted_all[["beep"]] == 1)  + 
    fixed[4] * (fitted_all[["beep"]] == 2) + fixed[5] *  (fitted_all[["beep"]] == 3)
  
  fitted_day[,varLabs[v]] <- fixed[1] + fixed[2] * fitted_day[["day"]]
  
  fitted_beep[,varLabs[v]] <- fixed[1] + fixed[2] * median(fitted_day[["day"]]) +  fixed[3] * (fitted_beep[["beep"]] == 1)  + 
    fixed[4] * (fitted_beep[["beep"]] == 2) + fixed[5] *  (fitted_beep[["beep"]] == 3)
  
  # Detrend data:
  data_detrended[,varLabs[v]] <- Data5b[,varLabs[v]] - (fixed[1] + fixed[2] * Data5b[["day"]]  +  fixed[3] * (Data5b[["beep"]] == 1)  + 
                                                          fixed[4] * (Data5b[["beep"]] == 2) + fixed[5] *  (Data5b[["beep"]] == 3))
  
  # Test statistic dataframe:   
  ### TODOSACHA: add cohen's D, for prepost I use 
  ### df6 %>% rstatix::cohens_d(value ~ type, paired = TRUE) 
  ids <- rownames(anova(lmRes))
  testStatistics[[v]] <- cbind(data.frame(var = varLabs[v], effect = ids), anova(lmRes))
  
  coefficients[[v]] <- data.frame(
    var = varLabs[v],
    type = names(coef(lmRes)),
    coef = coef(lmRes),
    std = coef(lm.beta(lmRes))
  )
}

# Make sense of significant trends for descriptives

testStatistics <- do.call(rbind, testStatistics)
coefficients <- do.call(rbind, coefficients)

testStatistics[,7] <- sprintf("%.5f", testStatistics[,7])

testStatistics_sig_day <- subset(testStatistics, effect == "day"); testStatistics_sig_day # this has the p-values for the standardized slopes one line below
coefficients_day <- subset(coefficients, coefficients[,2] == "day"); coefficients_day     # this shows standardized slopes over time; e.g. C19_occupied decreased by 0.18; home increased slightly by 0.03

testStatistics_sig_beep <- subset(testStatistics, testStatistics[,7] < 0.05 & effect == "factor(beep)"); testStatistics_sig_beep
coefficients_beep<- subset(coefficients, coefficients[,2] == "factor(beep)1" |  
                             coefficients[,2] == "factor(beep)2" | 
                             coefficients[,2] == "factor(beep)3"); coefficients_beep


# Plot versus averaged:
Data5b_avg <- Data5b %>% group_by(day,beep,conc,beepFactor,date,midTime) %>% summarize_at(varLabs, list(~mean(., na.rm=TRUE)))

# to long:
fitted_long <- fitted_all %>% pivot_longer(cols=all_of(varLabs), values_to = "fitted")
mean_long <- Data5b_avg %>% pivot_longer(cols=all_of(varLabs), values_to = "mean")
longDF <- left_join(fitted_long,mean_long) %>% pivot_longer(cols = c("fitted" , "mean"), names_to = "type")

# Dates:
Sys.setlocale("LC_ALL","C")
breaks <- as.POSIXct(paste0(as.Date("2020-03-15") + 1:14, " 00:01"), format = "%Y-%m-%d %H:%M")


### Plot the trends that are removed:
p1 <- ggplot(longDF, aes(x = midTime, y = value, lty = type, pch = type)) + 
  geom_line() + 
  facet_grid(name ~ ., scales = "free") + 
  scale_x_time(breaks = breaks, labels = format(breaks, "%B %d")) + 
  xlab("") + ylab("") + theme_bw() + 
  scale_linetype_discrete("") + 
  scale_colour_discrete("") + 
  ggtitle("Average time-series (dashed) and fitted fixed-effects (solid)",
          subtitle = "Effects not significant at alpha = 0.05 set to zero in fitted model") + 
  theme(legend.position = "none")

p2 <- ggplot(longDF, aes(x = midTime, y = value, lty = type, pch = type)) + 
  geom_line() + 
  facet_grid(name ~ .) + ylim(range(longDF$value)) + 
  scale_x_time(breaks = breaks, labels = format(breaks, "%B %d")) + 
  xlab("") + ylab("") + theme_bw() + 
  scale_linetype_discrete("") + 
  scale_colour_discrete("") + 
  ggtitle("Average time-series (dashed) and fitted fixed-effects (solid)",subtitle = "Effects not significant at alpha = 0.05 set to zero in fitted model") + 
  theme(legend.position = "none")


# Plot detrended data:
data_detrended_avg <- data_detrended %>% group_by(day,beep,conc,beepFactor,date,midTime) %>% summarize_at(varLabs, list(~mean(., na.rm=TRUE)))
detrended_long <- data_detrended_avg %>% pivot_longer(cols=all_of(varLabs), values_to = "mean")

### Plot the trends that are removed:
p3 <- ggplot(detrended_long, aes(x = midTime, y = mean)) + 
  geom_line() + 
  facet_grid(name ~ ., scales = "free") + 
  scale_x_time(breaks = breaks, labels = format(breaks, "%B %d")) + 
  xlab("") + ylab("") + theme_bw() + 
  scale_linetype_discrete("") + 
  scale_colour_discrete("")  + 
  ggtitle("Detrended data according to fixed-effects model",subtitle = "Effects not significant at alpha = 0.05 set to zero in fitted model")

p4 <- ggplot(detrended_long, aes(x = midTime, y = mean)) + 
  geom_line() + 
  facet_grid(name ~ .) + ylim(range(detrended_long$mean)) + 
  scale_x_time(breaks = breaks, labels = format(breaks, "%B %d")) + 
  xlab("") + ylab("") + theme_bw() + 
  scale_linetype_discrete("") + 
  scale_colour_discrete("")  + 
  ggtitle("Detrended data according to fixed-effects model (fixed y-axis)",subtitle = "Effects not significant at alpha = 0.05 set to zero in fitted model")


# Save detrended data for students to work with
# saveRDS(data_detrended, file=paste0(datapath, "data_students_network.RData"))



# -------------------------------------------------------------------------
# --------------- 3. Fit mlVAR network models: orthogonal estimation ------
# -------------------------------------------------------------------------

# load(paste0(datapath, "network_orthogonal.RData"))

# Run mlVAR (orthogonal):
res <- mlVAR(data_detrended,
             vars=varLabs,
             idvar="id",
             dayvar="day",
             beepvar="beep",
             lags = 1,
             temporal = "orthogonal",
             contemporaneous = "orthogonal",
             nCores = 8)
# save(res, file=paste0(datapath, "network_orthogonal.RData"))

# Plot :
names <- c("Relax","Irritable","Worry","Nervous","Future",
           "Anhedonia","Tired","Alone",
           "Social-offline", "Social-online", "Outdoors",
           "C19-occupied", "C19-worry", "Home")

gr <- list('Stress'=c(1:7), 'Social'=c(8:10), 'COVID-19'=c(11:14))



# Obtain the closest indistinguishable model ------

# Get the results: Between coefficient matrix as well as between innovation covariance.
A = res$results$Beta$mean[,,1]
Z = res$results$Theta$cov$mean

# Solve for the covariance
library(expm); library(Matrix)
Sigma_VAR = matrix(solve(diag(1, ncol = ncol(A)^2, nrow = nrow(A)^2) - fastmatrix::kronecker.prod(A)) %*% fastmatrix::vec(Z),
                   ncol = ncol(A), nrow = nrow(A))
K_VAR = function( Delta ) {
  
  A %^% Delta %*% Sigma_VAR
  }

eigen(Sigma_VAR) # Suggestive of there being one large component. 
# Capture Lambda as the normalized eigenvector of the largest eigenvalue.

L = eigen(Sigma_VAR)$vectors[,1] / c(sqrt(t(eigen(Sigma_VAR)$vectors[,1]) %*% eigen(Sigma_VAR)$vectors[,1]))
eigen(A)$vectors[,1] # The first vector of A is nearly in the same direction.
eigen(A)$vectors[,1] - (-1)*L 

  # Project A onto Lambda
A_0 =  L %*% t(L) %*% A %*% L %*% t(L)

  # Project the difference onto the orthogonal complement of Lambda
B_tilde = (A - A_0) %*% (diag(1, nrow = nrow(A), ncol = ncol(A)) - L %*% t(L))

# Create the VAR(1), indistinguishable from one dimensional D-CF(1) model.
A_tilde = A_0 + B_tilde
Z_tilde = L %*% t(L) %*% Z %*% L %*% t(L)# Project onto Lambda, as it must be proportional to LL^T


par(mfrow=c(2,2))
qgraph(Z, layout = "circle", labels = expression(X[1], X[2], X[3], X[4]))
title("VAR(1) and indistinguishable VAR(1).",outer = T)

qgraph(Z_tilde, layout = "circle", labels = expression(X[1], X[2], X[3], X[4]))
qgraph(A_tilde, layout = "circle", labels = expression(X[1], X[2], X[3], X[4]))
qgraph(A, layout = "circle", labels = expression(X[1], X[2], X[3], X[4]))
par(mfrow=c(1,1))

Sigma_DCF = matrix(solve(diag(1, ncol = ncol(A)^2, nrow = nrow(A)^2) - fastmatrix::kronecker.prod(A_tilde)) %*% fastmatrix::vec(Z_tilde),
                   ncol = ncol(A), nrow = nrow(A))
K_VAR = function( Delta ) {
  
  A_tilde %^% Delta %*% Sigma_DCF
}
