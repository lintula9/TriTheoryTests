# Intro plotting & statistics.
intro_simulation <- readRDS("Datas/intro_simulation.rds")
source("Libraries.R")
# Summary of sample n past 2, 5, 10, 20, 30, 40, 50 and 
# est_dimension - true_dimension ≠ 0.
n_cuts           <- c(2,5,10,20,30,40,50)
n_past_names     <- paste0("n_past_",n_cuts)
intro_simulation[,n_past_names] <- sapply(
  n_cuts, FUN = \(x) intro_simulation$sample_n >= x)
intro_simulation <- intro_simulation %>%
  mutate(obs_est_diff = true_dimension - est_dimension)
# ggplot with coloring
ggplot(intro_simulation |> filter(sample_n < 20))


sapply(n_past_names, 
       FUN = \(nam) intro_simulation[intro_simulation[[nam]],] |>
         group_by()
         pull(obs_est_diff) |> (\(x) sum(x == 0) / length(x))() )