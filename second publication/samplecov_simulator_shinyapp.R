# ==============================================================================
# Robust OU Process Simulator (Final Crash-Proof Version)
# ==============================================================================

library(shiny)
library(ggplot2)
library(dplyr)
library(tidyr)

# --- 1. NUMERICAL HELPERS ---

get_matrix_sqrt <- function(Q) {
  # Clean potential NAs or Infs from the input matrix
  Q[!is.finite(Q)] <- 0
  tryCatch({
    t(chol(Q)) 
  }, error = function(e) {
    e_Q <- eigen(Q, symmetric = TRUE)
    val <- pmax(Re(e_Q$values), 0) 
    return(Re(e_Q$vectors %*% diag(sqrt(val)) %*% t(e_Q$vectors)))
  })
}

mat_exp <- function(A) {
  e <- eigen(A)
  V <- e$vectors
  D <- diag(exp(e$values))
  res <- Re(V %*% D %*% solve(V))
  res[!is.finite(res)] <- 0
  return(res)
}

sigmoid <- function(t, t0, k) {
  1 / (1 + exp(-k * (t - t0)))
}

# --- 2. GENERATION LOGIC ---

generate_system_matrices <- function(K, Rank) {
  if(Rank > K) Rank <- K
  
  # 1. Diffusion Sigma
  rnd <- matrix(rnorm(K^2), K, K)
  Rot <- qr.Q(qr(rnd))
  eigs <- c(runif(Rank, 0.5, 1.5), rep(0, K - Rank))
  eigs <- sort(eigs, decreasing = TRUE)
  Sigma <- Rot %*% diag(eigs)
  
  # 2. Drift B (Ensure it is strictly stable)
  B_raw <- matrix(rnorm(K^2, mean = 0, sd = 0.5), K, K)
  # Symmetrize to make eigenvalues real and easier to stabilize
  B_sym <- (B_raw + t(B_raw)) / 2
  ev <- eigen(B_sym)$values
  # Shift eigenvalues so the minimum is at least 0.5
  B <- B_sym + diag(K) * (max(0.5, -min(ev) + 0.5))
  
  return(list(B = B, Sigma = Sigma))
}

# --- 3. SIMULATION CORE ---

generate_trajectory <- function(K, B, mu, Sigma, meas_noise_sd, scenario, steepness, steps=2000) {
  dt <- 0.1
  # Discrete-time transition: Phi = exp(-B * dt)
  Phi <- mat_exp(-B * dt)
  
  X <- matrix(0, nrow = steps, ncol = K)
  X[1, ] <- mu
  curr <- numeric(K)
  
  Q_base <- (Sigma %*% t(Sigma)) * dt
  L_base <- get_matrix_sqrt(Q_base)
  t_mid <- steps / 2
  
  for (t in 1:(steps - 1)) {
    L_t <- L_base
    if (scenario != "Stationary") {
      scale_factor <- if (scenario == "Decay") 1 - sigmoid(t, t_mid, steepness) else sigmoid(t, t_mid, steepness)
      L_t <- L_base * scale_factor
    }
    
    Z <- rnorm(K)
    curr <- Phi %*% curr + L_t %*% Z
    
    # Anti-explosion check
    curr[!is.finite(curr)] <- 0
    if(max(abs(curr)) > 1e6) curr <- (curr / max(abs(curr))) * 1e6
    
    X[t+1, ] <- curr + mu
  }
  
  Y <- X
  if (meas_noise_sd > 0) {
    Noise <- matrix(rnorm(steps * K, sd = meas_noise_sd), nrow = steps, ncol = K)
    Y <- Y + Noise
  }
  return(data.frame(Time = 1:steps, Y))
}

run_eigen_sim <- function(K, B, mu, Sigma, meas_noise_sd, scenario, steepness, reps, N_seq) {
  total_steps <- 2000
  results_list <- vector("list", reps * length(N_seq))
  idx <- 1
  for(r in 1:reps) {
    traj_df <- generate_trajectory(K, B, mu, Sigma, meas_noise_sd, scenario, steepness, total_steps)
    Y_mat <- as.matrix(traj_df[, -1])
    for(n in N_seq) {
      step_size <- max(floor(total_steps / n), 1)
      start_idx <- sample(1:step_size, 1)
      indices <- seq(start_idx, by = step_size, length.out = n)
      indices <- indices[indices <= total_steps]
      Y_sub <- Y_mat[indices, , drop=FALSE]
      
      # Final safety check before eigen
      cv <- cov(Y_sub)
      if(any(!is.finite(cv)) || nrow(Y_sub) < 2) {
        eigs <- rep(0, K)
      } else {
        eigs <- sort(Re(eigen(cv, only.values = TRUE)$values), decreasing = TRUE)
      }
      
      results_list[[idx]] <- data.frame(Rep = r, N = n, Rank = paste0("Eig", 1:K), Value = eigs)
      idx <- idx + 1
    }
  }
  return(do.call(rbind, results_list))
}

# --- 4. UI ---

ui <- fluidPage(
  theme = bslib::bs_theme(version = 4, bootswatch = "flatly"),
  titlePanel("OU Process: Subsampling Simulator"),
  sidebarLayout(
    sidebarPanel(
      width = 4,
      h4("1. System Design"),
      splitLayout(
        numericInput("K_in", "Dim (K)", value = 3, min = 2, max = 10),
        numericInput("R_in", "Rank (R)", value = 2, min = 1, max = 10)
      ),
      actionButton("gen_sys_btn", "Generate New System", class="btn-primary btn-block"),
      br(),
      textAreaInput("B_txt", "Drift Matrix (B)", rows = 4),
      textAreaInput("S_txt", "Diffusion Matrix (Sigma)", rows = 4),
      hr(),
      h4("2. Experiment"),
      selectInput("scenario", "Scenario", choices = c("Stationary", "Decay", "Growth")),
      sliderInput("steepness", "Transition Speed", min=0.01, max=0.5, value=0.1),
      numericInput("noise_sd", "Measurement Error (SD)", value = 0.2, step = 0.1),
      actionButton("run_btn", "Run Simulation", class = "btn-success btn-lg btn-block")
    ),
    mainPanel(
      width = 8,
      tabsetPanel(
        tabPanel("Eigenvalues",
                 br(),
                 plotOutput("eigPlot", height = "500px"),
                 wellPanel(verbatimTextOutput("theo_eigs"))
        ),
        tabPanel("Trajectory",
                 br(),
                 plotOutput("tsPlot", height = "500px"))
      )
    )
  )
)

# --- 5. SERVER ---

server <- function(input, output, session) {
  
  mat_to_str <- function(M) paste(apply(M, 1, function(x) paste(round(x, 3), collapse=", ")), collapse=",\n")
  
  # Triggered by Dimension change or Button
  observeEvent(c(input$gen_sys_btn, input$K_in), {
    req(input$K_in, input$R_in)
    sys <- generate_system_matrices(input$K_in, min(input$R_in, input$K_in))
    updateTextAreaInput(session, "B_txt", value = mat_to_str(sys$B))
    updateTextAreaInput(session, "S_txt", value = mat_to_str(sys$Sigma))
  }, ignoreInit = FALSE)
  
  get_matrices <- reactive({
    req(input$B_txt, input$S_txt, input$K_in)
    K <- input$K_in
    
    clean_txt <- function(s) {
      s <- gsub("[^0-9.-]+", ",", s) # Keep only numbers, dots, and minus
      vals <- as.numeric(unlist(strsplit(s, ",")))
      vals <- vals[!is.na(vals)]
      return(vals)
    }
    
    B_vals <- clean_txt(input$B_txt)
    S_vals <- clean_txt(input$S_txt)
    
    if(length(B_vals) != K*K || length(S_vals) != K*K) return(NULL)
    
    return(list(K=K, B=matrix(B_vals, K, K, byrow=TRUE), S=matrix(S_vals, K, K, byrow=TRUE), mu=rep(0, K)))
  })
  
  output$theo_eigs <- renderPrint({
    mats <- get_matrices()
    if(is.null(mats)) return("Parsing inputs... ensure matrix has K*K values.")
    
    BBt <- mats$S %*% t(mats$S)
    I_mat <- diag(mats$K)
    L_op <- kronecker(I_mat, mats$B) + kronecker(mats$B, I_mat)
    vec_C <- try(solve(L_op, as.vector(BBt)), silent=TRUE)
    
    if(inherits(vec_C, "try-error")) return("Lyapunov Solver Failed (System Unstable).")
    eigs <- sort(eigen(matrix(vec_C, mats$K))$values, decreasing=TRUE)
    cat("Stationary Eigenvalues: ", paste(round(Re(eigs), 4), collapse=", "))
  })
  
  sim_data <- reactiveValues(eigs = NULL, ts = NULL)
  
  observeEvent(input$run_btn, {
    mats <- get_matrices()
    req(mats)
    withProgress(message = "Processing...", {
      res <- run_eigen_sim(mats$K, mats$B, mats$mu, mats$S, input$noise_sd, input$scenario, input$steepness, 10, seq(50, 1000, 50))
      sim_data$eigs <- res
      
      ts_raw <- generate_trajectory(mats$K, mats$B, mats$mu, mats$S, input$noise_sd, input$scenario, input$steepness, 1000)
      sim_data$ts <- ts_raw[, 1:min(mats$K+1, 5)] %>% 
        tidyr::pivot_longer(cols = -Time, names_to = "Variable", values_to = "Value")
    })
  })
  
  output$eigPlot <- renderPlot({
    req(sim_data$eigs)
    mats <- get_matrices()
    
    summ <- sim_data$eigs %>% 
      group_by(N, Rank) %>% 
      summarise(Med=median(Value, na.rm=TRUE), Lo=quantile(Value, 0.1, na.rm=TRUE), Hi=quantile(Value, 0.9, na.rm=TRUE), .groups="drop")
    
    p <- ggplot(summ, aes(x=N, group=Rank)) +
      geom_ribbon(aes(ymin=Lo, ymax=Hi, fill=Rank), alpha=0.2) +
      geom_line(aes(y=Med, color=Rank), size = 1) +
      scale_color_discrete(type = hcl.colors(mats$K, "Viridis")) +
      scale_fill_discrete(type = hcl.colors(mats$K, "Viridis")) +
      labs(title = "Eigenvalue Convergence over N", subtitle = "Dashed lines = Stationary Theory", y="Value") +
      theme_minimal()
    
    # Overlay Theory
    BBt <- mats$S %*% t(mats$S)
    L_op <- kronecker(diag(mats$K), mats$B) + kronecker(mats$B, diag(mats$K))
    vec_C <- try(solve(L_op, as.vector(BBt)), silent=TRUE)
    if(!inherits(vec_C, "try-error")) {
      true_vals <- sort(Re(eigen(matrix(vec_C, mats$K))$values), decreasing=TRUE)
      true_df <- data.frame(Rank=paste0("Eig", 1:mats$K), Value=true_vals)
      p <- p + geom_hline(data=true_df, aes(yintercept=Value, color=Rank), linetype="dashed", alpha=0.6)
    }
    p
  })
  
  output$tsPlot <- renderPlot({
    req(sim_data$ts)
    ggplot(sim_data$ts, aes(x=Time, y=Value, color=Variable)) + 
      geom_line() + 
      labs(title = "Sample Trajectory (Max 4 Dims)") +
      theme_minimal()
  })
}

shinyApp(ui, server)