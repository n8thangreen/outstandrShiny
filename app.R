library(shiny)
library(bslib)
library(DT)
library(dplyr)
library(sandwich)

# Define UI
ui <- page_navbar(
  title = "Indirect Treatment Comparison Tool",
  
  # Data Input Tab
  nav_panel(
    title = "Data Input & Analysis",
    page_sidebar(
      sidebar = sidebar(
        width = 300,
        
        # File inputs
        card(
          card_header("Upload Data"),
          fileInput("ipd_file", "Upload Individual Patient Data (CSV)", accept = c(".csv")),
          fileInput("agg_file", "Upload Aggregate Data (CSV)", accept = c(".csv"))
        ),
        
        # Data type specification
        card(
          card_header("Specify Data Type"),
          radioButtons("data_type", "Data Type:",
                      choices = c("Binary", "Continuous"),
                      selected = "Binary")
        ),
        
        # Model selection
        card(
          card_header("Select Models to Fit"),
          checkboxGroupInput("models", "Models:",
                           choices = c("MAIC" = "maic", 
                                      "STC" = "stc"),
                           selected = "maic")
        ),
        
        # Outcome scale selection
        card(
          card_header("Outcome Scale"),
          radioButtons("outcome_scale", "Select Outcome Scale:",
                      choices = c("Log Odds Ratio" = "lor", 
                                 "Risk Ratio" = "rr"),
                      selected = "lor")
        ),
        
        # Variables selection UI will be generated dynamically based on data
        uiOutput("variable_selection"),
        
        # Run analysis button
        actionButton("run_analysis", "Run Analysis", class = "btn-primary")
      ),
      
      # Main panel
      card(
        card_header("Data Preview"),
        tabsetPanel(
          tabPanel("Individual Patient Data", DTOutput("ipd_preview")),
          tabPanel("Aggregate Data", DTOutput("agg_preview"))
        )
      )
    )
  ),
  
  # Results Tab
  nav_panel(
    title = "Results",
    card(
      card_header("Analysis Results"),
      DTOutput("results_table")
    )
  )
)

# Define server logic
server <- function(input, output, session) {
  
  # Reactive values to store the uploaded data
  ipd_data <- reactiveVal(NULL)
  agg_data <- reactiveVal(NULL)
  
  # Process the uploaded IPD file
  observeEvent(input$ipd_file, {
    req(input$ipd_file)
    tryCatch({
      data <- read.csv(input$ipd_file$datapath, stringsAsFactors = FALSE)
      ipd_data(data)
    }, error = function(e) {
      showNotification(paste("Error reading IPD file:", e$message), type = "error")
    })
  })
  
  # Process the uploaded aggregate data file
  observeEvent(input$agg_file, {
    req(input$agg_file)
    tryCatch({
      data <- read.csv(input$agg_file$datapath, stringsAsFactors = FALSE)
      agg_data(data)
    }, error = function(e) {
      showNotification(paste("Error reading aggregate file:", e$message), type = "error")
    })
  })
  
  # Dynamic UI for variable selection based on uploaded data
  output$variable_selection <- renderUI({
    req(ipd_data())
    ipd <- ipd_data()
    
    card(
      card_header("Variable Selection"),
      selectInput("outcome_var", "Outcome Variable:", choices = names(ipd)),
      selectInput("treatment_var", "Treatment Variable:", choices = names(ipd)),
      checkboxGroupInput("covariates", "Covariates for Adjustment:", 
                         choices = setdiff(names(ipd), c(input$outcome_var, input$treatment_var)))
    )
  })
  
  # Data previews
  output$ipd_preview <- renderDT({
    req(ipd_data())
    datatable(ipd_data(), options = list(pageLength = 5, scrollX = TRUE))
  })
  
  output$agg_preview <- renderDT({
    req(agg_data())
    datatable(agg_data(), options = list(pageLength = 5, scrollX = TRUE))
  })
  
  # Function to run MAIC
  run_maic <- function(ipd, agg, outcome_var, treatment_var, covariates, data_type, outcome_scale) {
    # Extract aggregate data means/proportions for matching
    agg_means <- as.numeric(agg[1, covariates])
    
    # Function to compute weights
    X <- as.matrix(ipd[, covariates])
    
    # Optimization function to find alpha
    loglik <- function(alpha) {
      sum(exp(X %*% alpha)) - sum(agg_means * alpha)
    }
    
    # Find optimal alpha using optimization
    res <- optim(rep(0, length(covariates)), loglik, method = "BFGS")
    alpha <- res$par
    
    # Calculate weights
    weights <- exp(X %*% alpha)
    
    # Normalize weights
    weights <- weights * (nrow(ipd) / sum(weights))
    
    # Apply weights based on data type and outcome scale
    if (data_type == "Binary") {
      # For binary outcomes, use weighted logistic regression
      formula <- as.formula(paste(outcome_var, "~", treatment_var))
      model <- glm(formula, data = ipd, family = binomial(), weights = weights)
      
      # Extract coefficient and SE
      coef <- coef(model)[treatment_var]
      vcov_robust <- sandwich::vcovHC(model, type = "HC0")
      se <- sqrt(diag(vcov_robust))[treatment_var]
      
      if (outcome_scale == "rr") {
        # Convert log odds ratio to log relative risk (approximation)
        p0 <- weighted.mean(ipd[[outcome_var]][ipd[[treatment_var]] == 0], 
                           weights[ipd[[treatment_var]] == 0])
        coef <- log(exp(coef) / (1 - p0 + p0 * exp(coef)))
        # Approximation for SE
        se <- se * abs(coef / coef(model)[treatment_var])
      }
    } else if (data_type == "Continuous") {
      # For continuous outcomes, use weighted linear regression
      formula <- as.formula(paste(outcome_var, "~", treatment_var))
      model <- lm(formula, data = ipd, weights = weights)
      
      coef <- coef(model)[treatment_var]
      vcov_robust <- sandwich::vcovHC(model, type = "HC0")
      se <- sqrt(diag(vcov_robust))[treatment_var]
    }
    
    # Calculate p-value and CI
    p_value <- 2 * (1 - pnorm(abs(coef / se)))
    lower_ci <- coef - 1.96 * se
    upper_ci <- coef + 1.96 * se
    
    return(list(
      Method = "MAIC",
      Estimate = coef,
      SE = se,
      Lower_CI = lower_ci,
      Upper_CI = upper_ci,
      P_value = p_value
    ))
  }
  
  # Run analysis when button is clicked
  analysis_results <- eventReactive(input$run_analysis, {
    req(ipd_data(), agg_data(), input$outcome_var, input$treatment_var, length(input$models) > 0)
    
    # Show notification
    showNotification("Running analysis...", type = "message", duration = NULL, id = "analysis")
    
    # Initialize results dataframe
    results <- data.frame(
      Method = character(),
      Estimate = numeric(),
      SE = numeric(),
      Lower_CI = numeric(),
      Upper_CI = numeric(),
      P_value = numeric(),
      stringsAsFactors = FALSE
    )
    
    # Run MAIC if selected
    if ("maic" %in% input$models) {
      maic_result <- run_maic(
        ipd_data(), 
        agg_data(), 
        input$outcome_var, 
        input$treatment_var, 
        input$covariates, 
        input$data_type, 
        input$outcome_scale
      )
      
      results <- rbind(results, data.frame(
        Method = maic_result$Method,
        Estimate = maic_result$Estimate,
        SE = maic_result$SE,
        Lower_CI = maic_result$Lower_CI,
        Upper_CI = maic_result$Upper_CI,
        P_value = maic_result$P_value
      ))
    }
    
    # Note: STC implementation will be added later
    if ("stc" %in% input$models) {
      results <- rbind(results, data.frame(
        Method = "STC (Not implemented yet)",
        Estimate = NA,
        SE = NA,
        Lower_CI = NA,
        Upper_CI = NA,
        P_value = NA
      ))
    }
    
    # Remove the notification when complete
    removeNotification(id = "analysis")
    
    return(results)
  })
  
  # Display results table
  output$results_table <- renderDT({
    req(analysis_results())
    
    results <- analysis_results()
    
    # Format the results table
    formatted_results <- results %>%
      mutate(
        `Estimate (95% CI)` = ifelse(
          is.na(Estimate),
          "Not implemented yet",
          paste0(
            sprintf("%.3f", Estimate), 
            " (", 
            sprintf("%.3f", Lower_CI), 
            " to ", 
            sprintf("%.3f", Upper_CI), 
            ")"
          )
        ),
        `P-value` = ifelse(is.na(P_value), "", sprintf("%.4f", P_value))
      ) %>%
      select(Method, `Estimate (95% CI)`, `P-value`)
    
    datatable(
      formatted_results,
      options = list(dom = 't'),
      rownames = FALSE,
      caption = paste0("Results using ", 
                      switch(input$outcome_scale,
                             "lor" = "Log Odds Ratio",
                             "rr" = "Risk Ratio"))
    )
  })
}

# Run the application
shinyApp(ui = ui, server = server)
