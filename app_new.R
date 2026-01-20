library(shiny)
library(bslib)
library(DT)
library(dplyr)
library(sandwich)
library(outstandR)
library(glue)
library(tidyr)   # Required for reshaping
library(stringr) # Required for reshaping

# --- Helper Function for Wide Data Conversion ---
# Included directly to ensure availability even if internal to the package
reshape_ald_to_long_helper <- function(df) {
  # Logic adapted from outstandR::reshape_ald_to_long
  
  # Separate the 'colname' into 'statistic', 'variable', and 'treatment' columns
  # Assumes columns not starting with y or N are simple statistics (e.g. mean.age)
  variable_cols <- df |> 
    dplyr::select(-starts_with("y"), -starts_with("N")) |> 
    tidyr::pivot_longer(cols = everything(),
                        names_to = "colname",
                        values_to = "value") |> 
    tidyr::separate(.data$colname,
                    into = c("statistic", "variable"), sep = "\\.") |> 
    dplyr::select(.data$variable, .data$statistic, .data$value) |> 
    dplyr::mutate(trt = NA)  # Set trt to NA for covariate stats
  
  # N columns: Expected format N.Trt (e.g., N.B)
  N_cols <- df |> 
    dplyr::select(starts_with("N")) |> 
    tidyr::pivot_longer(cols = everything(),
                        names_to = "colname",
                        values_to = "value") |> 
    tidyr::separate(.data$colname,
                    into = c("statistic", "trt"), sep = "\\.") |> 
    dplyr::select(.data$statistic, .data$trt, .data$value) |> 
    dplyr::mutate(variable = NA)
  
  # Y columns: Expected format y.Trt.Stat (e.g., y.B.sum or y.B.mean)
  # Note: The provided file logic separates into: variable, trt, statistic
  # Assuming 'y' is the variable name here based on file snippet: y.B.sum -> var:y, trt:B, stat:sum
  y_cols <- df |> 
    dplyr::select(starts_with("y")) |> 
    tidyr::pivot_longer(cols = everything(),
                        names_to = "colname",
                        values_to = "value") |> 
    tidyr::separate(.data$colname,
                    into = c("variable", "trt", "statistic"),
                    sep = "\\.") |> 
    dplyr::select(.data$variable, .data$statistic, .data$trt, .data$value) 
  
  dplyr::bind_rows(variable_cols, y_cols, N_cols) |> 
    dplyr::arrange(.data$variable, .data$statistic, .data$trt)
}


# --- UI ---
ui <- page_navbar(
  title = "Indirect Treatment Comparison Tool (Pro)",
  
  # Data Input Tab
  nav_panel(
    title = "Data Input & Analysis",
    page_sidebar(
      sidebar = sidebar(
        width = 400,
        
        # File inputs
        card(
          card_header("Upload Data"),
          fileInput("ipd_file", "Upload Individual Patient Data (CSV)", accept = c(".csv")),
          fileInput("agg_file", "Upload Aggregate Data (CSV)", accept = c(".csv")),
          
          # New: Aggregate Data Format Toggle
          radioButtons("ald_format", "Aggregate Data Format:",
                       choices = c("Long (Standard)" = "long", 
                                   "Wide (e.g. mean.age, N.B)" = "wide"),
                       selected = "long")
        ),
        
        # Data type specification
        card(
          card_header("Specify Data Type"),
          radioButtons("data_type", "Data Type:",
                       choices = c("Binary", "Continuous", "Count"),
                       selected = "Binary")
        ),
        
        # Model selection
        card(
          card_header("Select Models to Fit"),
          checkboxGroupInput("models", "Models:",
                             choices = c("MAIC" = "maic", 
                                         "STC" = "stc",
                                         "G-Comp (ML)" = "gcomp_ml",
                                         "G-Comp (Bayes)" = "gcomp_bayes",
                                         "MIM (Multiple Imputation)" = "mim"), # Added MIM
                             selected = "maic")
        ),
        
        # Simulation Settings
        conditionalPanel(
          condition = "input.models.includes('gcomp_ml') || input.models.includes('gcomp_bayes') || input.models.includes('mim')",
          card(
            card_header("Simulation Settings"),
            numericInput("sim_n", "Simulation/Synthesis Sample Size (N):", value = 1000, min = 100, step = 100)
          )
        ),
        
        # Outcome scale selection
        card(
          card_header("Outcome Scale"),
          radioButtons("outcome_scale", "Select Outcome Scale:",
                       choices = c("Log Odds Ratio" = "log_odds", 
                                   "Risk Ratio" = "log_relative_risk",
                                   "Risk Difference" = "risk_difference",
                                   "Mean Difference" = "mean_difference",
                                   "Rate Difference" = "rate_difference"),
                       selected = "log_odds")
        ),
        
        # Variables selection UI
        uiOutput("variable_selection"),
        
        actionButton("run_analysis", "Run Analysis", class = "btn-primary")
      ),
      
      # Main panel
      card(
        card_header("Data Preview"),
        tabsetPanel(
          tabPanel("Individual Patient Data", DTOutput("ipd_preview")),
          tabPanel("Aggregate Data (Processed)", DTOutput("agg_preview"))
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

# --- Server ---
server <- function(input, output, session) {
  
  ipd_data <- reactiveVal(NULL)
  agg_data <- reactiveVal(NULL)
  
  # Process IPD
  observeEvent(input$ipd_file, {
    req(input$ipd_file)
    tryCatch({
      data <- read.csv(input$ipd_file$datapath, stringsAsFactors = FALSE)
      ipd_data(data)
    }, error = function(e) showNotification(paste("Error reading IPD:", e$message), type = "error"))
  })
  
  # Process Aggregate Data (with Reshape Logic)
  observeEvent(c(input$agg_file, input$ald_format), {
    req(input$agg_file)
    tryCatch({
      data <- read.csv(input$agg_file$datapath, stringsAsFactors = FALSE)
      
      # Handle Wide format
      if (input$ald_format == "wide") {
        data <- reshape_ald_to_long_helper(data)
        showNotification("Reshaped Wide data to Long format.", type = "message")
      }
      
      agg_data(data)
    }, error = function(e) showNotification(paste("Error reading Aggregate file:", e$message), type = "error"))
  })
  
  # Dynamic Variable Selection
  output$variable_selection <- renderUI({
    req(ipd_data())
    ipd <- ipd_data()
    card(
      card_header("Variable Selection"),
      selectInput("outcome_var", "Outcome Variable:", choices = names(ipd)),
      selectInput("treatment_var", "Treatment Variable:", choices = names(ipd)),
      checkboxGroupInput("progfactors", "Prognostic variables:", choices = names(ipd)),
      checkboxGroupInput("effmodifier", "Effect modifiers:", choices = names(ipd))
    )
  })
  
  output$ipd_preview <- renderDT({ req(ipd_data()); datatable(ipd_data(), options = list(scrollX = TRUE)) })
  output$agg_preview <- renderDT({ req(agg_data()); datatable(agg_data(), options = list(scrollX = TRUE)) })
  
  # Run Analysis
  analysis_results <- eventReactive(input$run_analysis, {
    req(ipd_data(), agg_data(), input$data_type, input$outcome_var, input$treatment_var, length(input$models) > 0)
    
    showNotification("Running analysis...", type = "message", id = "analysis")
    
    results_df <- data.frame(
      Method = character(), Treatments = character(), Estimate = numeric(),
      Std.Error = numeric(), lower.0.95 = numeric(), upper.0.95 = numeric(),
      stringsAsFactors = FALSE
    )
    
    # Family setup
    ffamily <- switch(input$data_type,
                      "Binary" = binomial(),
                      "Continuous" = gaussian(),
                      "Count" = poisson(),
                      stop("Unsupported data_type"))
    force(ffamily) 
    
    # Formula setup
    progfactors <- input$progfactors %||% ""
    effmodifier <- input$effmodifier %||% ""
    interaction_term <- if (length(effmodifier) > 0 && any(effmodifier != "")) {
      paste0(input$treatment_var, ":", effmodifier)
    } else { "" }
    
    formula_parts <- c(progfactors, input$treatment_var, interaction_term)
    rhs <- paste(formula_parts[formula_parts != ""], collapse = " + ")
    form <- as.formula(glue::glue("{input$outcome_var} ~ {rhs}"))
    
    # Helper for results
    append_res <- function(df, method_name, res_obj) {
      contrasts <- res_obj$results$contrasts
      rbind(df, data.frame(
        Method = method_name,
        Treatments = names(contrasts$means),
        Estimate = unlist(contrasts$means),
        Std.Error = sqrt(unlist(contrasts$variances)),
        lower.0.95 = sapply(contrasts$CI, \(x) x[1]),
        upper.0.95 = sapply(contrasts$CI, \(x) x[2])
      ))
    }
    
    # --- MAIC ---
    if ("maic" %in% input$models) {
      tryCatch({
        strat <- strategy_maic(formula = form, family = ffamily, trt_var = input$treatment_var)
        res <- outstandR::outstandR(ipd_trial = ipd_data(), ald_trial = agg_data(), strategy = strat, scale = input$outcome_scale)
        results_df <- append_res(results_df, "MAIC", res)
      }, error = function(e) showNotification(paste("MAIC Error:", e$message), type = "error"))
    }
    
    # --- STC ---
    if ("stc" %in% input$models) {
      tryCatch({
        strat <- strategy_stc(formula = form, family = ffamily, trt_var = input$treatment_var)
        res <- outstandR::outstandR(ipd_trial = ipd_data(), ald_trial = agg_data(), strategy = strat, scale = input$outcome_scale)
        results_df <- append_res(results_df, "STC", res)
      }, error = function(e) showNotification(paste("STC Error:", e$message), type = "error"))
    }
    
    # --- G-Comp ML ---
    if ("gcomp_ml" %in% input$models) {
      tryCatch({
        strat <- strategy_gcomp_ml(formula = form, family = ffamily, trt_var = input$treatment_var, N = as.integer(input$sim_n))
        res <- outstandR::outstandR(ipd_trial = ipd_data(), ald_trial = agg_data(), strategy = strat, scale = input$outcome_scale)
        results_df <- append_res(results_df, "G-Comp (ML)", res)
      }, error = function(e) showNotification(paste("G-Comp ML Error:", e$message), type = "error"))
    }
    
    # --- G-Comp Bayes ---
    if ("gcomp_bayes" %in% input$models) {
      tryCatch({
        showNotification("Running G-Comp Bayes...", type = "message", id = "bayes")
        strat <- strategy_gcomp_bayes(formula = form, family = ffamily, trt_var = input$treatment_var, N = as.integer(input$sim_n))
        res <- outstandR::outstandR(ipd_trial = ipd_data(), ald_trial = agg_data(), strategy = strat, scale = input$outcome_scale)
        results_df <- append_res(results_df, "G-Comp (Bayes)", res)
        removeNotification(id = "bayes")
      }, error = function(e) showNotification(paste("G-Comp Bayes Error:", e$message), type = "error"))
    }
    
    # --- MIM (Multiple Imputation) ---
    if ("mim" %in% input$models) {
      tryCatch({
        showNotification("Running MIM...", type = "message", id = "mim")
        # strategy_mim takes similar args to gcomp, N is number of simulations
        strat <- strategy_mim(formula = form, family = ffamily, trt_var = input$treatment_var, N = as.integer(input$sim_n))
        res <- outstandR::outstandR(ipd_trial = ipd_data(), ald_trial = agg_data(), strategy = strat, scale = input$outcome_scale)
        results_df <- append_res(results_df, "MIM", res)
        removeNotification(id = "mim")
      }, error = function(e) showNotification(paste("MIM Error:", e$message), type = "error"))
    }
    
    removeNotification(id = "analysis")
    return(results_df)
  })
  
  output$results_table <- renderDT({
    req(analysis_results())
    results <- analysis_results()
    if (nrow(results) == 0) return(NULL)
    
    formatted <- results %>%
      mutate(
        `Estimate (95% CI)` = paste0(sprintf("%.3f", Estimate), " (", sprintf("%.3f", lower.0.95), " to ", sprintf("%.3f", upper.0.95), ")"),
        `Std. Error` = sprintf("%.3f", Std.Error)
      ) %>%
      select(Method, Treatments, `Estimate (95% CI)`, `Std. Error`)
    
    datatable(formatted, options = list(dom = 't'), rownames = FALSE, caption = paste0("Results scale: ", input$outcome_scale))
  })
}

shinyApp(ui, server)