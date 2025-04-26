
##TODO: how does this compare to my function?
#
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
