#' Subgroup Analysis for ALPS Mortality Study
#' 
#' This script contains functions for conducting subgroup analyses

library(survival)
library(dplyr)

#' Perform subgroup analysis for categorical and continuous variables
#' @param data Main dataset
#' @param covariates Vector of covariate names
#' @param groups Vector of grouping variables
#' @param predictor_vars Vector of predictor variables
#' @param outcome_vars Vector of outcome variables
#' @param time_var Time variable name
#' @return List of subgroup results
perform_subgroup_analysis <- function(data, covariates, groups, predictor_vars, 
                                    outcome_vars, time_var) {
  subgroup_results <- data.frame()
  
  for (grp in groups) {
    # Handle categorical variables
    if (is.factor(data[[grp]])) {
      levels <- unique(data[[grp]])
      for (level in levels) {
        # Create subset
        subset_data <- data[data[[grp]] == level,]
        
        # Adjust covariates (remove current grouping variable)
        adj_covariates <- setdiff(covariates, grp)
        
        # Run Cox model
        results <- cox2023_death(
          x = subset_data[predictor_vars],
          y = subset_data[outcome_vars],
          t = subset_data[time_var],
          cov1 = adj_covariates,
          ID = subset_data["eid"],
          d = subset_data
        )
        
        # Add subgroup information
        results$subgroup <- paste0(grp, "_", level)
        subgroup_results <- rbind(subgroup_results, results)
      }
    } else {
      # Handle continuous variables (split by median)
      median_val <- median(data[[grp]], na.rm = TRUE)
      
      # Create subsets
      subset_low <- data[data[[grp]] <= median_val,]
      subset_high <- data[data[[grp]] > median_val,]
      
      # Adjust covariates
      adj_covariates <- setdiff(covariates, grp)
      
      # Run Cox models for both groups
      results_low <- cox2023_death(
        x = subset_low[predictor_vars],
        y = subset_low[outcome_vars],
        t = subset_low[time_var],
        cov1 = adj_covariates,
        ID = subset_low["eid"],
        d = subset_low
      )
      results_low$subgroup <- paste0(grp, "_<=", round(median_val, 2))
      
      results_high <- cox2023_death(
        x = subset_high[predictor_vars],
        y = subset_high[outcome_vars],
        t = subset_high[time_var],
        cov1 = adj_covariates,
        ID = subset_high["eid"],
        d = subset_high
      )
      results_high$subgroup <- paste0(grp, "_>", round(median_val, 2))
      
      subgroup_results <- rbind(subgroup_results, results_low, results_high)
    }
  }
  
  return(subgroup_results)
}

#' Calculate interaction terms for subgroup analysis
#' @param data Dataset
#' @param predictors Vector of predictor variables
#' @param covariates Vector of covariate names
#' @param interaction_vars Vector of variables to test interactions with
#' @return Data frame of interaction test results
test_interactions <- function(data, predictors, covariates, interaction_vars) {
  results <- data.frame(
    predictor = character(),
    interaction_var = character(),
    adj_covariates = character(),
    p_value = numeric()
  )
  
  for (pred in predictors) {
    for (var in interaction_vars) {
      # Remove interaction variable from covariates
      adj_covs <- setdiff(covariates, var)
      
      # Create formula with interaction term
      formula <- as.formula(paste0(
        "Surv(death_days_ins2, all_death==1) ~ scale(", pred, ") * ", var,
        " + ", paste0(adj_covs, collapse = "+")
      ))
      
      # Fit model
      fit <- coxph(formula, data = data)
      
      # Extract interaction p-value
      p_value <- summary(fit)$coefficients[grep(":", rownames(summary(fit)$coefficients)), 5]
      
      # Store results
      current_result <- data.frame(
        predictor = pred,
        interaction_var = var,
        adj_covariates = paste(adj_covs, collapse = "+"),
        p_value = p_value
      )
      
      results <- rbind(results, current_result)
    }
  }
  
  return(results)
} 