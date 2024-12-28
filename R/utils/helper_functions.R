#' Helper Functions for UK Biobank ALPS Mortality Analysis
#' 
#' This script contains common utility functions used across analyses

library(survival)
library(survminer)

#' Create formula for Cox models
#' @param time_var Name of time variable
#' @param event_var Name of event variable
#' @param predictor Predictor variable
#' @param covariates Vector of covariate names
#' @return Formula object
create_cox_formula <- function(time_var, event_var, predictor, covariates = NULL) {
  if (is.null(covariates)) {
    formula_str <- paste0("Surv(", time_var, ",", event_var, "==1)~scale(", predictor, ")")
  } else {
    formula_str <- paste0("Surv(", time_var, ",", event_var, "==1)~scale(", 
                         predictor, ")+", paste0(covariates, collapse = "+"))
  }
  return(as.formula(formula_str))
}

#' Extract Cox model results
#' @param model Fitted Cox model
#' @return Data frame with HR, CI, and p-value
extract_cox_results <- function(model) {
  summary_model <- summary(model)
  
  data.frame(
    HR = summary_model$conf.int[1,1],
    LCI = summary_model$conf.int[1,3],
    UCI = summary_model$conf.int[1,4],
    z = summary_model$coefficients[1,"z"],
    p_value = summary_model$coefficients[1,"Pr(>|z|)"],
    n_obs = summary_model$n
  )
} 