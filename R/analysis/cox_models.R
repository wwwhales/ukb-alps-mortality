#' Cox Proportional Hazards Models for ALPS Mortality Analysis
#' 
#' This script contains functions for running Cox regression analyses

source("R/utils/helper_functions.R")

#' Run Cox regression analysis
#' @param x Predictor variables
#' @param y Outcome variables
#' @param t Time variable
#' @param cov1 Covariates
#' @param ID ID variable
#' @param d Dataset
#' @return Results data frame
cox2023_death <- function(x, y, t, cov1 = NULL, ID, d) {
  library(survival)
  library(survminer)
  
  dat2 <- data.frame()
  
  for (i in 1:ncol(x)) {
    for (j in 1:ncol(y)) {
      nx <- names(x)[[i]]
      ny <- names(y)[[j]]
      nt <- names(t)
      id <- names(ID)
      
      # Create formula and fit model
      FML <- create_cox_formula(nt, ny, nx, cov1)
      cox_a <- coxph(FML, data = d)
      
      # Extract results
      results <- extract_cox_results(cox_a)
      
      # Additional calculations
      person_years <- sum(d[nt])
      
      # Combine results
      current_results <- data.frame(
        predictor = nx,
        outcome = ny,
        results,
        person_years = person_years,
        stringsAsFactors = FALSE
      )
      
      dat2 <- rbind(dat2, current_results)
    }
  }
  
  return(dat2)
} 