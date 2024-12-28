#' Cox Proportional Hazards Models for ALPS Mortality Analysis
#' 
#' This script contains functions for running Cox regression analyses

library(survival)
library(survminer)

#' Run Cox regression analysis for death outcomes
#' @param x Predictor variables
#' @param y Outcome variables
#' @param t Time variable
#' @param cov1 Covariates (optional)
#' @param ID ID variable
#' @param d Dataset
#' @return Results data frame
cox2023_death <- function(x, y, t, cov1 = NULL, ID, d) {
  dat2 <- data.frame()
  
  for (i in 1:ncol(x)) {
    for (j in 1:ncol(y)) {
      nx <- names(x)[[i]]
      ny <- names(y)[[j]]
      nt <- names(t)
      id <- names(ID)
      
      # Create formula and fit model
      if (is.null(cov1)) {
        FML <- as.formula(paste0("Surv(", nt, ",", ny, "==1)~scale(", nx, ")"))
      } else {
        covariates <- cov1
        FML <- as.formula(paste0("Surv(", nt, ",", ny, "==1)~scale(", nx, ")+", 
                                paste0(covariates, collapse = "+")))
      }
      
      cox_a <- coxph(FML, data = d)
      cox_b <- summary(cox_a)
      
      # Extract results
      HR <- cox_b$conf.int[1,1]
      LCI <- cox_b$conf.int[1,3]
      UCI <- cox_b$conf.int[1,4]
      z <- cox_b$coefficients[1,"z"]
      p <- cox_b$coefficients[1,"Pr(>|z|)"]
      person_years <- sum(d[nt])
      obs <- cox_b$n
      event <- cox_b$nevent
      control <- cox_b$n - cox_b$nevent
      
      # Combine results
      current_results <- data.frame(
        exposure = nx,
        outcome = ny,
        cov = ifelse(is.null(covariates), "", paste0(covariates, collapse = "+")),
        HR = HR,
        LCI = LCI,
        UCI = UCI,
        z = z,
        p = p,
        person_years = person_years,
        obs = obs,
        event = event,
        control = control,
        stringsAsFactors = FALSE
      )
      
      dat2 <- rbind(dat2, current_results)
    }
  }
  
  return(dat2)
}

#' Run Cox regression analysis for incidence outcomes
#' @param x Predictor variables
#' @param y Outcome variable
#' @param cov1 Covariates (optional)
#' @param t Time variable
#' @param d Dataset
#' @return Results data frame
cox2023_incidence <- function(x, y, cov1 = NULL, t, d) {
  dat2 <- data.frame()
  
  for (i in 1:ncol(x)) {
    nx <- names(x)[[i]]
    ny <- names(y)
    nt <- names(t)
    
    if (is.null(cov1)) {
      FML <- as.formula(paste0("Surv(", nt, ",", ny, "==1)~scale(", nx, ")"))
    } else {
      covariates <- cov1
      FML <- as.formula(paste0("Surv(", nt, ",", ny, "==1)~scale(", nx, ")+",
                              paste0(covariates, collapse = "+")))
    }
    
    cox_a <- coxph(FML, data = d)
    cox_b <- summary(cox_a)
    
    # Extract results
    HR <- cox_b$conf.int[1,1]
    LCI <- cox_b$conf.int[1,3]
    UCI <- cox_b$conf.int[1,4]
    z <- cox_b$coefficients[1,"z"]
    p <- cox_b$coefficients[1,"Pr(>|z|)"]
    person_years <- sum(d[nt])
    obs <- cox_b$n
    event <- cox_b$nevent
    control <- cox_b$n - cox_b$nevent
    
    # Combine results
    current_results <- data.frame(
      exposure = nx,
      outcome = ny,
      HR = HR,
      LCI = LCI,
      UCI = UCI,
      z = z,
      p = p,
      person_years = person_years,
      obs = obs,
      event = event,
      control = control,
      cov = ifelse(is.null(covariates), "", paste0(covariates, collapse = "+")),
      stringsAsFactors = FALSE
    )
    
    dat2 <- rbind(dat2, current_results)
  }
  
  return(dat2)
}

#' Run competing risk analysis
#' @param data Dataset
#' @param time_var Time variable
#' @param event_var Event variable
#' @param predictor Predictor variable
#' @param covariates Vector of covariate names
#' @param competing_event Competing event indicator
#' @return Results of competing risk analysis
run_competing_risk_analysis <- function(data, time_var, event_var, predictor, 
                                      covariates = NULL, competing_event) {
  library(cmprsk)
  
  # Prepare formula
  if (is.null(covariates)) {
    formula_str <- paste0("Surv(", time_var, ",", event_var, ")~scale(", predictor, ")")
  } else {
    formula_str <- paste0("Surv(", time_var, ",", event_var, ")~scale(", predictor, ")+",
                         paste0(covariates, collapse = "+"))
  }
  
  # Fit competing risk model
  cph_model <- finegray(as.formula(formula_str), data = data, etype = competing_event)
  fit <- coxph(Surv(fgstart, fgstop, fgstatus) ~ ., data = cph_model, weights = fgwt)
  
  return(summary(fit))
} 