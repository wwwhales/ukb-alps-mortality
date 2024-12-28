#' Table Generation Functions for ALPS Mortality Analysis
#' 
#' This script contains functions for generating tables

library(tableone)
library(dplyr)
library(tidyr)

#' Create Table 1 with baseline characteristics
#' @param data Dataset
#' @param strata_var Stratification variable
#' @param vars Variables to include
#' @param cat_vars Categorical variables
#' @return Table 1 object
create_table_one <- function(data, strata_var, vars, cat_vars) {
  tab1 <- CreateTableOne(
    vars = vars,
    strata = strata_var,
    data = data,
    factorVars = cat_vars,
    addOverall = TRUE
  )
  
  return(print(tab1, 
               catDigits = 1,
               contDigits = 1,
               pDigits = 3,
               showAllLevels = TRUE,
               quote = FALSE,
               noSpaces = TRUE))
}

#' Format Cox model results for publication
#' @param results Cox model results
#' @param exposure_order Order of exposure variables
#' @param outcome_order Order of outcome variables
#' @return Formatted results table
format_cox_results <- function(results, exposure_order, outcome_order) {
  results %>%
    mutate(
      HR_CI = paste0(round(HR, 2), " (", round(LCI, 2), "-", round(UCI, 2), ")"),
      p = case_when(
        p < 0.001 ~ "< .001",
        p >= 0.001 & p < 0.01 ~ sub("0", "", sprintf("%.3f", p)),
        p >= 0.01 ~ sub("^0?", "", sprintf("%.2f", p))
      )
    ) %>%
    select(exposure, outcome, HR_CI, p) %>%
    arrange(match(exposure, exposure_order)) %>%
    arrange(match(outcome, outcome_order)) %>%
    pivot_wider(
      names_from = exposure,
      values_from = c(HR_CI, p),
      names_vary = "slowest"
    )
}

#' Calculate person-years and event rates
#' @param data Dataset
#' @param disease_list List of diseases
#' @return Data frame with person-years and event rates
calculate_person_years <- function(data, disease_list) {
  results <- lapply(disease_list, function(x) {
    ind <- paste0(x, "_status")
    time <- paste0(x, "_days_ins2")
    
    data.frame(
      outcome = x,
      events_number = table(data[[ind]])[2],
      median_followup = median(data[[time]], na.rm = TRUE) / 365,
      person_years = sum(data[[time]], na.rm = TRUE) / 365
    )
  })
  
  do.call(rbind, results)
}

#' Format subgroup analysis results
#' @param results Subgroup analysis results
#' @param group_var Grouping variable
#' @return Formatted subgroup results
format_subgroup_results <- function(results, group_var) {
  results %>%
    mutate(
      HR_CI = paste0(round(HR, 2), " (", round(LCI, 2), "-", round(UCI, 2), ")"),
      p = ifelse(p < 0.001, "< .001", sprintf("%.3f", p))
    ) %>%
    select(subgroup, exposure, outcome, HR_CI, p) %>%
    pivot_wider(
      names_from = exposure,
      values_from = c(HR_CI, p),
      names_vary = "slowest"
    )
}

#' Format competing risk analysis results
#' @param results Competing risk analysis results
#' @return Formatted competing risk results
format_competing_risk_results <- function(results) {
  results %>%
    mutate(
      HR_CI = paste0(round(HR, 2), " (", round(LCI, 2), "-", round(UCI, 2), ")"),
      p = ifelse(p < 0.001, "< .001", sprintf("%.3f", p))
    ) %>%
    select(exposure, outcome, HR_CI, p) %>%
    pivot_wider(
      names_from = exposure,
      values_from = c(HR_CI, p),
      names_vary = "slowest"
    )
} 