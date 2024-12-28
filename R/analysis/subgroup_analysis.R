#' Subgroup Analysis for ALPS Mortality Study
#' 
#' This script contains functions for conducting subgroup analyses

source("R/utils/helper_functions.R")

#' Perform subgroup analysis
#' @param data Main dataset
#' @param group_var Grouping variable
#' @param continuous Logical indicating if group_var is continuous
#' @return List of subgroup results
perform_subgroup_analysis <- function(data, group_var, continuous = FALSE) {
  if (!continuous) {
    # Categorical subgroup analysis
    unique_groups <- unique(data[[group_var]])
    subgroup_results <- lapply(unique_groups, function(g) {
      subset_data <- data[data[[group_var]] == g,]
      # Add your analysis for each subgroup
    })
  } else {
    # Continuous variable - split by median
    median_val <- median(data[[group_var]], na.rm = TRUE)
    subgroup1 <- data[data[[group_var]] <= median_val,]
    subgroup2 <- data[data[[group_var]] > median_val,]
    # Add your analysis for each subgroup
  }
  
  return(subgroup_results)
} 