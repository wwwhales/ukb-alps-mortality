#' Visualization Functions for ALPS Mortality Analysis
#' 
#' This script contains functions for creating plots and figures

library(ggplot2)
library(survminer)

#' Create Kaplan-Meier survival curves
#' @param data Dataset
#' @param time_var Time variable
#' @param event_var Event variable
#' @param group_var Grouping variable
#' @return ggplot object
plot_km_curves <- function(data, time_var, event_var, group_var) {
  fit <- survfit(as.formula(paste0("Surv(", time_var, ",", event_var, ") ~", group_var)), 
                 data = data)
  
  ggsurvplot(fit,
             data = data,
             risk.table = TRUE,
             pval = TRUE,
             conf.int = TRUE,
             xlab = "Time in years",
             ylab = "Survival probability",
             risk.table.height = 0.25,
             ggtheme = theme_minimal())
}

#' Create forest plot for subgroup analyses
#' @param results Results data frame
#' @return ggplot object
create_forest_plot <- function(results) {
  # Add forest plot implementation
} 