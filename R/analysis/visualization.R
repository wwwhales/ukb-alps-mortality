#' Visualization Functions for ALPS Mortality Analysis
#' 
#' This script contains functions for creating plots and figures

library(ggplot2)
library(survminer)
library(forestplot)
library(patchwork)

#' Create Kaplan-Meier survival curves
#' @param data Dataset
#' @param time_var Time variable
#' @param event_var Event variable
#' @param group_var Grouping variable
#' @param title Plot title
#' @return ggplot object
plot_km_curves <- function(data, time_var, event_var, group_var, title = "") {
  fit <- survfit(as.formula(paste0("Surv(", time_var, ",", event_var, ") ~", group_var)), 
                 data = data)
  
  ggsurvplot(fit,
             data = data,
             risk.table = TRUE,
             pval = TRUE,
             conf.int = TRUE,
             xlab = "Time in years",
             ylab = "Survival probability",
             title = title,
             risk.table.height = 0.25,
             ggtheme = theme_minimal())
}

#' Create forest plot for mortality or incidence results
#' @param results Results data frame
#' @param outcome_order Vector specifying order of outcomes
#' @param exposure_var Name of exposure variable to plot
#' @return forestplot object
create_forest_plot <- function(results, outcome_order, exposure_var) {
  # Filter and order data
  plot_data <- results[results$exposure == exposure_var,]
  plot_data <- plot_data[match(outcome_order, plot_data$outcome),]
  
  # Format HR and CI
  plot_data$HR_CI <- paste0(round(plot_data$HR, 2), " (", 
                           round(plot_data$LCI, 2), "-", 
                           round(plot_data$UCI, 2), ")")
  
  # Create table text
  tabletext <- cbind(
    c("Disease (ICD-10 code)", "\n", plot_data$outcome),
    c("Total participants", "\n", plot_data$obs),
    c("HR (95% CI)", "\n", plot_data$HR_CI),
    c("P value", "\n", round(plot_data$p, 2))
  )
  
  # Create forest plot
  forestplot(
    labeltext = tabletext,
    mean = c(NA, NA, plot_data$HR),
    lower = c(NA, NA, plot_data$LCI),
    upper = c(NA, NA, plot_data$UCI),
    zero = 1.0,
    boxsize = 0.3,
    col = fpColors(box = "#374E55FF", zero = "grey", lines = "black"),
    xlab = "Hazard Ratio (95% CI)",
    txt_gp = fpTxtGp(label = gpar(cex = 1)),
    hrzl_lines = TRUE,
    lineheight = unit(1, "cm"),
    graph.pos = 4
  )
}

#' Create restricted cubic spline plot
#' @param data Dataset
#' @param var Variable for spline
#' @param time_var Time variable
#' @param event_var Event variable
#' @param covariates Vector of covariate names
#' @param knots Number of knots for spline
#' @return ggplot object
plot_rcs_curve <- function(data, var, time_var, event_var, covariates, knots = 3) {
  library(rms)
  
  # Set up datadist
  dd <- datadist(data)
  options(datadist = "dd")
  
  # Create formula
  formula_str <- paste0("Surv(", time_var, ",", event_var, "==1) ~ rcs(", var, 
                       ", ", knots, ") + ", paste0(covariates, collapse = "+"))
  
  # Fit model
  fit <- cph(as.formula(formula_str), data = data)
  
  # Create prediction plot
  pred <- Predict(fit, var, fun = exp)
  
  # Convert to ggplot
  ggplot(pred) +
    theme_classic() +
    scale_color_manual(values = c("#374E55B2")) +
    theme(axis.title = element_text(size = 15, color = "black")) +
    labs(x = paste(var, "value"), y = "HR") +
    theme(legend.position = "none", 
          axis.text = element_text(size = 12, color = "black")) +
    geom_hline(yintercept = 1, linetype = 2)
}

#' Create correlation plot between ALPS variables
#' @param data Dataset containing ALPS variables
#' @return ggplot object
plot_alps_correlation <- function(data) {
  # Calculate correlation
  correlation <- cor.test(data$DTI_ALPS_th60_L, data$DTI_ALPS_th60_R, 
                         method = "pearson")
  
  ggplot(data, aes(x = DTI_ALPS_th60_L, y = DTI_ALPS_th60_R)) +
    geom_point(size = 1, alpha = 0.7, shape = 21, 
               color = "black", fill = "#374E55B2") +
    stat_smooth(method = "lm", se = TRUE, color = "grey40") +
    theme_classic() +
    annotate("text", 
             x = min(data$DTI_ALPS_th60_L), 
             y = max(data$DTI_ALPS_th60_R),
             label = paste("Pearson r =", 
                          round(correlation$estimate, 2),
                          "\np < 2.2e-16"),
             hjust = 0, vjust = 1) +
    theme(axis.title = element_text(size = 15, color = "black")) +
    labs(x = "Left DTI-ALPS index", y = "Right DTI-ALPS index") +
    theme(legend.position = "none", 
          axis.text = element_text(size = 12, color = "black"))
}

#' Create quartile HR plot
#' @param data Dataset with quartile results
#' @param var_name Variable name for plotting
#' @return ggplot object
plot_quartile_hr <- function(data, var_name) {
  ggplot(data, aes(x = factor(alps_levels1), y = HR, fill = factor(alps_levels1))) +
    geom_point(aes(color = factor(alps_levels1)), 
               stat = "identity", pch = 21, size = 3) +
    geom_errorbar(aes(ymin = LCI, ymax = UCI, color = factor(alps_levels1)),
                 size = 0.75, width = 0.08) +
    scale_fill_manual(values = c("#374E55B2", "#DF8F44B2", "#00A1D5B2", "#B24745B2")) +
    scale_color_manual(values = c("#374E55B2", "#DF8F44B2", "#00A1D5B2", "#B24745B2")) +
    scale_x_discrete(labels = c("Q4", "Q3", "Q2", "Q1")) +
    theme_classic() +
    theme(axis.title = element_text(size = 15, color = "black")) +
    labs(x = paste("Quartiles of", var_name), 
         y = "HR of all-cause mortality") +
    theme(legend.position = "none", 
          axis.text = element_text(size = 12, color = "black"))
} 