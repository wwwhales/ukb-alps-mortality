#' Main Analysis Script for ALPS Mortality Study
#' 
#' This script orchestrates the entire analysis workflow

# Load required packages and source files
library(survival)
library(survminer)
library(dplyr)
library(tidyr)
library(ggplot2)
library(tableone)
library(forestplot)
library(rms)
library(xlsx)

source("R/utils/data_preprocessing.R")
source("R/utils/helper_functions.R")
source("R/analysis/cox_models.R")
source("R/analysis/subgroup_analysis.R")
source("R/analysis/visualization.R")
source("R/analysis/table_generation.R")

# Set working directory and read data
data_path <- "data/alps_death.csv"
alps_data <- read_clean_alps_data(data_path)

# Create quartiles and quintiles
alps_data <- create_alps_quartiles(alps_data)
alps_data_quintiles <- create_alps_quintiles(alps_data)

# Define variables for analysis
predictor_vars <- c("DTI_ALPS_th60_mean", "DTI_ALPS_th60_L", "DTI_ALPS_th60_R")
outcome_vars <- c("all_death", "cancer_death", "cv_death", "digestive_death", 
                 "neuro_death", "resp_death")
covariates <- c("Age", "Sex", "Townsand", "BMI", "college", "Sleep", 
                "Smoking", "Alcohal", "ethinic")
disease_list <- c("neuro_disease", "cancer_disease", "cvd_disease", 
                 "digestive_disease", "respiratory_disease")

# Create Table 1
table1_vars <- c("Age", "Sex", "Sleep", "college", "Townsand", "BMI", 
                 "Smoking", "Alcohal", "ethinic")
cat_vars <- c("Sex", "college", "Sleep", "Smoking", "Alcohal", "ethinic")
baseline_table <- create_table_one(alps_data, "mean_levels", table1_vars, cat_vars)
write.csv(baseline_table, "output/table1.csv")

# Run main Cox analyses
# Model 1: Unadjusted
model0_results <- cox2023_death(
  x = alps_data[predictor_vars],
  y = alps_data[outcome_vars],
  t = alps_data["death_days_ins2"],
  ID = alps_data["eid"],
  d = alps_data
)

# Model 2: Age and sex adjusted
model1_results <- cox2023_death(
  x = alps_data[predictor_vars],
  y = alps_data[outcome_vars],
  t = alps_data["death_days_ins2"],
  cov1 = c("Age", "Sex"),
  ID = alps_data["eid"],
  d = alps_data
)

# Model 3: Fully adjusted
model2_results <- cox2023_death(
  x = alps_data[predictor_vars],
  y = alps_data[outcome_vars],
  t = alps_data["death_days_ins2"],
  cov1 = covariates,
  ID = alps_data["eid"],
  d = alps_data
)

# Run incidence analyses
incidence_results <- list()
for (disease in disease_list) {
  status_var <- paste0(disease, "_status")
  time_var <- paste0(disease, "_days_ins2")
  
  # Subset data for those without missing status
  subset_data <- alps_data[!is.na(alps_data[[status_var]]),]
  
  results <- cox2023_incidence(
    x = subset_data[predictor_vars],
    y = subset_data[status_var],
    t = subset_data[time_var],
    cov1 = covariates,
    d = subset_data
  )
  
  incidence_results[[disease]] <- results
}

# Run subgroup analyses
subgroup_vars <- c("Age", "Sex", "BMI", "Smoking")
subgroup_results <- perform_subgroup_analysis(
  data = alps_data,
  covariates = covariates,
  groups = subgroup_vars,
  predictor_vars = predictor_vars,
  outcome_vars = outcome_vars,
  time_var = "death_days_ins2"
)

# Create visualizations
# Kaplan-Meier curves
km_plot <- plot_km_curves(
  data = alps_data,
  time_var = "death_days_ins2",
  event_var = "all_death",
  group_var = "alps_levels1",
  title = "Survival by ALPS Quartiles"
)

# Forest plots
forest_plot <- create_forest_plot(
  results = model2_results,
  outcome_order = c("all_death", "cancer_death", "cv_death", 
                   "digestive_death", "neuro_death", "resp_death"),
  exposure_var = "DTI_ALPS_th60_mean"
)

# RCS curves
rcs_plot <- plot_rcs_curve(
  data = alps_data,
  var = "DTI_ALPS_th60_mean",
  time_var = "death_days_ins2",
  event_var = "all_death",
  covariates = covariates
)

# Correlation plot
corr_plot <- plot_alps_correlation(alps_data)

# Save results
save(model0_results, model1_results, model2_results, 
     incidence_results, subgroup_results,
     file = "output/analysis_results.RData")

# Save plots
ggsave("output/km_curves.pdf", km_plot)
ggsave("output/rcs_curves.pdf", rcs_plot)
ggsave("output/correlation.pdf", corr_plot)

# Format and save tables
formatted_cox_results <- format_cox_results(
  model2_results,
  exposure_order = predictor_vars,
  outcome_order = outcome_vars
)

formatted_subgroup_results <- format_subgroup_results(
  subgroup_results,
  group_var = subgroup_vars
)

write.xlsx(formatted_cox_results, 
           "output/cox_results.xlsx", 
           sheetName = "Main Results")

write.xlsx(formatted_subgroup_results, 
           "output/subgroup_results.xlsx", 
           sheetName = "Subgroup Analysis") 