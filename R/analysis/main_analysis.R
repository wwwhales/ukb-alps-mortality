#' Main Analysis Script for ALPS Mortality Study
#' 
#' This script orchestrates the entire analysis workflow

# Load required packages and source files
source("R/utils/data_preprocessing.R")
source("R/utils/helper_functions.R")
source("R/analysis/cox_models.R")
source("R/analysis/subgroup_analysis.R")
source("R/analysis/visualization.R")

# Set working directory and read data
data_path <- "data/alps_death.csv"
alps_data <- read_clean_alps_data(data_path)

# Prepare data for analysis
analysis_data <- prepare_survival_data(alps_data)

# Run main Cox analysis
cox_results <- cox2023_death(
  x = analysis_data[predictor_vars],
  y = analysis_data[outcome_vars],
  t = analysis_data[time_var],
  cov1 = covariates,
  ID = analysis_data$id,
  d = analysis_data
)

# Run subgroup analyses
subgroup_results <- perform_subgroup_analysis(
  data = analysis_data,
  group_var = "Age",
  continuous = TRUE
)

# Create visualizations
km_plot <- plot_km_curves(
  data = analysis_data,
  time_var = "follow_up_time",
  event_var = "death",
  group_var = "alps_group"
)

# Save results
save(cox_results, subgroup_results, file = "output/analysis_results.RData")
ggsave("output/survival_curves.pdf", km_plot) 