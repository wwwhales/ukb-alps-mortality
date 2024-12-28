# UK Biobank ALPS Mortality Analysis

This repository contains code for analyzing the relationship between diffusion tensor image analysis along the perivascular space (DTI-ALPS) and mortality risk using UK Biobank data.

## Project Structure

- `R/`: R scripts for analysis
  - `analysis/`: Main analysis scripts
  - `utils/`: Helper functions and data preprocessing
- `data/`: Directory for input data (not tracked by git)
- `output/`: Directory for analysis outputs (not tracked by git)

## Key Files

### Analysis Scripts
- `baseline_analysis.R`: Baseline characteristics and Table 1 analysis
- `cox_models.R`: Cox proportional hazards models for mortality analysis
- `competing_risk.R`: Competing risk analysis
- `subgroup_analysis.R`: Subgroup and interaction analyses
- `visualization.R`: Plotting functions for figures

### Utility Scripts
- `data_preprocessing.R`: Data cleaning and preprocessing functions
- `helper_functions.R`: Common helper functions used across analyses

## Requirements

- R >= 4.0.0
- Required R packages:
  - survival
  - survminer
  - tableone
  - dplyr
  - tidyr
  - ggplot2
  - rms
  - forestplot
  - xlsx

## Usage

1. Place your UK Biobank data in the `data/` directory
2. Run scripts in the following order:
   ```R
   source("R/utils/data_preprocessing.R")
   source("R/analysis/baseline_analysis.R")
   source("R/analysis/cox_models.R")
   source("R/analysis/competing_risk.R")
   source("R/analysis/subgroup_analysis.R")
   source("R/analysis/visualization.R")
   ```

## Output

Analysis results will be saved in the `output/` directory, including:
- Tables in Excel format
- Forest plots
- Survival curves
- Other figures

## License

None

## Contact

None
