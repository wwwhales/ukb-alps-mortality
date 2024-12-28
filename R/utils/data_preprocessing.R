#' Data Preprocessing Functions for UK Biobank ALPS Mortality Analysis
#' 
#' This script contains functions for data cleaning and preprocessing

library(dplyr)
library(tidyr)

#' Read and clean ALPS mortality data
#' @param file_path Path to the raw data file
#' @return Cleaned dataset
read_clean_alps_data <- function(file_path) {
  alps_death <- read.csv(file_path)
  
  # Clean NA values for key variables
  alps_death1 <- alps_death[!is.na(alps_death$Sex) & 
                           !is.na(alps_death$Age) & 
                           !is.na(alps_death$Sleep) & 
                           !is.na(alps_death$college) & 
                           !is.na(alps_death$Townsand) & 
                           !is.na(alps_death$BMI) & 
                           !is.na(alps_death$Smoking) & 
                           !is.na(alps_death$Alcohal),]
  
  # Create derived variables
  alps_death1$friuts_vegetables <- with(alps_death1,
    ifelse(friuts < 5 | vegetables < 5, 1,
           ifelse(friuts > 8 | vegetables > 8, 3, 2)))
  
  return(alps_death1)
}

#' Prepare data for survival analysis
#' @param data Cleaned dataset
#' @return Dataset ready for survival analysis
prepare_survival_data <- function(data) {
  # Add any additional data preparation steps here
  return(data)
} 