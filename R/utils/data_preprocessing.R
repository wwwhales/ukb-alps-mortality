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

#' Create quartile categories for ALPS variables
#' @param data Dataset containing ALPS variables
#' @return Dataset with added quartile categories
create_alps_quartiles <- function(data) {
  data$alps_levels1 <- cut(data$DTI_ALPS_th60_mean, 
                          breaks = quantile(data$DTI_ALPS_th60_mean),
                          labels = c("1", "2", "3", "4"))
  data$alps_levels1[is.na(data$alps_levels1)] <- "1"
  data$alps_levels1 <- factor(data$alps_levels1, levels = c(4,3,2,1))
  
  # Create L and R quartiles
  data$L_levels <- cut(data$DTI_ALPS_th60_L, 
                      breaks = quantile(data$DTI_ALPS_th60_L),
                      labels = c("1", "2", "3", "4"))
  data$L_levels[is.na(data$L_levels)] <- "1"
  data$L_levels <- factor(data$L_levels, levels = c(4,3,2,1))
  
  data$R_levels <- cut(data$DTI_ALPS_th60_R, 
                      breaks = quantile(data$DTI_ALPS_th60_R),
                      labels = c("1", "2", "3", "4"))
  data$R_levels[is.na(data$R_levels)] <- "1"
  data$R_levels <- factor(data$R_levels, levels = c(4,3,2,1))
  
  return(data)
}

#' Create quintile categories for ALPS variables
#' @param data Dataset containing ALPS variables
#' @return Dataset with added quintile categories
create_alps_quintiles <- function(data) {
  breaks <- c(0, 0.2, 0.4, 0.6, 0.8, 1)
  labels <- c("1", "2", "3", "4", "5")
  
  data$mean_levels <- cut(data$DTI_ALPS_th60_mean, 
                         breaks = quantile(data$DTI_ALPS_th60_mean, breaks),
                         labels = labels)
  data$mean_levels[is.na(data$mean_levels)] <- "1"
  data$mean_levels <- factor(data$mean_levels, levels = rev(labels))
  
  data$L_levels <- cut(data$DTI_ALPS_th60_L, 
                      breaks = quantile(data$DTI_ALPS_th60_L, breaks),
                      labels = labels)
  data$L_levels[is.na(data$L_levels)] <- "1"
  data$L_levels <- factor(data$L_levels, levels = rev(labels))
  
  data$R_levels <- cut(data$DTI_ALPS_th60_R, 
                      breaks = quantile(data$DTI_ALPS_th60_R, breaks),
                      labels = labels)
  data$R_levels[is.na(data$R_levels)] <- "1"
  data$R_levels <- factor(data$R_levels, levels = rev(labels))
  
  return(data)
}

#' Prepare data for competing risk analysis
#' @param data Input dataset
#' @return Dataset prepared for competing risk analysis
prepare_competing_risk_data <- function(data) {
  data %>% mutate(
    cvd_death = ifelse(neuro_death == 1, 2, cv_death)
  ) %>%
    mutate(cvd_death = as.factor(cvd_death))
} 