### Regression Analysis of qPCR to determine primer efficiency
# author: Vera Laub
# last edited: 2025-02-03
# input: .xls of Ct values
# output: Excel file containing summary statistics

# Load libraries
library(readxl)
library(dplyr)
library(ggplot2)
library(writexl)

# Load data
data <- read_excel("~/Documents/postdoc/protocols/qPCR/2025-01-30_genomicPrimertrial_qPCR/2025-01-30_Primertrial_qPCR_1-36.xlsx", 
                   sheet = "results_mc")
data <- data %>% filter(`Sample Name` != "H2O")

# Set up cell content in variables for analysis
samples <- data$`Sample Name`
targets <- data$`Target Name`
ct_values <- data$CT

# Prepare an empty list to store results for each target
results_list <- list()

# Create an empty data frame to store all summary statistics
summary_df <- data.frame(Target = character(),
                         R_squared = numeric(),
                         Intercept = numeric(),
                         Intercept_p_value = numeric(),
                         Slope = numeric(),
                         Slope_p_value = numeric(),
                         stringsAsFactors = FALSE)

# Loop through each unique target
unique_targets <- unique(targets)
for (target in unique_targets) {
  
  # Filter data for the current target
  target_data <- data %>% filter(`Target Name` == target)
  
  # Extract ct_values and samples for the current target
  ct_values_target <- target_data$CT
  samples_target <- target_data$`Sample Name`
  
  # Perform linear regression
  lm_model <- lm(ct_values_target ~ samples_target)
  
  # Summarize the model
  model_summary <- summary(lm_model)
  
  # Extract the relevant statistics
  r_squared <- model_summary$r.squared
  intercept <- coef(lm_model)[1]
  slope <- coef(lm_model)[2]
  intercept_p_value <- model_summary$coefficients[1, 4]
  slope_p_value <- model_summary$coefficients[2, 4]
  
  # Add this information to the summary data frame
  summary_df <- rbind(summary_df, data.frame(
    Target = target,
    R_squared = r_squared,
    Intercept = intercept,
    Intercept_p_value = intercept_p_value,
    Slope = slope,
    Slope_p_value = slope_p_value
  ))
}

# Sort the data frame by Target
summary_df <- summary_df %>% arrange(Target)

# Define the output Excel file path
output_file <- "~/Documents/postdoc/protocols/qPCR/2025-01-30_genomicPrimertrial_qPCR/regression_summary.xlsx"

# Write the summary data to an Excel file
write_xlsx(summary_df, path = output_file)

# Print the summary dataframe to the console (optional)
print(summary_df)
