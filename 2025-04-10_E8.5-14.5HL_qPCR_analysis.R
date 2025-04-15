### Analysis of RT-qPCR to determine expression levels relative to control
# author: Vera Laub
# last edited: 2025-04-14
# input: .xlsx of Ct values
# output: Excel file containing summary statistics and plots

# Load required packages
library(readxl)
library(dplyr)
library(tidyr)

# ---- READ DATA ----
qpcr_data <- read_excel("Documents/postdoc/results/pioneering_activity_project/qPCR/2025-04-10_E8.5-14.5HLs_qPCR/2025-04-10_E8.5-14.5HL_qPCR_Vera_forR.xlsx", 
                        sheet = "results_for_R")
qpcr_data <- qpcr_data %>% filter(`Sample Name` != "H2O")
ref_gene <- "Actb"                  # Housekeeping gene
ref_sample <- "Control"             # Reference sample

# Define output folder
output <- "~/Documents/postdoc/results/pioneering_activity_project/qPCR/2025-04-10_E8.5-14.5HLs_qPCR/analysis/"



