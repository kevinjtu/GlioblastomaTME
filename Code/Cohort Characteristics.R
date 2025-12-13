#libraries####
library(data.table)
library(readxl)
library(ggplot2)
library(dplyr)
library(tidyr)

# Load the cell type data####
setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME")
combined_standard <- fread("GBM Clinical-TME Cell States.csv")

# Assuming your data is stored in `combined_standard` as a data.table
dt <- as.data.table(combined_standard)

# Define total number of rows (excluding NA rows per variable later)
total_n <- nrow(dt)

# Bin age into 10-year intervals
dt[, age_bin := cut(age, breaks = seq(0, 100, by = 10), right = FALSE, include.lowest = TRUE)]

# Define variables to summarize
vars_to_summarize <- c("age_bin", "sex", "race", "recurrence", 
                       "KarnPerfScore", "MGMT_status", "IDH_status", "platform")

# Function to generate % and count strings
summarize_var <- function(data, var) {
  var_data <- data[[var]]
  var_dt <- data[!is.na(var_data), .N, by = var]
  var_total <- sum(var_dt$N)
  var_dt[, pct := round((N / var_total) * 100, 1)]
  var_dt[, summary := sprintf("%s%% (%d/%d)", pct, N, var_total)]
  var_dt[, variable := var]
  setnames(var_dt, var, "value")
  var_dt[, .(variable, value, summary)]
}

# Apply summary function to each variable and bind results
summary_table <- rbindlist(lapply(vars_to_summarize, function(v) summarize_var(dt, v)))

# Print the summary table
print(summary_table)

#export ot excel
write.csv(summary_table, "Cohort_Characteristics_Summary.csv", row.names = FALSE)
