#libraries####
library(data.table)
library(Biobase)
library(dplyr)
library(textshape)
library(tibble)
library(readxl)
library(survival)
library(gridExtra)
library(ggplot2)
library(rms)
library(forestmodel)
library(stringr)
library(patchwork)

surv_metric <- "OS" # options are OS, DFS

# Load the cell type data####
setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME")
combined_standard <- read_excel("GBM TME Cell Type-Clinical Data.xlsx")
print(names(combined_standard))
new_names <- names(combined_standard)
new_names <- str_replace_all(new_names, "\\-", "")       # Remove hyphens
new_names <- str_replace_all(new_names, "\\/", "_")       # Remove slashes
new_names <- str_replace_all(new_names, " ", "_")         # Replace spaces with underscores
new_names <- str_replace_all(new_names, "\\+", "")         # Remove plus signs
combined_standard <- setNames(combined_standard, new_names)
print(names(combined_standard))

combined_standard <- combined_standard %>%
  filter(!is.na(Subtype))
combined_standard <- combined_standard %>% mutate_at(vars(7:10), as.numeric) #set survival data to numeric

#genearte kaplan meier curves
cell_types <- unique(colnames(combined_standard[, 14:23]))
subtypes <- unique(combined_standard$Subtype)

for (cell_type in cell_types) {
  # Add quartile classification for the current cell type
  combined_standard <- combined_standard %>%
    mutate(Quartile = ntile(!!sym(cell_type), 2)) %>%
    drop_na(OS_years, OS_event, Quartile, Subtype)
  
  plots <- list()  # Initialize a list to store plots
  
  for (subtype in subtypes) {
    # Subset data for the current subtype
    subset_data <- combined_standard %>% filter(Subtype == subtype)
    
    if (nrow(subset_data) < 2) next  # Skip if there isn't enough data to perform survival analysis
    
    # Create the survival object
    surv_obj <- Surv(subset_data$OS_years, subset_data$OS_event)
    
    # Fit Kaplan-Meier survival curves
    fit <- survfit(surv_obj ~ Quartile, data = subset_data)
    
    # Compute log-rank p-value
    logrank <- survdiff(surv_obj ~ Quartile, data = subset_data)
    p_value <- logrank$pvalue
    
    # Convert survival object to data frame for ggplot
    surv_data <- as_tibble(surv_summary(fit, data = subset_data))
    
    # Generate Kaplan-Meier plot
    plot <- ggplot(surv_data, aes(x = time, y = surv, color = factor(Quartile))) +
      geom_step() +
      scale_color_brewer(palette = "Set1") +
      labs(
        title = paste("Subtype:", subtype),
        subtitle = paste("Log-rank p =", signif(p_value, 3)),
        x = "Time (years)",
        y = "Survival Probability",
        color = "Quartile"
      ) +
      theme_minimal() +
      theme(legend.position = "bottom")
    
    # Add plot to the list
    plots[[subtype]] <- plot
  }
  
  # Combine plots for all subtypes for the current cell type
  if (length(plots) > 0) {
    combined_plot <- wrap_plots(plots, ncol = 3) +
      plot_annotation(title = paste("Kaplan-Meier Curves for", cell_type))
    
    # Display the combined plot
    print(combined_plot)
  }
}
