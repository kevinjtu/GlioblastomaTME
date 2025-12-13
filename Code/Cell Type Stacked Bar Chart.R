#libraries####
library(data.table)
library(readxl)
library(ggplot2)
library(dplyr)
library(tidyr)
library(GGally)

# Load the cell type data####
setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME")
combined_standard <- read_excel("GBM TME Cell Type-Clinical Data.xlsx")
print(names(combined_standard))
new_names <- names(combined_standard)

#plot stacked bar chart####
# Define cell type columns
cell_type_columns <- c("Oligodendrocyte", "Pericyte", "Neuron", "Astrocyte", 
                       "Endothelial", "Myeloid", "Lymphoid", "OPC", 
                       "radial glial cell", "Malignant")

# Separate Malignant and recalculate proportions
combined_standard <- combined_standard %>%
  arrange(desc(Malignant)) %>% # Order by Malignant
  mutate(sample_index = row_number(), # Add numeric index
         non_malignant_total = rowSums(across(all_of(setdiff(cell_type_columns, "Malignant")))),
         across(all_of(setdiff(cell_type_columns, "Malignant")), 
                ~ . / non_malignant_total, .names = "scaled_{col}")) # Rescale other proportions

# Prepare data for the main plot (non-malignant cell types)
main_data <- combined_standard %>%
  select(sample_index, starts_with("scaled_")) %>%
  pivot_longer(cols = starts_with("scaled_"), 
               names_to = "Cell_Type", 
               values_to = "Proportion") %>%
  mutate(Cell_Type = gsub("scaled_", "", Cell_Type)) # Clean cell type names

# Prepare data for the Malignant bar plot
malignant_data <- combined_standard %>%
  select(sample_index, Malignant)

# Main stacked bar chart (non-malignant populations)
main_plot <- ggplot(main_data, aes(x = sample_index, y = Proportion, fill = Cell_Type)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Proportion of Non-Malignant Cell Types",
       x = NULL,
       y = "Recalculated Proportion",
       fill = "Cell Type") +
  theme(axis.text.x = element_blank(), # Hide x-axis labels
        axis.ticks.x = element_blank())

# Malignant bar chart
malignant_plot <- ggplot(malignant_data, aes(x = sample_index, y = Malignant)) +
  geom_bar(stat = "identity", fill = "red") +
  theme_minimal() +
  labs(title = "Proportion of Malignant Cell Types",
       x = "Sample Index (Ordered by Malignant)",
       y = "Proportion") +
  theme(axis.text.x = element_blank(), # Hide x-axis labels
        axis.ticks.x = element_blank())

# Combine plots
combined_plot <- main_plot / malignant_plot + plot_layout(heights = c(3, 1)) # Main plot larger

# Display the combined plot
print(combined_plot)


