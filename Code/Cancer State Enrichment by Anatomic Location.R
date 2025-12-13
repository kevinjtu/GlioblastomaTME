#libraries####
library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(RColorBrewer)
library(stringr)

#import & clean data####
subtype_cols <- c("OPC_like_Cancer", "AC_like_Cancer", "NPC_like_Cancer", "MES_like_Cancer")
setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME")
data <- as.data.frame(fread("GBM Clinical-TME Cell States.csv")) %>%
  mutate(subtype = subtype_cols[
    max.col( as.matrix( across(all_of(subtype_cols)) ),
             ties.method = "first")]) %>%
  filter(study == "IVYGAP")

#get TME fractions accounting for cellularity
# define your six tumor‐cell columns
tumor_cells <- c(
  "OPC_like_Cancer",
  "AC_like_Cancer",
  "NPC_like_Cancer",
  "MES_like_Cancer"
)

data2 <- data %>%
  # 1) compute tumor purity as the row‐sum of those six columns
  rowwise() %>%
  mutate(tumor_purity = sum(c_across(all_of(tumor_cells)))) %>%
  ungroup() %>%
  
  # 2) renormalize all non‐tumor cell fractions so they sum to 1
  #    by dividing each by (1 – tumor_purity)
  mutate(
    across(
      `Oligodendrocyte`:`Radial_Glia`,
      ~ . / (1 - tumor_purity)
    )
  )

# sanity check: the non‐tumor fractions should now sum to 1
data2 %>%
  mutate(check = rowSums(across(`Oligodendrocyte`:`Radial_Glia`))) %>%
  dplyr::select(sample_id, tumor_purity, check) %>%
  head(10)

df <- as.data.frame(data2) %>%
  select(sample_id, location, OPC_like_Cancer, AC_like_Cancer, NPC_like_Cancer, MES_like_Cancer) %>%
  mutate(sum = rowSums(across(OPC_like_Cancer:MES_like_Cancer), na.rm = TRUE)) %>%
  mutate(across(
    .cols = OPC_like_Cancer:MES_like_Cancer,
    .fns = ~ .x / sum,
    .names = "{.col}_prop"
  ))

# make plots ####
# Define the custom order for the locations on the x-axis
location_order <- c("Perinecrotic Zone", "Pseudopalisading Cells", "Cellular Tumor", "Infiltrating Tumor", "Leading Edge")

# Define the desired stacking order for the cell types (bottom to top)
cell_type_order <- c("MES-like", "AC-like", "OPC-like", "NPC-like")

# Process the data for plotting:
# 1. Calculate the average proportion for each cell type, grouped by location.
# 2. Reshape the data from a "wide" to a "long" format, which is required by ggplot.
# 3. Clean up the cell type names and set the factor order for both location and cell_type.
plot_data <- df %>%
  # Select only the location and the proportion columns
  select(location, ends_with("_prop")) %>%
  # Ensure we only include the locations specified in our order
  # This also filters out any other locations you may not want to plot
  filter(location %in% location_order) %>%
  # Group by the location to calculate averages
  group_by(location) %>%
  # Calculate the mean for each cell type within each location
  summarise(across(everything(), mean, na.rm = TRUE), .groups = 'drop') %>%
  # Pivot the data from wide to long format
  pivot_longer(
    cols = -location,
    names_to = "cell_type",
    values_to = "average_proportion"
  ) %>%
  # Perform some string cleaning on the cell type names for the legend
  mutate(
    cell_type = str_replace(cell_type, "_prop", ""),
    cell_type = str_replace(cell_type, "_like_Cancer", "-like"),
    # Convert the 'location' column to a factor to enforce the custom order
    location = factor(location, levels = location_order),
    # Convert 'cell_type' to a factor to enforce the stacking order
    cell_type = factor(cell_type, levels = cell_type_order)
  )


# --- Plotting ---

# Define the custom color palette based on the image provided
# AC=Red, MES=Green, OPC=Blue/Teal, NPC=Purple
color_palette <- c(
  "AC-like" = "#F8766D",
  "MES-like" = "#7CAE00",
  "OPC-like" = "#00BFC4",
  "NPC-like" = "#C77CFF"
)

# Create the 100% stacked area chart using ggplot2
ggplot(plot_data, aes(x = location, y = average_proportion, fill = cell_type, group = cell_type)) +
  # Use geom_area with position="fill", add black borders, and set transparency
  geom_area(position = "fill", color = "black", alpha = 0.8) +
  
  # Manually set the fill colors using our defined palette
  scale_fill_manual(values = color_palette) +
  
  # Format the y-axis labels to be percentages (e.g., 25%, 50%)
  scale_y_continuous(labels = scales::percent_format()) +
  
  # Set the plot title, axis labels, and legend title
  labs(
    title = "Tumor Composition by Anatomic Location",
    x = "Tumor Location",
    y = "Cancer Cell Composition",
    fill = "Cell Type"
  ) +
  
  # Apply a clean, minimal theme
  theme_bw() +
  
  # Customize theme elements for better readability
  theme(
    # Rotate the x-axis labels to prevent them from overlapping
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 11),
    axis.title = element_text(size = 12),
    # Center the plot title
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    # Customize the legend
    legend.position = "right",
    legend.title = NULL,
    # Remove panel grid lines
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

