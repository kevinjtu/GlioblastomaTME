#libraries####
library(data.table)
library(dplyr)
library(ggplot2)
library(forcats)
library(compositions)
library(tidyr)
library(ggpubr)

#import & clean data####
subtype_cols <- c("OPC_like_Cancer", "AC_like_Cancer", "NPC_like_Cancer", "MES_like_Cancer")
setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME")
data <- as.data.frame(fread("GBM Clinical-TME Cell States.csv")) %>%
  mutate(subtype = subtype_cols[
    max.col( as.matrix( across(all_of(subtype_cols)) ),
             ties.method = "first")])

df <- data #%>% filter(recurrence == "primary")

cell_cols <- names(df)[ which(names(df) == "Oligodendrocyte"):
                          which(names(df) == "MES_like_Cancer") ]

#CLR Transformation####
# 1) extract the raw proportions matrix
df_cells <- df %>% dplyr::select(all_of(cell_cols))

# 2) add a tiny pseudocount to avoid zeros
#pcount <- min(df_cells[df_cells > 0], na.rm = TRUE) * 0.5
#df_cells[df_cells == 0] <- pcount

# 3) compute the CLR
clr_mat <- clr(as.matrix(df_cells))

# 4) rename columns to mark them as CLR‐transformed
colnames(clr_mat) <- paste0(cell_cols, "_clr")

# 5) rebuild your analysis data frame, replacing raw with CLR
df_clr <- df %>%
  dplyr::select(-all_of(cell_cols)) %>%
  bind_cols(as.data.frame(clr_mat))

# 6) update the list of columns you’ll loop over
cell_cols_clr <- paste0(cell_cols, "_clr")


#PCA Analysis####
pca_col_names <- names(df_clr)[ which(names(df_clr) == "Oligodendrocyte_clr"):
                          which(names(df_clr) == "MES_like_Cancer_clr") ]


perform_pca_and_plot_by_study <- function(data, platform_name, pca_columns) {
  
  # Ensure the data for the platform is not empty
  if (nrow(data) < 2) {
    cat(paste("Skipping", platform_name, "due to insufficient data.\n"))
    return(NULL)
  }
  
  # --- PCA Calculation ---
  # Select only the numeric columns for PCA.
  # We scale the data so that variables with larger variance don't dominate the analysis.
  pca_result <- prcomp(data[, pca_columns], center = TRUE, scale. = TRUE)
  
  # Extract the variance explained by the first two principal components
  pca_summary <- summary(pca_result)
  pc1_variance <- round(pca_summary$importance[2, 1] * 100, 1)
  pc2_variance <- round(pca_summary$importance[2, 2] * 100, 1)
  
  # Combine PCA results (the new coordinates) with the original metadata (study info)
  pca_plot_df <- as.data.frame(pca_result$x)
  pca_plot_df$study <- data$study
  pca_plot_df$age <- data$age
  pca_plot_df$treatment <- data$treatment
  pca_plot_df$patient_id <- data$patient_id
  pca_plot_df$recurrence <- data$recurrence
  pca_plot_df$idh <- data$IDH_status
  pca_plot_df$mgmt <- data$MGMT_status
  pca_plot_df$KarnPerfScore <- data$KarnPerfScore
  pca_plot_df$sex <- data$sex
  pca_plot_df$resection <- data$resection
  
  
  
  
  
  
  
  # Get a list of unique studies within this dataset
  unique_studies <- unique(pca_plot_df$study)
  
  # Define the directory to save plots
  save_dir <- "/Users/kevintu/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Manuscript/Figures/PCA"
  # Create the directory if it doesn't exist
  if (!dir.exists(save_dir)) {
    dir.create(save_dir, recursive = TRUE)
  }
  
  # --- Plotting Loop ---
  # Loop through each unique study and create a dedicated plot
  for (current_study in unique_studies) {
    
    # Generate the plot using ggplot2. We use two geom_point layers. 
    # The second layer is drawn on top of the first.
    p <- ggplot(pca_plot_df, aes(x = PC1, y = PC2)) +
      # Layer 1: Plot all other studies (the background points) in black.
      geom_point(data = . %>% filter(study != current_study), color = "black", size = 0.5, alpha = 0.6) +
      # Layer 2: Plot the highlighted study in blue on top.
      geom_point(data = . %>% filter(study == current_study), color = "#0F9ED5", size = 1, alpha = 1.0) +
      
      # Add a simplified title
      labs(
        title = current_study
      ) +
      theme_bw() + # Use a black & white theme
      theme(
        legend.position = "none", # Hide the legend as per the request
        plot.title = element_text(hjust = 0.5, face = "bold", size = 10), # Center the title
        
        # Remove all axis text, titles, and tick marks
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        
        # Remove gridlines
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      )
    
    # Print the plot to the plots pane
    print(p)
    
    # --- Save the Plot ---
    # Construct a clean filename
    file_name <- paste0("PCA_", gsub(" ", "_", platform_name), "_Highlighting_", current_study, ".png")
    # Save the plot using ggsave
    ggsave(
      filename = file_name,
      plot = p,
      path = save_dir,
      width = 1.2,
      height = 1.6,
      units = "in",
      dpi = 300,
      bg = "transparent" # Set background to transparent for saving
    )
  }
}

# --- 4. Execute the Analysis ---

# **IMPORTANT**: Make sure your `df_clr` is loaded and ready before this step.
# Also, ensure `pca_col_names` correctly lists all columns from
# `Oligodendrocyte_clr` through `MES_like_Cancer_clr`.

# Split the main dataframe by the 'platform' column
df_rnaseq <- df_clr #%>% filter(platform == "rnaseq")
#df_microarray <- df_clr %>% filter(platform == "microarray")

# Generate plots for the RNA-Seq data
cat("--- Generating PCA plots for RNA-Seq data ---\n")
perform_pca_and_plot_by_study(df_rnaseq, "RNA-Seq", pca_col_names)

# Generate plots for the Microarray data
cat("\n--- Generating PCA plots for Microarray data ---\n")
perform_pca_and_plot_by_study(df_microarray, "Microarray", pca_col_names)
