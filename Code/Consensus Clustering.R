#libraries####
library(data.table)
library(dplyr)
library(ggplot2)
library(stringr)
library(compositions)
library(M3C)
library(sva)
library(survival)
library(survminer)
library(caret)
library(randomForest)
library(pheatmap)
library(RColorBrewer)






#Import & clean data####
subtype_cols <- c("OPC_like_Cancer", "AC_like_Cancer", "NPC_like_Cancer", "MES_like_Cancer")
setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME")
data <- as.data.frame(fread("GBM Clinical-TME Cell States.csv")) %>%
  dplyr::mutate(subtype = subtype_cols[
    max.col( as.matrix( across(all_of(subtype_cols)) ),
             ties.method = "first")])

df <- data %>% filter(recurrence == "primary") %>%
  filter(platform == "rnaseq")

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


# ComBat batch correction####
df_clr_rnaseq <- df_clr %>%
  filter(platform == "rnaseq")

start_col <- "Oligodendrocyte_clr"
end_col <- "MES_like_Cancer_clr"
cell_type_cols <- names(df_clr_rnaseq)[which(names(df_clr_rnaseq) == start_col):which(names(df_clr_rnaseq) == end_col)]

# Prepare the data matrix for ComBat.
# ComBat requires a matrix where features (cell types) are rows and samples are columns.
# We select the relevant columns and transpose the resulting matrix.
expression_matrix <- t(df_clr_rnaseq[, cell_type_cols])

# Define the batch variable. This is the factor that you want to correct for.
# In this case, it's the 'study' column.
batch_info <- df_clr_rnaseq$study


# --- 3. Running ComBat ---
# ComBat can optionally preserve biological variation while removing batch effects.
# We create a model matrix for any covariates we want to protect from correction.
# Using model.matrix(~1, ...) creates a model with just an intercept,
# which is appropriate when you don't have other biological variables to preserve.
modcombat <- model.matrix(~1, data = df_clr_rnaseq)

# Execute the ComBat function.
# 'par.prior=TRUE' specifies that the parametric empirical Bayes method should be used.
# This is generally recommended for datasets with more than a handful of samples per batch.
combat_corrected_matrix <- ComBat(dat = expression_matrix, 
                                  batch = batch_info, 
                                  mod = modcombat, 
                                  par.prior = TRUE)


# --- 4. Finalizing the Corrected Dataframe ---
# Transpose the corrected matrix back to the original orientation (samples as rows).
combat_corrected_df <- as.data.frame(t(combat_corrected_matrix))

# Extract the metadata (all columns except the ones that were corrected).
metadata <- df_clr_rnaseq %>%
  dplyr::select(-one_of(cell_type_cols))

# Combine the metadata with the newly corrected cell type data.
df_clr_corrected <- cbind(metadata, combat_corrected_df)

# --- 5. Verification ---
cat("\nDimensions of original filtered data:", dim(df_clr_rnaseq), "\n")
cat("Dimensions of corrected data:", dim(df_clr_corrected), "\n")

write.csv(df_clr_corrected, "GBM_TME_CellStates_CLR_ComBat_Corrected.csv", row.names = FALSE)


#consensus Clustering####
tempus_df <- df_clr %>%
  filter(study == "TEMPUS")

# Select the columns for clustering. We are selecting all columns from
# 'Oligodendrocyte_clr' to 'Radial_Glia_clr'.
cols_for_clustering <- colnames(tempus_df)
start_col <- which(cols_for_clustering == "Oligodendrocyte_clr")
end_col <- which(cols_for_clustering == "MES_like_Cancer_clr")

# Create a matrix of the data to be clustered.
# ConsensusClusterPlus requires a matrix where rows are features (cell types)
# and columns are samples. Therefore, we select the columns and then transpose the matrix.
data_matrix <- t(tempus_df[, start_col:end_col])

# Check for and handle missing values (NAs) if any.
# For this example, we're assuming no NAs. If there are, you might need to
# remove them or use an imputation method.
if(any(is.na(data_matrix))){
  warning("NA values detected in the clustering data. This may cause errors. Consider imputation or removal.")
}

# --- 4. Run Consensus Clustering ---
#
# This function will repeatedly cluster the data and assess the stability
# of the clusters.
#
# Parameters:
# - d: The data matrix (features in rows, samples in columns)
# - maxK: The maximum number of clusters (k) to evaluate. We'll check from k=2 to k=6.
# - reps: Number of bootstrap repetitions (1000 is a good standard for robust results).
# - pItem: The proportion of samples to include in each bootstrap iteration.
# - pFeature: The proportion of features (cell types) to include.
# - clusterAlg: The clustering algorithm. 'pam' (Partitioning Around Medoids) is robust.
# - distance: The distance metric. 'pearson' is often used for gene expression-like data.
# - title: The directory where the output plots will be saved.
#
results_dir <- "consensus_clustering_output"
if (!dir.exists(results_dir)) {
  dir.create(results_dir)
}

# Execute the clustering
# Note: This step can be computationally intensive and may take some time.
consensus_results <- ConsensusClusterPlus(
  d = data_matrix,
  maxK = 10,
  reps = 1000,
  pItem = 0.8,
  pFeature = 1,
  clusterAlg = "pam",
  distance = "pearson",
  title = results_dir,
  seed = 123456, # for reproducibility
  plot = "pdf"
)


# --- 5. Interpreting the Results and Assigning Ecotypes ---
#
# After running the code, check the 'consensus_clustering_output' directory.
#
# Key plots to examine:
# 1. Delta Area Plot (delta_area.pdf): This plot shows the relative change in the
#    area under the CDF curve. The "elbow" of this plot indicates the k
#    value where the stability of the clusters is maximized. This is often the
#    best choice for the number of clusters.
# 2. Consensus CDF Plot (consensus_cdf.pdf): Shows the cumulative distribution
#    function for each k. The flattest middle segment indicates the most stable k.
#
# Let's assume after reviewing the plots, you decide that k=3 is the optimal number
# of clusters.

# Set your chosen k
optimal_k <- 2

# Extract the cluster assignments for the chosen k
ecotype_assignments <- consensus_results[[optimal_k]][["consensusClass"]]

# Create a new dataframe with patient_id/sample_id and their assigned ecotype
ecotype_df <- data.frame(
  sample_id = tempus_df$sample_id,
  ecotype = paste0("E", ecotype_assignments) # Naming them E1, E2, E3...
)

# Merge the ecotype assignments back into your original filtered dataframe
tempus_df_with_ecotypes <- tempus_df %>%
  left_join(ecotype_df, by = "sample_id")



#Monte Carlo Concensus Clustering####
tempus_df <- df_clr %>%
  filter(study == "GSAM")

cols_for_clustering <- colnames(tempus_df)
start_col <- which(cols_for_clustering == "Oligodendrocyte_clr")
end_col <- which(cols_for_clustering == "MES_like_Cancer_clr")

# Create a matrix where rows are features (cell types) and columns are samples.
data_matrix <- t(tempus_df[, start_col:end_col])
colnames(data_matrix) <- tempus_df$sample_id

data_matrix <- data.frame(data_matrix)
sanitized_names <- colnames(data_matrix)
tempus_df$sample_id <- sanitized_names

# --- 3. Run Monte Carlo Reference-based Consensus Clustering (M3C) ---
# M3C evaluates cluster stability against a null distribution generated from
# reference datasets (created by permuting the input data).
# The key output is the Relative Cluster Stability Index (RCSI), where the
# peak value indicates the optimal number of clusters (k).

# Note: This step can be computationally intensive.
# M3C automatically tests a range of k values.
res <- M3C(mydata = data_matrix, cores = 4, iters = 25, maxK = 10, des = NULL,
           pItem = 0.8, seed = 123)
           
# --- 4. Interpreting the M3C Results ---
# The primary output to check is the RCSI plot. M3C provides a built-in
# plotting function that shows the RCSI score, the consensus score, and the
# reference (null distribution) score for each k. The optimal k is where the
# RCSI score peaks.

# The optimal k is the one with the highest RCSI score.
optimal_k <- res$scores$K[which.max(res$scores$RCSI)]

# --- 5. Assign Ecotypes Based on Optimal K ---
# Extract the cluster assignments for the chosen optimal k.
# The assignments are stored within the results object.
ecotype_assignments <- res$assignments

# Create a new dataframe with sample_id and their assigned ecotype
ecotype_df <- data.frame(
  sample_id = tempus_df$sample_id,
  ecotype = paste0("E", ecotype_assignments) # Naming them E1, E2, E3...
)

# Merge the ecotype assignments back into your original dataframe
tempus_df_with_ecotypes <- tempus_df %>%
  left_join(ecotype_df, by = "sample_id")

#heatmap
tme_data$ecotype <- as.factor(tme_data$ecotype)

matrix_data <- tme_data[, -which(names(tme_data) == "ecotype")]
rownames(matrix_data) <- paste0("Sample_", 1:nrow(matrix_data))

annotation_data <- data.frame(
  Ecotype = tme_data$ecotype
)
rownames(annotation_data) <- rownames(matrix_data)

ecotype_levels <- levels(tme_data$ecotype)
ecotype_colors <- brewer.pal(length(ecotype_levels), "Set1")
names(ecotype_colors) <- ecotype_levels

annotation_colors <- list(
  Ecotype = ecotype_colors
)

pheatmap(
  t(matrix_data)[, order(tme_data$ecotype)], # Transpose and sort matrix by ecotype
  cluster_rows = TRUE,                   # Cluster the rows (variables)
  cluster_cols = FALSE,                  # Do not cluster columns to keep ecotype grouping
  show_colnames = FALSE,                 # Hide sample names for a cleaner look
  annotation_col = annotation_data,      # Use the ecotype annotation for columns
  annotation_colors = annotation_colors, # Specify the colors for the annotation
  scale = "none",                      # Scale values by column (sample Z-score)
  main = "Transposed Heatmap of TME Variables Grouped by Ecotype" # Title for the plot
)

# Ecotype classifier####
# Select the predictor variables (all TME components ending in '_clr')
# and the outcome variable ('ecotype').
tme_data <- tempus_df_with_ecotypes %>%
  select(ecotype, ends_with("_clr"))

# Convert the outcome variable 'ecotype' to a factor, which is required
# for classification models in R.
tme_data$ecotype <- as.factor(tme_data$ecotype)

# Check for missing values in the selected data.
# Random Forest cannot handle NAs, so we will remove rows with any missing values.
tme_data <- na.omit(tme_data)

# --- 4. Data Splitting ---
#
# We will split the data into a training set (80%) and a testing set (20%).
# The model will be built on the training set and evaluated on the unseen testing set.
#
set.seed(123) # for reproducibility
train_index <- createDataPartition(tme_data$ecotype, p = 0.8, list = FALSE)
train_data <- tme_data[train_index, ]
test_data  <- tme_data[-train_index, ]

# --- 5. Model Training ---
#
# We will train the Random Forest classifier using the 'caret' package.
# 'trainControl' specifies 10-fold cross-validation to ensure the model is robust.
#
set.seed(123) # for reproducibility
fit_control <- trainControl(method = "cv", number = 10)

# Train the model
# This may take a moment to run.
rf_model <- train(
  ecotype ~ .,          # Formula: predict ecotype using all other variables
  data = train_data,
  method = "rf",        # Specify the random forest algorithm
  trControl = fit_control,
  verbose = FALSE       # Suppress verbose output
)

# Print the model summary
print("Random Forest Model Summary:")
print(rf_model)


# --- 6. Model Evaluation ---
#
# Now, we use the trained model to make predictions on our test set.
#
predictions <- predict(rf_model, newdata = test_data)

# Create a confusion matrix to see the performance in detail.
# The rows show the actual classes, and the columns show the predicted classes.
confusion_mat <- confusionMatrix(predictions, test_data$ecotype)

print("Confusion Matrix and Performance Statistics:")
print(confusion_mat)


# --- 7. Feature Importance ---
#
# A major benefit of Random Forest is its ability to rank variables by importance.
# We will extract and plot the top 20 most important TME components.
#
var_imp <- varImp(rf_model, scale = TRUE)

# Create a plot of variable importance
imp_plot <- plot(var_imp, top = 20, main = "Top 20 Most Important TME Components")

print("Variable Importance Plot:")
print(imp_plot)


# Apply classifier to the rest of the data####
other_cohort_data <- df_clr %>%
  filter(study != "TEMPUS")

# Select the same predictor columns the model was trained on.
# We keep 'sample_id' to merge the results back later.
other_cohort_predictors <- other_cohort_data %>%
  select(sample_id, ends_with("_clr"))

# The model cannot predict on data with missing values.
# We'll remove rows with NAs. For a more advanced analysis,
# you could consider imputation methods instead.
other_cohort_predictors_clean <- na.omit(other_cohort_predictors)

# Use the trained model to predict ecotypes for the clean data
other_cohort_predictions <- predict(rf_model, newdata = other_cohort_predictors_clean)

# Combine the predictions with the sample IDs
other_cohort_classified <- data.frame(
  sample_id = other_cohort_predictors_clean$sample_id,
  ecotype_predicted = other_cohort_predictions
)

# Merge these new predictions back into the full 'df_clr' dataframe
df_clr_with_all_ecotypes <- df_clr %>%
  left_join(
    select(tempus_df_with_ecotypes, sample_id, ecotype), by = "sample_id"
  ) %>%
  left_join(other_cohort_classified, by = "sample_id") %>%
  mutate(
    ecotype_final = ifelse(!is.na(ecotype), as.character(ecotype), as.character(ecotype_predicted))
  ) %>%
  select(-ecotype, -ecotype_predicted) # Clean up intermediate columns

# test ecotypes and survival: KM curve####
studies_to_plot <- c("TCGA")

for (study_name in studies_to_plot) {
  
  # Filter the data for the current study
  survival_data <- tempus_df_with_ecotypes %>%
    filter(subtype == "MES_like_Cancer")
  
  # Ensure the necessary columns exist and are not all NA
  if (nrow(survival_data) > 0 && "OS_years" %in% names(survival_data) && "OS_event" %in% names(survival_data)) {
    
    # Remove rows where survival information or ecotype is missing
    survival_data <- survival_data %>%
      filter(!is.na(OS_years) & !is.na(OS_event) & !is.na(ecotype))
    
    # Ensure there is still data to plot
    if (nrow(survival_data) > 0) {
      
      # Create the survival object
      surv_obj <- Surv(time = survival_data$OS_years, event = survival_data$OS_event)
      
      # Create the survival fit object
      fit <- survfit(surv_obj ~ ecotype, data = survival_data)
      
      # Generate the Kaplan-Meier plot
      km_plot <- ggsurvplot(
        fit,
        data = survival_data,
        title = paste("Kaplan-Meier Curve"),
        pval = TRUE,
        conf.int = TRUE,
        risk.table = TRUE,
        xlab = "Time in Years",
        ylab = "Overall Survival Probability",
        legend.title = "Ecotype",
        palette = "jco" # Journal of Clinical Oncology color palette
      )
      
      # Print the plot
      print(km_plot)
      
    } else {
      print(paste("No complete survival data available for study:", study_name))
    }
  } else {
    print(paste("Skipping study:", study_name, "- not found or missing required survival columns."))
  }
}
# test ecotypes and survival: Cox####
cox_data <- tempus_df_with_ecotypes

cox_data <- cox_data %>%
  filter(!is.na(OS_years) & !is.na(OS_event) & !is.na(ecotype_final) & !is.na(study)) %>%
  mutate(
    # Ensure ecotype and study are factors for the model
    ecotype_final = as.factor(ecotype_final),
    study = as.factor(study)
  )

# Set a reference level for the ecotype if desired (e.g., E1)
# This makes the hazard ratios interpretable relative to this reference group.
if("E1" %in% levels(cox_data$ecotype_final)) {
  cox_data$ecotype_final <- relevel(cox_data$ecotype_final, ref = "E1")
}


# Check if there is enough data to build the model
if (nrow(cox_data) > 0) {
  
  # Build the stratified Cox Proportional Hazards model
  # The formula models survival based on ecotype, stratifying by study
  cox_model <- coxph(
    Surv(OS_years, OS_event) ~ ecotype,
    data = cox_data
  )
  
  # Print the summary of the Cox model
  # This provides hazard ratios (exp(coef)), confidence intervals, and p-values.
  print("Summary of Stratified Cox Proportional Hazards Model:")
  print(summary(cox_model))
  
  # Generate a forest plot to visualize the hazard ratios for each ecotype
  # relative to the reference level.
  print("Forest Plot of Hazard Ratios:")
  forest_plot <- ggforest(cox_model, data = cox_data)
  print(forest_plot)
  
} else {
  print("Not enough complete data across the specified studies to build a Cox model.")
}
