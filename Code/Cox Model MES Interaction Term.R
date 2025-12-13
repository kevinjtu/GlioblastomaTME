#libraries####
library(data.table)
library(dplyr)
library(ggplot2)
library(stringr)
library(compositions)
library(survival)
library(survminer)
library(ggeffects)
library(patchwork)

#import & clean data - CELL TYPES####
setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME")
data <- as.data.frame(fread("GBM Clinical-TME Cell Types.csv"))

data2 <- data %>%
  mutate(tumor_purity = Malignant) %>%
  dplyr::select(-c(Malignant)) %>%
  
  # 2) renormalize all non‐tumor cell fractions so they sum to 1
  #    by dividing each by (1 – tumor_purity)
  mutate(
    across(
      `Oligodendrocyte`:`radial_glial_cell`,
      ~ . / (1 - tumor_purity)
    )
  )

# sanity check: the non‐tumor fractions should now sum to 1
data2 %>%
  mutate(check = rowSums(across(`Oligodendrocyte`:`radial_glial_cell`))) %>%
  dplyr::select(sample_id, tumor_purity, check) %>%
  head(10)

data2 <- as.data.frame(data2)

df <- data2 %>%
  filter(recurrence != "recurrent") %>%
  filter(is.na(location))

cell_cols <- names(df)[ which(names(df) == "Oligodendrocyte"):
                          which(names(df) == "radial_glial_cell") ]

#Subtypes ####
subtype_cols <- c("OPC_like_Cancer", "AC_like_Cancer", "NPC_like_Cancer", "MES_like_Cancer")
setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME")
subtype <- as.data.frame(fread("GBM Clinical-TME Cell States.csv")) %>%
  dplyr::mutate(subtype = subtype_cols[
    max.col( as.matrix( across(all_of(subtype_cols)) ),
             ties.method = "first")]) %>%
  dplyr::select(patient_id, sample_id, subtype)

df <- df %>%
  left_join(subtype, by = c("patient_id", "sample_id"))

#CLR Transformation####
# 1) extract the raw proportions matrix
df_cells <- df %>% dplyr::select(all_of(cell_cols))

# 2) add a tiny pseudocount to avoid zeros
#pcount <- min(df_cells[df_cells > 0], na.rm = TRUE) * 0.5
#df_cells[df_cells == 0] <- pcount

# 3) compute the CLR
clr_mat <- clr(as.matrix(df_cells))

# 4) rename columns to mark them as CLR‐transformed
#colnames(clr_mat) <- paste0(cell_cols, "_clr")

# 5) rebuild your analysis data frame, replacing raw with CLR
df_clr <- df %>%
  dplyr::select(-all_of(cell_cols)) %>%
  bind_cols(as.data.frame(clr_mat))

#Cell type matching & Cox Model####
#quartilize Oligodendrocyte:radial_glial_cell
columns_to_quartilize <- c("Oligodendrocyte", "Pericyte", "Neuron", "Astrocyte",
                           "Endothelial", "Myeloid", "Lymphoid", "OPC", 
                           "radial_glial_cell")

# Loop through each column and create a new column with quartile bins
for (col in columns_to_quartilize) {
  quartile_col <- paste0(col, "_quartile")
  
  df_clr <- df_clr %>%
    group_by(study) %>%
    mutate(
      !!quartile_col := as.integer(cut(
        !!sym(col),
        breaks = quantile(!!sym(col), probs = seq(0, 1, 0.25), na.rm = TRUE),
        labels = 1:4,
        include.lowest = TRUE
      ))
    ) %>%
    ungroup()
}

#cox models
cell_types <- c(
  "Oligodendrocyte", "Pericyte", "Neuron", "Astrocyte",
  "Endothelial", "Myeloid", "Lymphoid", "OPC", "radial_glial_cell"
)

# Step 4: Initialize lists to store models and plots
cox_models <- list()
interaction_plots <- list()

# Step 5: Loop through each cell type to run analysis and generate plots
for (cell in cell_types) {
  # Define the specific quartile column for the current cell type
  quartile_col <- paste0(cell, "_quartile")
  
  # Construct the formula for the Cox model
  # The formula is: survival ~ cell_quartile * subtype + strata(study)
  df_clr$subtype <- relevel(factor(df_clr$subtype), ref = "MES_like_Cancer")
  formula_str <- paste0("Surv(OS_years, OS_event) ~ ", quartile_col, " * subtype + strata(study)")
  cox_formula <- as.formula(formula_str)
  
  # Run the Cox proportional hazards model
  # The model automatically handles missing data by omitting rows with NA
  # in any of the model's variables (OS_years, OS_event, quartile_col, subtype, study).
  fit <- coxph(cox_formula, data = df_clr)
  
  # Store the model summary
  cox_models[[cell]] <- summary(fit)
  
  # Generate predictions for the interaction term using ggeffects
  # We predict the survival probability for each subtype across the cell quartiles
  # Note: Use `ggpredict` for survival probabilities
  pred_effects <- ggpredict(fit, terms = c(quartile_col, "subtype"))
  
  # Create the plot for the current cell type
  p <- plot(pred_effects) +
    labs(
      title = paste(cell),
      x = paste("Quartile"),
      y = "Hazard Ratio",
      colour = "Subtype"
    ) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 12)) +
    theme(
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
      )
  
  # Store the plot
  interaction_plots[[cell]] <- p
}

# Step 6: Print the summaries of the Cox models
# This will show the coefficients, hazard ratios, and p-values for each model
for (cell in cell_types) {
  cat("==============================================================\n")
  cat("Cox Model Summary for:", cell, "\n")
  cat("==============================================================\n")
  print(cox_models[[cell]])
  cat("\n\n")
}

# Define the desired order
selected_cells <- c("Lymphoid", "Myeloid", "Neuron", "Pericyte")

# Subset and reorder interaction_plots
interaction_plots <- interaction_plots[selected_cells]

# Step 7: Display the combined plot of all interaction effects
# The `patchwork` package arranges the ggplot objects into a grid
combined_plot <- wrap_plots(interaction_plots, ncol = 4) +
  plot_annotation(
    title = " Interaction Between Cell Type  and Subtype",
    theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
  )

# Show the final combined plot
print(combined_plot)

#import & clean data - CELL STATES, select for significant####
subtype_cols <- c("OPC_like_Cancer", "AC_like_Cancer", "NPC_like_Cancer", "MES_like_Cancer")
setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME")
data <- as.data.frame(fread("GBM Clinical-TME Cell States.csv")) %>%
  mutate(subtype = subtype_cols[
    max.col( as.matrix( across(all_of(subtype_cols)) ),
             ties.method = "first")])

tumor_cells <- c("OPC_like_Cancer","AC_like_Cancer",
                 "NPC_like_Cancer","MES_like_Cancer")

data2 <- data %>%
  rowwise() %>%
  mutate(tumor_purity = sum(c_across(all_of(tumor_cells)))) %>%
  ungroup() %>%
  mutate(across(
    # only non-tumor up through Radial_Glia
    Oligodendrocyte:Radial_Glia, 
    ~ . / (1 - tumor_purity)
  )) %>%
  dplyr::select(-all_of(tumor_cells))

df <- data2 %>%
  filter(recurrence == "primary") %>%
  mutate(
    study_f = factor(study),
    STING_s = as.numeric(REACTOME_STING)
  )

cell_cols <- names(df)[
  which(names(df)=="Oligodendrocyte") :
    which(names(df)=="Radial_Glia")
]

#select lineages
lineages <- c(
  CD8_Effector_Memory_T_Cells  = "Lymphoid",
  Pro_inflammatory_Microglia   = "Myeloid",
  Radial_Glia                  = "Normal Tissue",
  Astrocyte                    = "Normal Tissue",
  OPC                          = "Normal Tissue",
  Neuron                       = "Normal Tissue",
  Oligodendrocyte              = "Normal Tissue",
  Arterial_Vessel              = "Endothelial",
  Capilary_Vessel              = "Endothelial",
  Tiplike_Vessel               = "Endothelial",
  Pericyte                     = "Endothelial",
  Hypoxia_associated_Monocytes = "Myeloid",
  Naive_Monocyte               = "Myeloid",
  BDM_Hypoxia_and_MES_like     = "Myeloid",
  Anti_inflammatory_Monocyte   = "Myeloid",
  Anti_inflammatory_BDM        = "Myeloid",
  Microglia_Aging_Signature    = "Myeloid",
  BDM_IFN_Signature            = "Myeloid",
  Proliferative_Microglia      = "Myeloid",
  BDM_MHC_Expression           = "Myeloid",
  Conventional_DC_1            = "Myeloid",
  Conventional_DC_2            = "Myeloid",
  Dendritic_like_Cell          = "Myeloid",
  Plasmacytoid_DC              = "Myeloid",
  Mast_Cell                    = "Myeloid",
  T_Cells_Stress_Signature     = "Lymphoid",
  Proliferating_T_Cells        = "Lymphoid",
  NK_like_CD8_T_Cells          = "Lymphoid",
  Cytotoxic_CD8_T_Cell         = "Lymphoid",
  Plasma_Cell                  = "Lymphoid",
  B_Cell                       = "Lymphoid",
  Resting_CD4_T_Cells          = "Lymphoid",
  CD4_T_Cell_IFN_Signature     = "Lymphoid",
  Regulatory_T_Cells           = "Lymphoid",
  Natural_Killer_Cells         = "Lymphoid"
)

#CLR Transformation####
# 1) extract the raw proportions matrix
df_cells <- df %>% dplyr::select(all_of(cell_cols))

# 2) add a tiny pseudocount to avoid zeros
#pcount <- min(df_cells[df_cells > 0], na.rm = TRUE) * 0.5
#df_cells[df_cells == 0] <- pcount

# 3) compute the CLR
clr_mat <- clr(as.matrix(df_cells))

# 4) rename columns to mark them as CLR‐transformed
#colnames(clr_mat) <- paste0(cell_cols, "_clr")

# 5) rebuild your analysis data frame, replacing raw with CLR
df_clr <- df %>%
  dplyr::select(-all_of(cell_cols)) %>%
  bind_cols(as.data.frame(clr_mat))

# Cell State Cox Model, Lymphoid####
#quartilize Oligodendrocyte:radial_glial_cell
# First, extract only lymphoid cell names from the lineage map
lymphoid_cells <- names(lineages[lineages %in% c("Normal Tissue", "Lymphoid", "Myeloid", "Endothelial")])

# Then filter those that are actually present as columns in df
columns_to_quartilize <- intersect(lymphoid_cells, names(df_clr))

# Loop through each column and create a new column with quartile bins
for (col in columns_to_quartilize) {
  quartile_col <- paste0(col, "_quartile")
  
  df_clr <- df_clr %>%
    group_by(study) %>%
    mutate(
      !!quartile_col := as.integer(cut(
        !!sym(col),
        breaks = quantile(!!sym(col), probs = seq(0, 1, 0.25), na.rm = TRUE),
        labels = 1:4,
        include.lowest = TRUE
      ))
    ) %>%
    ungroup()
}

# Step 4: Initialize lists to store models and plots
cox_models <- list()
interaction_plots <- list()

# Step 5: Loop through each cell type to run analysis and generate plots
for (cell in cell_types) {
  # Define the specific quartile column for the current cell type
  quartile_col <- paste0(cell, "_quartile")
  
  # Construct the formula for the Cox model
  # The formula is: survival ~ cell_quartile * subtype + strata(study)
  df_clr$subtype <- relevel(factor(df_clr$subtype), ref = "MES_like_Cancer")
  formula_str <- paste0("Surv(OS_years, OS_event) ~ ", quartile_col, " * subtype + strata(study)")
  cox_formula <- as.formula(formula_str)
  
  # Run the Cox proportional hazards model
  # The model automatically handles missing data by omitting rows with NA
  # in any of the model's variables (OS_years, OS_event, quartile_col, subtype, study).
  fit <- coxph(cox_formula, data = df_clr)
  
  # Store the model summary
  cox_models[[cell]] <- summary(fit)
  
  # Generate predictions for the interaction term using ggeffects
  # We predict the survival probability for each subtype across the cell quartiles
  # Note: Use `ggpredict` for survival probabilities
  pred_effects <- ggpredict(fit, terms = c(quartile_col, "subtype"))
  
  # Create the plot for the current cell type
  p <- plot(pred_effects) +
    labs(
      title = paste(cell),
      x = paste("Quartile"),
      y = "Hazard Ratio",
      colour = "Subtype"
    ) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 12)) +
    theme(
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  # Store the plot
  interaction_plots[[cell]] <- p
}

model_summary <- summary(cox_model)
model_coefs <- coef(cox_model)
vcov_matrix <- vcov(cox_model)

# --- Step 3: Define the combinations for which to predict the hazard ratio ---
# We create a dataframe representing all combinations of the interacting variables.
# Note: The reference level for TME_signature is not explicitly named in your output,
# so we'll call it 'Immunosuppressed' for clarity.
plot_data <- expand.grid(
  subtype = levels(cox_data$subtype),
  TME_signature = levels(cox_data$TME_signature)
)

# Filter out the NA interaction from the model (NPC_like_Cancer:Immunoactive)
# This combination doesn't have enough data to produce a result in your model.
plot_data <- plot_data %>%
  filter(!(subtype == "NPC_like_Cancer" & TME_signature == "Immunoactive"))

# --- Step 4: Calculate Hazard Ratios and Confidence Intervals ---

# The 'predict' function in R with type="terms" can help us get the linear predictors
# for each combination. We use a bit of a trick by creating a new data frame
# with our desired combinations and predicting on it.
# We also need to get the standard errors for these predictions to build CIs.

# Using predict() simplifies getting the correct linear combination of coefficients.
# 'se.fit=TRUE' gives us the standard errors for these combinations directly.
predictions <- predict(cox_model, newdata = plot_data, type = "terms", se.fit = TRUE)

# The 'predict' function with type="terms" returns a matrix for the value of each term.
# We need to sum the main effects and interaction effects for our plot.
term_names <- c(
    paste0("subtype", plot_data$subtype),
    paste0("TME_signature", plot_data$TME_signature),
    paste0("subtype", plot_data$subtype, ":TME_signature", plot_data$TME_signature)
)
# Match the term names from the model's coefficients
model_term_names <- names(coef(cox_model))

# Calculate the combined log-hazard for each row in plot_data
log_hazard <- numeric(nrow(plot_data))
se_log_hazard <- numeric(nrow(plot_data))

for (i in 1:nrow(plot_data)) {
    # Identify which terms from the model apply to this row
    is_ref_subtype <- plot_data$subtype[i] == "AC_like_Cancer"
    is_ref_tme <- plot_data$TME_signature[i] != "Immunoactive" # Assuming the other level is the ref

    terms_to_sum <- c()
    if (!is_ref_subtype) {
        terms_to_sum <- c(terms_to_sum, paste0("subtype", plot_data$subtype[i]))
    }
    if (!is_ref_tme) {
        terms_to_sum <- c(terms_to_sum, "TME_signatureImmunoactive")
    }
    if (!is_ref_subtype && !is_ref_tme) {
        interaction_term <- paste0("subtype", plot_data$subtype[i], ":TME_signatureImmunoactive")
        # Check if the interaction term exists in the model
        if(interaction_term %in% names(model_coefs)) {
          terms_to_sum <- c(terms_to_sum, interaction_term)
        }
    }

    # For the reference group (AC_like_Cancer & Immunosuppressed), log-hazard is 0
    if (length(terms_to_sum) == 0) {
        log_hazard[i] <- 0
        se_log_hazard[i] <- 0 # No variance for the reference
    } else {
        # Calculate log-hazard by summing relevant coefficients
        log_hazard[i] <- sum(model_coefs[terms_to_sum], na.rm = TRUE)

        # Calculate the variance of the summed coefficients using the vcov matrix
        # Var(A+B) = Var(A) + Var(B) + 2*Cov(A,B)
        sub_vcov <- vcov_matrix[terms_to_sum, terms_to_sum]
        se_log_hazard[i] <- sqrt(sum(sub_vcov))
    }
}


# Combine all calculated data
plot_data$log_hr <- log_hazard
plot_data$log_hr_se <- se_log_hazard

# Calculate Hazard Ratio and 95% Confidence Intervals
plot_data <- plot_data %>%
  mutate(
    HR = exp(log_hr),
    lower_ci = exp(log_hr - 1.96 * log_hr_se),
    upper_ci = exp(log_hr + 1.96 * log_hr_se)
  )

# --- Step 5: Generate the Plot with ggplot2 ---
interaction_plot <- ggplot(
    data = plot_data,
    aes(
        x = subtype,
        y = HR,
        group = TME_signature,
        color = TME_signature
    )
) +
  # Draw error bars for the confidence intervals
  geom_errorbar(
    aes(ymin = lower_ci, ymax = upper_ci),
    width = 0.2,
    position = position_dodge(0.3),
    linewidth = 0.8
  ) +
  # Draw lines connecting the points for each TME signature
  geom_line(
    position = position_dodge(0.3),
    linewidth = 1
  ) +
  # Draw points for the Hazard Ratios
  geom_point(
    position = position_dodge(0.3),
    size = 3.5,
    stroke = 1.5,
    aes(shape = TME_signature)
  ) +
  # Add a horizontal line at HR=1 for reference (no effect)
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    color = "grey40"
  ) +
  # Use a log scale for the y-axis, which is standard for Hazard Ratios
  scale_y_log10(
    breaks = c(0.25, 0.5, 1.0, 2.0, 4.0),
    labels = scales::label_number(accuracy = 0.1)
  ) +
  # Customize colors and labels for clarity
  scale_color_manual(
      values = c("Immunoactive" = "#0072B2", "Immunosuppressed" = "#D55E00"),
      name = "TME Signature"
  ) +
  scale_shape_manual(
      values = c("Immunoactive" = 16, "Immunosuppressed" = 17), # circle and triangle
      name = "TME Signature"
  ) +
  labs(
    title = "Interaction of Subtype and TME Signature on Survival",
    subtitle = "Hazard Ratios from Cox Proportional Hazards Model",
    x = "Cancer Subtype",
    y = "Hazard Ratio (log scale)",
    color = "TME Signature",
    shape = "TME Signature"
  ) +
  # Apply a clean theme
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
    axis.title = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    legend.position = "top",
    legend.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

# Display the plot
print(interaction_plot)
