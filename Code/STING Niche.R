#libraries####
library(data.table)
library(dplyr)
library(survival)
library(survminer)
library(ggplot2)
library(forcats)
library(compositions)
library(tidyr)
library(ggpubr)
library(purrr)
library(ppcor)
library(ggrepel)
library(ggeffects)

#Set  survival metric here: "OS", "PFS", or "DFS"
Surv_metric <- "OS"
time_column <- paste0(Surv_metric, "_years")
event_column <- paste0(Surv_metric, "_event")

#import & clean data####
subtype_cols <- c("OPC_like_Cancer", "AC_like_Cancer", "NPC_like_Cancer", "MES_like_Cancer")
setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME")
data <- as.data.frame(fread("GBM Clinical-TME Cell States.csv")) %>%
  mutate(subtype = subtype_cols[
    max.col( as.matrix( across(all_of(subtype_cols)) ),
             ties.method = "first")])

df <- data %>% filter(recurrence == "primary")

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

# STING Niche by location####
df_long <- df_clr %>%
  filter(study == "IVYGAP") %>%
  dplyr::select(sample_id, location, niche_all) %>%
  filter(!location %in% c(
    "Hyperplastic blood vessels in cellular tumor",
    "Microvascular Proliferation"
  )) %>%
  mutate(location = recode(location,
                           "Perinecrotic Zone" = "PZ",
                           "Pseudopalisading Cells" = "PC",
                           "Cellular Tumor" = "CT",
                           "Infiltrating Tumor" = "IT",
                           "Leading Edge" = "LE"
  ))

# Reorder the levels in desired order
df_long$location <- factor(df_long$location, levels = c("PZ", "PC", "CT", "IT", "LE"))

comparisons <- list(
  c("PZ", "PC"),
  c("PZ", "CT"),
  c("PZ", "IT"),
  c("PZ", "LE")
  )

# Create the plot
# Plotting code in the style of the attached image
ggplot(df_long, aes(x = location, y = niche_all)) +
  geom_jitter(shape = 16, color = "black", 
              position = position_jitter(width = 0.25, height = 0)) +
  stat_summary(fun = "mean", geom = "crossbar", 
               width = 0.5, fatten = 1.5, color = "red") +
  theme_classic() +
  labs(
    title = "STING Niche GSEA", # Title from the example image
    x = NULL,        # The image has no x-axis label
    y = "Enrichment Score"           # Set an empty string for the y-axis label
  ) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5), 
    axis.ticks = element_blank()
  ) +
    stat_compare_means(
    comparisons = comparisons, 
    method = "wilcox.test", 
    label = "p.signif"
  )
# STING Niche associtions with cell types####

# 1. Identify the columns for the cell types you want to analyze.
# This code dynamically finds all columns between 'Oligodendrocyte_clr' and 'MES_like_Cancer_clr'.
all_colnames <- names(df_clr)
start_col <- which(all_colnames == "Oligodendrocyte_clr")
end_col <- which(all_colnames == "MES_like_Cancer_clr")
cell_type_columns <- all_colnames[start_col:end_col]

# 2. Prepare the control variable.
# The 'study' column needs to be converted to a numeric factor for the correlation test.
df_clr$study_numeric <- as.numeric(as.factor(df_clr$study))

# 3. Initialize an empty list to store the results from each loop iteration.
results_list <- list()

# 4. Loop through each cell type column and perform the analysis.
for (cell_type in cell_type_columns) {
  
  # Print a message to track progress (optional)
  # message("Processing: ", cell_type)
  
  # Use tryCatch to handle potential errors during the correlation test,
  # for instance if a column has insufficient data.
  test_result <- tryCatch({
    
    # Extract the column data for the current cell type.
    # Using [[cell_type]] is a robust way to select a column by its name as a string.
    y_variable <- df_clr[[cell_type]]
    
    # Perform the partial spearman correlation test.
    # This tests the correlation between 'niche_all' and the current 'cell_type' column,
    # while controlling for the 'study' variable.
    ppcor::spcor.test(x = df_clr$niche_all, y = y_variable, z = df_clr$study_numeric, method = "spearman")
    
  }, error = function(e) {
    # If the test fails, return a list with NA values to prevent the loop from stopping.
    warning("Could not compute correlation for '", cell_type, "'. Error: ", e$message)
    return(list(estimate = NA, p.value = NA))
  })
  
  # Store the results for the current cell type in the list.
  # We create a small data frame for each result.
  results_list[[cell_type]] <- data.frame(
    cell_type = cell_type,
    partial_spearman_R = test_result$estimate,
    p_value = test_result$p.value
  )
}

# 5. Combine the list of data frames into a single, final data frame.
# bind_rows is an efficient way to stack the data frames from the list.
correlation_results <- bind_rows(results_list)

# Reset row names to be sequential (1, 2, 3, ...)
rownames(correlation_results) <- NULL

# --- Display the Final Dataframe ---
print("Partial Spearman Correlation Results (Controlling for Study):")
print(correlation_results)

#volcano plot
correlation_results <- correlation_results %>%
  mutate(
    log10_p = -log10(p_value),
    significant = p_value < 0.05
  )

ggplot(correlation_results, aes(x = partial_spearman_R, y = log10_p)) +
  geom_point(aes(color = significant), size = 2) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey")) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted") +
  theme_bw() +
  labs(
    x = "Partial Spearman R",
    y = "-log10(p-value)",
    title = "Partial Spearman Correlations"
  ) +
  # Add labels only for significant points
  geom_text_repel(
    data = correlation_results %>% filter(significant),
    aes(label = cell_type),
    size = 3,
    max.overlaps = 20
  )

# STING Niche Cox Model####
# cox model
table(df_clr$subtype, df_clr$study)

df_idh_wt <- df_clr %>%
  filter(IDH_status == "wt") %>%
  filter(platform == "rnaseq") %>%
  filter(is.na(location)) %>%
  mutate(
    mes_group = ifelse(subtype == "MES_like_Cancer", "MES_like", "not_MES_like"),
    niche_all_quartile = ntile(niche_all, 4)  # Creates 4 quartile groups (1 = lowest, 4 = highest)
  )

model_formula <- as.formula(
  paste("Surv(", time_column, ",", event_column, ") ~ niche_all + age + KarnPerfScore + strata(study)")
)

cox_model <- coxph(model_formula, data = df_idh_wt)
summary(cox_model)

gg_niche <- ggpredict(cox_model, terms = "niche_all [all]")
ggplot(gg_niche, aes(x = x, y = predicted)) +
  geom_line(color = "#0072B2", size = 1.2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, fill = "#0072B2") +
  labs(
    x = "STING Niche Enrichment",
    y = "Hazard Ratio",
    title = "Cox Model: Hazard Ratio across niche_all"
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(hjust = c(0, 1))
  )

#subtype interaction term
model_formula <- as.formula(
  paste("Surv(", time_column, ",", event_column, ") ~ niche_all * MES_like_Cancer_clr + strata(study)")
)

cox_model <- coxph(model_formula, data = df_idh_wt)
summary(cox_model)


# STING niche KM Survival Curves by subtype####
#km curve, mes-like
df_km <- df_clr %>%
  filter(IDH_status == "wt") %>%
  filter(study == "CGGA Batch 2") %>%
  mutate(
    mes_group = ifelse(subtype == "MES_like_Cancer", "MES_like", "not_MES_like")
  )

df_km_mes <- df_km %>%
  #filter(mes_group == "MES_like") %>%
  mutate(
    niche_all_group = ntile(niche_all, 2),  # 1 = lowest third, 3 = highest third
  ) %>%
  mutate(OS_years = as.numeric(OS_years)) %>%
  mutate(OS_event = as.numeric(OS_event))

surv_fit_object <- survfit(Surv(OS_years, OS_event) ~ niche_all_group, data = df_km_mes)

ggsurvplot(
  surv_fit_object,
  data = df_km_mes,
  title = "MES-like, TCGA",
  xlab = "Years",
  legend.title = "Niche Enrichment",
  legend.labs = c("Low", "High"),
  pval = TRUE,             
  conf.int = TRUE,         
  risk.table = TRUE,      
  ggtheme = theme_minimal()
)

#km curve, not mes-like
df_km_notmes <- df_km %>%
  filter(mes_group == "not_MES_like") %>%
  mutate(
    niche_all_group = ntile(niche_all, 2),  
  ) %>%
  mutate(OS_years = as.numeric(OS_years)) %>%
  mutate(OS_event = as.numeric(OS_event))

surv_fit_object <- survfit(Surv(OS_years, OS_event) ~ niche_all_group, data = df_km_notmes)

ggsurvplot(
  surv_fit_object,
  data = df_km_notmes,
  title = "Not MES like, TCGA",
  xlab = "Years",
  legend.title = "Niche Enrichment",
  legend.labs = c("Low", "High"),
  pval = TRUE,            
  conf.int = TRUE,         
  risk.table = TRUE,      
  ggtheme = theme_minimal() 
)


