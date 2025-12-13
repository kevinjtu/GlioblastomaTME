#libraries####
library(data.table)
library(dplyr)
library(survival)
library(ggplot2)
library(forcats)
library(compositions)

#Set  survival metric here: "OS", "PFS", or "DFS"
Surv_metric <- "OS"

#import & clean data####
subtype_cols <- c("OPC_like_Cancer", "AC_like_Cancer", "NPC_like_Cancer", "MES_like_Cancer")
setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME")
data <- as.data.frame(fread("GBM Clinical-TME Cell States.csv")) %>%
  mutate(subtype = subtype_cols[
    max.col( as.matrix( across(all_of(subtype_cols)) ),
             ties.method = "first")])

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
      `Oligodendrocyte`:`Natural_Killer_Cells`,
      ~ . / (1 - tumor_purity)
    )
  ) %>%
  
  # 3) drop the original tumor‐cell columns
  dplyr::select(-all_of(tumor_cells))

# sanity check: the non‐tumor fractions should now sum to 1
data2 %>%
  mutate(check = rowSums(across(`Oligodendrocyte`:`Natural_Killer_Cells`))) %>%
  dplyr::select(sample_id, tumor_purity, check) %>%
  head(10)

data2 <- as.data.frame(data2)

df <- data2 %>%
  filter(recurrence == "primary")

cell_cols <- names(df)[ which(names(df) == "Oligodendrocyte"):
                          which(names(df) == "OPC") ]

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

#Cox Model####
df_idhwt <- df_clr %>%
  filter(IDH_status == "wt",
         !is.na(OS_years),
         !is.na(OS_event),
         !is.na(STING1_tpm),
         !is.na(subtype))

df_idhwt <- df_idhwt %>%
  mutate(subtype_recoded = ifelse(subtype == "MES_like_Cancer", "MES_like_Cancer", "Not-MES"))
df_idhwt$subtype_recoded <- factor(df_idhwt$subtype_recoded, levels = c("Not-MES", "MES_like_Cancer"))

if (nrow(df_idhwt) > 0 && length(unique(df_idhwt$subtype_recoded)) > 1) {
  # The reference level was set during factor creation above.
  # If "Not-MES" is the reference, the interaction term will be STING1_tpm:subtype_recodedMES_like_Cancer
  cat("Reference level for subtype_recoded in interaction model is:", levels(df_idhwt$subtype_recoded)[1], "\n")
  
  tryCatch({
    cox_model_interaction <- coxph(
      Surv(OS_years, OS_event) ~ REACTOME_STING * subtype_recoded + KarnPerfScore + strata(study), # Using subtype_recoded here
      data = df_idhwt
    )
    
    cat("\nSummary of the Cox model with interaction term (using recoded subtype):\n")
    summary_interaction_model <- summary(cox_model_interaction)
    print(summary_interaction_model)
    
    # Interpretation for Approach 2 (with recoded subtypes):
    # - The coefficient for 'STING1_tpm' is its effect in the *reference* recoded subtype (e.g., "Not-MES").
    # - The coefficient for 'subtype_recodedMES_like_Cancer' (if "Not-MES" is ref) represents the
    #   difference in the log-hazard for "MES_like_Cancer" compared to "Not-MES",
    #   when STING1_tpm is zero.
    # - The coefficient for 'STING1_tpm:subtype_recodedMES_like_Cancer' is the key interaction term.
    #   A significant p-value for this term means that the effect of STING1_tpm on survival
    #   is significantly *different* in 'MES_like_Cancer' compared to its effect in "Not-MES".
    #
    # To achieve your goal:
    # - Look for a significant p-value for 'STING1_tpm:subtype_recodedMES_like_Cancer'.
    # This would indicate that the relationship between STING1_tpm and survival specifically
    # changes for the MES_like_Cancer subtype relative to all other combined subtypes ("Not-MES").
    
  }, error = function(e) {
    cat("Could not fit interaction model with recoded subtype:", e$message, "\n")
  })
} else {
  cat("Not enough data or distinct recoded subtypes to fit an interaction model after filtering.\n",
      "Need at least two recoded subtype levels for interaction.\n",
      "Samples:", nrow(df_idhwt), "Unique recoded subtypes:", length(unique(df_idhwt$subtype_recoded)), "\n")
}