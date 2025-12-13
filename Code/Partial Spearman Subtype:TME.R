#libraries####
library(data.table)
library(dplyr)
library(survival)
library(ggplot2)
library(forcats)
library(compositions)
library(ppcor)
library(openxlsx)

#import & clean data####
setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME")
data <- as.data.frame(fread("GBM Clinical-TME Cell States.csv"))

#get TME fractions accounting for cellularity
# define your six tumor‐cell columns
tumor_cells <- c(
  "OPC_like_Cancer",
  "AC_like_Cancer",
  "NPC_like_Cancer",
  "MES_like_Cancer"
)

df <- data %>%
  filter(recurrence == "primary") %>%
  filter(is.na(location)) %>%
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


#spearman correlation, partial  moeling####
# Select tumor and non-tumor cell type columns
tumor_cols <- c("OPC_like_Cancer_clr", "AC_like_Cancer_clr", 
                "NPC_like_Cancer_clr", "MES_like_Cancer_clr")

nontumor_cols <- names(df_clr)[which(names(df_clr) == "Oligodendrocyte_clr"):
                                 which(names(df_clr) == "Radial_Glia_clr")]

# Initialize results dataframe
results <- data.frame(
  Tumor = character(),
  NonTumor = character(),
  Rho = numeric(),
  Pval = numeric(),
  stringsAsFactors = FALSE
)

# Loop over all combinations
for (tumor in tumor_cols) {
  for (non in nontumor_cols) {
    # Drop NA rows
    complete_idx <- complete.cases(df_clr[, c(tumor, non, "study")])
    
    x <- df_clr[[tumor]][complete_idx]
    y <- df_clr[[non]][complete_idx]
    z <- as.data.frame(model.matrix(~ study, data = df_clr[complete_idx, ]))[ , -1, drop = FALSE]
    
    # Compute partial Spearman correlation
    test <- tryCatch(
      pcor.test(x, y, z, method = "spearman"),
      error = function(e) list(estimate = NA, p.value = NA)
    )
    
    # Store results
    results <- rbind(results, data.frame(
      Tumor = tumor,
      NonTumor = non,
      Rho = test$estimate,
      Pval = test$p.value,
      stringsAsFactors = FALSE
    ))
  }
}

# Sort results by p-value
results$FDR <- p.adjust(results$Pval, method = "BH")
results <- results[order(results$FDR), ]

#print out into morhpeus heatmap####
# 1. Remove "_clr" from Tumor and NonTumor names
results <- results %>%
  mutate(
    Tumor = str_remove(Tumor, "_clr$"),
    NonTumor = str_remove(NonTumor, "_clr$")
  )

# 2. Define lineages
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

# 3. Add lineage column
results$Lineage <- lineages[results$NonTumor]

# 4. Rho matrix
rho_mat <- acast(results, NonTumor ~ Tumor, value.var = "Rho")
rho_mat <- as.data.frame(rho_mat)
rho_mat <- cbind(Lineage = lineages[rownames(rho_mat)], rho_mat)
write.xlsx(rho_mat, "partial_spearman_rho.xlsx", rowNames = TRUE)

# 5. FDR matrix (as -log10)
fdr_mat_raw <- acast(results, NonTumor ~ Tumor, value.var = "FDR")
fdr_mat_log <- -log10(fdr_mat_raw)
fdr_mat_log[is.infinite(fdr_mat_log)] <- 300  # Cap extreme -log10(FDR)
fdr_mat_log <- as.data.frame(fdr_mat_log)
fdr_mat_log <- cbind(Lineage = lineages[rownames(fdr_mat_log)], fdr_mat_log)
write.xlsx(fdr_mat_log, "partial_spearman_fdr_log10.xlsx", rowNames = TRUE)
