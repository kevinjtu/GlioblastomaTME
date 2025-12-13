#libraries####
library(data.table)
library(dplyr)
library(tidyr)
library(tibble)
library(msigdbr)
library(ggplot2)
library(ggpubr)
library(scales)
library(mediation)
library(purrr)
library(openxlsx)
library(compositions)

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
  ) %>%
  
  # 3) drop the original tumor‐cell columns
  dplyr::select(-all_of(tumor_cells))

# sanity check: the non‐tumor fractions should now sum to 1
data2 %>%
  mutate(check = rowSums(across(`Oligodendrocyte`:`Radial_Glia`))) %>%
  dplyr::select(sample_id, tumor_purity, check) %>%
  head(10)

data2 <- as.data.frame(data2)

df <- data2 %>%
  filter(recurrence == "primary")

cell_cols <- names(df)[ which(names(df) == "Oligodendrocyte"):
                          which(names(df) == "Radial_Glia") ]

#CLR Transformation####
# 1) extract the raw proportions matrix
df_cells <- df %>% dplyr::select(all_of(cell_cols))

# 2) compute the CLR
clr_mat <- clr(as.matrix(df_cells))

# 3) rename columns to mark them as CLR‐transformed#
colnames(clr_mat) <- paste0(cell_cols, "_clr")

# 4) rebuild your analysis data frame, replacing raw with CLR
df_clr <- df %>%
  dplyr::select(-all_of(cell_cols)) %>%
  bind_cols(as.data.frame(clr_mat))

# 5) update the list of columns you’ll loop over
cell_cols_clr <- paste0(cell_cols, "_clr")


# ComBat batch correction####
df_rnaseq <- df_clr %>%
  filter(platform == "rnaseq")

start_col <- "Oligodendrocyte"
end_col <- "Radial_Glia"
cell_type_cols <- names(df_rnaseq)[which(names(df_rnaseq) == start_col):which(names(df_rnaseq) == end_col)]

# ComBat requires a matrix where features (cell types) are rows and samples are columns.
expression_matrix <- t(df_rnaseq[, cell_type_cols])
batch_info <- df_rnaseq$study

#Running ComBat
modcombat <- model.matrix(~1, data = df_clr_rnaseq)
combat_corrected_matrix <- ComBat(dat = expression_matrix, 
                                  batch = batch_info, 
                                  mod = modcombat, 
                                  par.prior = TRUE)

combat_corrected_df <- as.data.frame(t(combat_corrected_matrix))

metadata <- df_rnaseq %>%
  dplyr::select(-one_of(cell_type_cols))
df_corrected <- cbind(metadata, combat_corrected_df)

# --- 5. Verification ---
cat("\nDimensions of original filtered data:", dim(df_clr_rnaseq), "\n")
cat("Dimensions of corrected data:", dim(df_clr_corrected), "\n")


# heatmap matrix for Morpheus####
df_gsam <- df_clr %>%
  filter(study == "TEMPUS") %>%
  dplyr::select(patient_id, tumor_purity, all_of(cell_cols_clr))

tumor_purity <- df_gsam %>%
  dplyr::select(patient_id, tumor_purity)

cell_scaled <- df_gsam %>%
  dplyr::select(patient_id, all_of(cell_cols_clr)) %>%
  column_to_rownames("patient_id")

cell_scaled <- scale(cell_scaled)

matrix <- cbind(
  tumor_purity = tumor_purity$tumor_purity[match(rownames(cell_scaled), tumor_purity$patient_id)],
  cell_scaled
)
matrix <- t(matrix)

# Create a clean version of rownames without "_clr"
clean_rownames <- gsub("_clr$", "", rownames(matrix))

#mark lineage
# Define the lineage mapping
lineages <- c(
  Subtype                         = NA,
  CD8_Effector_Memory_T_Cells     = "Lymphoid",
  Pro_inflammatory_Microglia      = "Myeloid",
  Radial_Glia                     = "Normal Tissue",
  Astrocyte                       = "Normal Tissue",
  Arterial_Vessel                 = "Endothelial",
  Hypoxia_associated_Monocytes    = "Myeloid",
  OPC                             = "Normal Tissue",
  Pericyte                        = "Endothelial",
  Mast_Cell                       = "Myeloid",
  Cytotoxic_CD8_T_Cell            = "Lymphoid",
  Microglia_Aging_Signature       = "Myeloid",
  Naive_Monocyte                  = "Myeloid",
  Neuron                          = "Normal Tissue",
  Plasma_Cell                     = "Lymphoid",
  T_Cells_Stress_Signature        = "Lymphoid",
  Anti_inflammatory_Monocyte      = "Myeloid",
  BDM_Hypoxia_and_MES_like        = "Myeloid",
  Anti_inflammatory_BDM           = "Myeloid",
  Natural_Killer_Cells            = "Lymphoid",
  BDM_IFN_Signature               = "Myeloid",
  Plasmacytoid_DC                 = "Myeloid",
  B_Cell                          = "Lymphoid",
  Proliferating_T_Cells           = "Lymphoid",
  Capilary_Vessel                 = "Endothelial",
  Tiplike_Vessel                  = "Endothelial",
  Conventional_DC_1               = "Myeloid",
  CD4_T_Cell_IFN_Signature        = "Lymphoid",
  Proliferative_Microglia         = "Myeloid",
  BDM_MHC_Expression              = "Myeloid",
  Dendritic_like_Cell             = "Myeloid",
  NK_like_CD8_T_Cells             = "Lymphoid",
  Resting_CD4_T_Cells             = "Lymphoid",
  Conventional_DC_2               = "Myeloid",
  Regulatory_T_Cells              = "Lymphoid",
  Oligodendrocyte                 = "Normal Tissue"
)

# Assign new rownames and extract corresponding lineages
rownames(matrix) <- clean_rownames
lineage_column <- lineages[rownames(matrix)]

# Bind the lineage as the first column
matrix_annotated <- cbind(Lineage = lineage_column, matrix)

# Optional: view result
head(matrix_annotated)
rownames(matrix_annotated)[1] <- "Tumor Cellularity"

setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Heatmap - TME Signatures")
write.csv(matrix_annotated, "Morpheus_GBM_Ecotype_Heatmap.csv")
