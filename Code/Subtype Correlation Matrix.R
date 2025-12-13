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
library(corrplot)

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

#Spearman correlation between cell types####
cell_type_columns <- names(df_clr)[
  which(names(df_clr) == "Oligodendrocyte"):
    which(names(df_clr) == "Radial_Glia")
]
subtypes <- unique(df_clr$subtype)

final_corr_matrices <- list()

for (st in subtypes) {
  message(paste("Processing subtype:", st))
  df_subtype <- df_clr %>% filter(subtype == st)
  studies_in_subtype <- unique(df_subtype$study)
  study_corr_matrices <- list()
  
  for (s in studies_in_subtype) {
    df_study <- df_subtype %>% filter(study == s)
    
    if (nrow(df_study) > 1) {
      cor_matrix <- cor(df_study[, cell_type_columns], method = "spearman", use = "pairwise.complete.obs")
      study_corr_matrices[[s]] <- cor_matrix
    }
  }
  
  if (length(study_corr_matrices) > 0) {
    avg_corr_matrix <- Reduce("+", study_corr_matrices) / length(study_corr_matrices)
    final_corr_matrices[[st]] <- avg_corr_matrix
  }
}

message("Determining consistent cell type order...")
if (length(final_corr_matrices) > 0) {
  reference_subtype <- names(final_corr_matrices)[1]
  reference_matrix <- final_corr_matrices[[reference_subtype]]
  hclust_order <- hclust(as.dist(1 - reference_matrix))$order
  ordered_cell_types <- colnames(reference_matrix)[hclust_order]
} else {
  ordered_cell_types <- cell_type_columns
  message("Warning: No correlation matrices were generated. Using default order.")
}

message("Plotting correlation matrices for each subtype separately...")

for (st in names(final_corr_matrices)) {
  corr_matrix_to_plot <- final_corr_matrices[[st]]
  corr_matrix_ordered <- corr_matrix_to_plot[ordered_cell_types, ordered_cell_types]
  
  corrplot(corr_matrix_ordered,
           method = "color",
           type = "upper",
           order = "original",
           tl.col = "black",
           tl.cex = 0.7,
           #tl.pos = "l",         # Only show labels on the left
           #cl.lim = c(-1, 1),
           addgrid.col = NA,
           title = st,
           mar = c(0, 0, 2, 0)
  )
}