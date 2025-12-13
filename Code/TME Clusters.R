#libraries####
library(data.table)
library(dplyr)
library(forcats)
library(compositions)
library(textshape)
library(NbClust)
library(survival)
library(survminer)
library(readxl)
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
                          which(names(df) == "Natural_Killer_Cells") ]

#CLR Transformation####
# 1) extract the raw proportions matrix
df_cells <- df %>% dplyr::select(all_of(cell_cols))

# 2) add a tiny pseudocount to avoid zeros
#pcount <- min(df_cells[df_cells > 0], na.rm = TRUE) * 0.5
#df_cells[df_cells == 0] <- pcount

# 3) compute the CLR
clr_mat <- clr(as.matrix(df_cells))

# 4) rename columns to mark them as CLR‐transformed#
colnames(clr_mat) <- paste0(cell_cols, "_clr")

# 5) rebuild your analysis data frame, replacing raw with CLR
df_clr <- df %>%
  dplyr::select(-all_of(cell_cols)) %>%
  bind_cols(as.data.frame(clr_mat))

# 6) update the list of columns you’ll loop over
cell_cols_clr <- paste0(cell_cols, "_clr")

#heatmap matrix for Morpheus####
df_gsam <- df_clr %>%
  filter(study == "TCGA", IDH_status == "wt", recurrence == "primary") %>%
  dplyr::select(patient_id, tumor_purity, all_of(cell_cols_clr))
tumor_purity <- df_gsam %>%
  dplyr::select(patient_id, tumor_purity)
cell_scaled <- df_gsam %>%
  dplyr::select(patient_id, all_of(cell_cols_clr)) %>%
  column_to_rownames("patient_id") %>%
  as.matrix() %>%
  scale()
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
write.csv(matrix_annotated, "Morpheus_GBM_TME_Heatmap.csv")

#find optimal number of clusters - patients
d <- df_clr %>%
  filter(study == "TCGA") %>%
  dplyr::select(patient_id, all_of(cell_cols_clr)) %>%
  column_to_rownames("patient_id") %>%
  scale() #convert to z scores

# List of clustering indices
index_list <- c("kl", "ch", "hartigan", "cindex", "db", "silhouette", 
                "duda", "pseudot2", "beale", "ratkowsky", "ball", "ptbiserial")
results <- list()
for (index in index_list) {
  message("Running index: ", index)
  try({
    res <- NbClust(d, distance = "euclidean", min.nc = 2, max.nc = 15, 
                   method = "kmeans", index = index)
    
    best_nc <- if (is.null(res$Best.nc)) NA else res$Best.nc[1]
    results[[index]] <- best_nc
  }, silent = TRUE)
}
cluster_results <- tibble(
  index = names(results),
  best_n_clusters = unlist(results)
)

print(cluster_results)

#find optimal number of clusters - cell clusters
d <- df_clr %>%
  filter(study == "GSAM") %>%
  dplyr::select(patient_id, all_of(cell_cols_clr)) %>%
  column_to_rownames("patient_id") %>%
  scale() %>% #convert to z scores
  t()

res <- NbClust(d, distance = "euclidean", min.nc=2, max.nc=15, 
               method = "kmeans", index = "sdbw")
res$Best.nc


#4 clusters optimal: use k-means method with db, ptbiserial, sdindex as indices. Distance is euclidean. 100,000 iterations

#Cluster kaplan meier curves with survival####
setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Heatmap - TME Signatures")
clusters <- read_excel("GBM TME Cluster Assignments.xlsx")

subtype_cols <- c("OPC_like_Cancer", "AC_like_Cancer", "NPC_like_Cancer", "MES_like_Cancer")

data <- data %>%
  mutate(
    subtype = subtype_cols[
      max.col( as.matrix( across(all_of(subtype_cols)) ),
               ties.method = "first")])

gsam_clusters <- data %>%
  filter(study == "TCGA", recurrence == "primary", subtype == "MES_like_Cancer") %>%
  dplyr::select(patient_id, OS_event, OS_years, DFS_event, DFS_years, subtype) %>%
  left_join(clusters, by = "patient_id") %>%
  dplyr::mutate(TME_signature = factor(TME_signature), subtype = factor(subtype))

#overal survival plot
fit_os <- survfit(
  Surv(OS_years, OS_event) ~ TME_signature,
  data = gsam_clusters)

os_plot <- ggsurvplot(
  fit_os,
  data         = gsam_clusters,
  palette      = "Set1",
  
  legend.title = "TME Cluster",
  xlab         = "Time (years)",
  ylab         = "Overall Survival Probability",
  pval         = TRUE,        # show log-rank p-value
  pval.method  = FALSE,       # suppress "Log-rank" label
  pval.size    = 4)

print(os_plot)