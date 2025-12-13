#libraries####
library(data.table)
library(dplyr)
library(tidyr)
library(tibble)
library(msigdbr)
library(ggplot2)
library(ggpubr)
library(scales)
library(compositions)
library(mclust)
library(GSVA)

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


#separate samples by low & high STING reactome enrichment####
setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Expression Data")
df <- fread("IVYGAP fpkm.csv")
setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/IVYGAP")
genes <- fread("rows-genes.csv")
genes <- genes[,c(1,4)]
df[,1] <- genes$gene_symbol
colnames(df)[1] <- "Hugo_Symbol"

df <- df %>%
  filter(!is.na(Hugo_Symbol)) %>%
  distinct(Hugo_Symbol, .keep_all = TRUE)
df <- column_to_rownames(df, var = "Hugo_Symbol")

df <- data.frame(lapply(df, function(x) as.numeric(as.character(x))),
                 row.names = rownames(df))

tpm_df <- as.data.frame(apply(df, 2, function(col) {
  (col / sum(col, na.rm = TRUE)) * 1e6
}))
colSums(tpm_df)

tpm_df <- as.matrix(tpm_df)

all_c2 <- msigdbr(species = "Homo sapiens", collection = "C2")
reactome_sting_genes <- all_c2 %>%
  filter(gs_name == "REACTOME_STING_MEDIATED_INDUCTION_OF_HOST_IMMUNE_RESPONSES") %>%
  pull(gene_symbol)

param_react   <- ssgseaParam(exprData = tpm_df,
                             geneSets = list(REACTOME_STING = reactome_sting_genes),
                             normalize = TRUE)

res_reactome <- gsva(param_react)

res_reactome_t <- t(res_reactome) %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "sample_id")
res_reactome_t$sample_id <- sub("^X", "", res_reactome_t$sample_id)

STING_clusters <- Mclust(res_reactome_t$REACTOME_STING, G = 2, modelNames = "V")
STING_clusters_df <- data.frame(
  sample_id = res_reactome_t$sample_id,
  STING_cluster = STING_clusters$classification
)

res_reactome_t <- merge(res_reactome_t, STING_clusters_df, by = "sample_id")

df_clr <- df_clr %>%
  filter(!is.na(location)) %>%
  merge(res_reactome_t, by = "sample_id") %>%
  group_by(subtype) %>%
  mutate(
    STING_status = if_else(
      REACTOME_STING >= median(REACTOME_STING, na.rm = TRUE),
      "High STING",
      "Low STING"
    )
  ) %>%
  ungroup()
  
# Matricies for heatmaps####
location_order <- c(
  "Perinecrotic Zone", "Pseudopalisading Cells", "Cellular Tumor",
  "Infiltrating Tumor", "Leading Edge")

# 2. Define cell type columns (adjust this if needed)
cell_cols <- names(df_clr)[grepl("_clr$", names(df_clr))]

# 3. Define a function to create matrix
make_heatmap_matrix <- function(data, subtype_filter, sting_status_filter) {
  data %>%
    filter(
      subtype == subtype_filter,
      STING_status == sting_status_filter,
      location %in% location_order
    ) %>%
    group_by(location) %>%
    summarise(across(all_of(cell_cols), mean, na.rm = TRUE), .groups = "drop") %>%
    mutate(location = factor(location, levels = location_order)) %>%
    arrange(location) %>%
    column_to_rownames("location") %>%
    t()  # Transpose to have cell types as rows and locations as columns
}

# 4. Generate all 4 matrices
heatmap_high_mes <- make_heatmap_matrix(df_clr, "MES_like_Cancer", "High STING")
heatmap_low_mes  <- make_heatmap_matrix(df_clr, "MES_like_Cancer", "Low STING")
heatmap_high_npc <- make_heatmap_matrix(df_clr, "NPC_like_Cancer", "High STING")
heatmap_low_npc  <- make_heatmap_matrix(df_clr, "NPC_like_Cancer", "Low STING")

setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/")

write.csv(heatmap_high_mes, "heatmap_high_STING_MES.csv")
write.csv(heatmap_low_mes,  "heatmap_low_STING_MES.csv")
write.csv(heatmap_high_npc, "heatmap_high_STING_NPC.csv")
write.csv(heatmap_low_npc,  "heatmap_low_STING_NPC.csv")
