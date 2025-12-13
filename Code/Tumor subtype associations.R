# libraries ####
library(compositions)   
library(dplyr)
library(tidyr)
library(pheatmap)
library(ggplot2)
library(ALDEx2)        
library(stats)
library(purrr)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(openxlsx)

# 1. Import data & DEFINE YOUR COLUMN SETS ---------------------------------------------
setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME")
data <- as.data.frame(fread("GBM Clinical-TME Cell States.csv"))

tumor.cols <- c("OPC_like_Cancer","AC_like_Cancer","NPC_like_Cancer","MES_like_Cancer")

# select all other columns that are TME (here by negative selection; adjust as needed)
tme.cols   <- setdiff(colnames(data), c(tumor.cols, 
                                        # non-cell metadata
                                        "sample_id","patient_id","age","sex","race","grade",
                                        "OS_event","OS_years","DFS_event","DFS_years",
                                        "PFS_event","PFS_years","recurrence","location",
                                        "IDH_status","study","ARID1A_tpm","STING1_tpm"
))

tumor.mat <- data[tumor.cols]
tme.mat   <- data[tme.cols]

# 2. CLR TRANSFORM -------------------------------------
clr.tumor <- as.data.frame(clr(acomp(tumor.mat)))
colnames(clr.tumor) <- paste0("clr_", tumor.cols)

clr.tme   <- as.data.frame(clr(acomp(tme.mat)))
colnames(clr.tme)   <- paste0("clr_", tme.cols)

data$dominant <- tumor.cols[ max.col(data[tumor.cols], ties.method = "first") ]

set.seed(42)
km <- kmeans(data[tumor.cols], centers = 4, nstart = 50)
data$cluster  <- factor(km$cluster)

# combine back into a single df
df.clr <- bind_cols(data, clr.tumor, clr.tme)

# 3. PAIRWISE Spearman CORRELATIONS -----------------------------------------------
cors <- expand.grid(tumor = tumor.cols, tme = tme.cols, stringsAsFactors = FALSE) %>%
  rowwise() %>%
  mutate(
    rho = cor(
      df.clr[[paste0("clr_", tumor)]],
      df.clr[[paste0("clr_", tme)]],
      method = "spearman"
    ),
    pval = cor.test(
      df.clr[[paste0("clr_", tumor)]],
      df.clr[[paste0("clr_", tme)]],
      method = "spearman"
    )$p.value
  ) %>%
  ungroup() %>%
  mutate(padj = p.adjust(pval, method = "BH")) %>%
  arrange(padj)

# View top significant correlations
head(cors)

# 4. MULTIVARIATE LINEAR MODELS ------------------------------------------
# regress each TME cell on all four tumor‐subtypes + optional covariates
covs <- c("age","sex","race")    # adjust as desired
lm.results <- lapply(tme.cols, function(cell) {
  formula <- as.formula(
    paste0("clr_", cell, " ~ ",
           paste(paste0("clr_", tumor.cols), collapse = " + "),
           if(length(covs)>0) paste0(" + ", paste(covs, collapse=" + ")) else "")
  )
  fit <- lm(formula, data = df.clr)
  coefs <- summary(fit)$coefficients
  data.frame(
    tme = cell,
    term = rownames(coefs),
    estimate = coefs[, "Estimate"],
    pval      = coefs[, "Pr(>|t|)"],
    row.names = NULL
  )
}) %>% bind_rows() %>%
  group_by(term) %>%
  mutate(padj = p.adjust(pval, method="BH")) %>%
  ungroup()

# 5. DIFFERENTIAL‐ABUNDANCE TESTING ---------------------------------------
# e.g. Wilcoxon between clusters for each TME cell
da.wilcox <- expand.grid(
  cluster1 = levels(df.clr$cluster),
  cluster2 = levels(df.clr$cluster),
  tme      = tme.cols,
  stringsAsFactors = FALSE
) %>%
  # only unique pairs
  filter(cluster1 < cluster2) %>%
  # for each row, compute one p-value
  mutate(
    pval = pmap_dbl(
      list(cluster1, cluster2, tme),
      function(c1, c2, tme_cell) {
        x <- df.clr %>%
          filter(cluster == c1) %>%
          pull(paste0("clr_", tme_cell))
        y <- df.clr %>%
          filter(cluster == c2) %>%
          pull(paste0("clr_", tme_cell))
        if (length(x) > 0 && length(y) > 0) {
          wilcox.test(x, y)$p.value
        } else {
          NA_real_
        }
      }
    )
  ) %>%
  ungroup() %>%
  mutate(padj = p.adjust(pval, method = "BH"))

head(da.wilcox)

# or use ALDEx2 on the CLR‐matrix for a more compositional‐aware test
#aldex.res <- aldex.clr(data[tme.cols], data$cluster, mc.samples=128)
#aldex.tt   <- aldex.ttest(aldex.res)
#aldex.effect <- aldex.effect(aldex.res)

# 7. VISUALIZATION mean CLR(TME) --------------------------------------------------------
# (a) heatmap of mean CLR‐TME by dominant subtype
heat.df <- df.clr %>%
  group_by(dominant) %>%
  summarise_at(vars(starts_with("clr_")), mean) %>%
  column_to_rownames("dominant")
pheatmap(as.matrix(heat.df), scale = "row", 
         main = "Mean CLR(TME) per dominant subtype")

# (b) boxplot for a top TME cell
top.tme <- cors$tme[1]
ggplot(df.clr, aes(x = dominant, y = .data[[paste0("clr_", top.tme)]])) +
  geom_boxplot() +
  labs(y = paste0("CLR(", top.tme,")"), x = "Dominant subtype") +
  theme_minimal()

# 8. VISUALIZATION spearman correlations - Broad Morop#####
# (a) heatmap of mean CLR‐TME by dominant subtype
setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Heatmap - Subtype Associations")
rho_df <- as.data.frame(rho_mat) %>%
  rownames_to_column(var = "Subtype")
lineages <- c(
  Subtype                         = NA,  # blank under the “Subtype” column
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
wb <- createWorkbook()
addWorksheet(wb, "Spearman Associations")
writeData(
  wb, "Spearman Associations",
  x = as.list(names(rho_df)),
  startRow = 1, colNames = FALSE
)
writeData(
  wb, "Spearman Associations",
  x = as.list(lineages[names(rho_df)]),
  startRow = 2, colNames = FALSE
)
writeData(
  wb, "Spearman Associations",
  x = rho_df,
  startRow = 3, colNames = FALSE)
freezePane(wb, "Spearman Associations", firstActiveRow = 3)
saveWorkbook(wb, "spearman_subtype_associations_with_lineage.xlsx", overwrite = TRUE)

#p-values
sig_mat <- -log10(padj_mat)
as.data.frame(sig_mat) %>%
  rownames_to_column(var = "Subtype") %>%
  writexl::write_xlsx("spearman_subtype_associations_padj.xlsx")

#go to morpheus from BROAD institute, 1 minus pearson average heirarchivcal clustering for both rows and columns