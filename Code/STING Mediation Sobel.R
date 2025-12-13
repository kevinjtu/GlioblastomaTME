#libraries####
library(data.table)
library(dplyr)
library(tidyr)
library(tibble)
library(msigdbr)
library(ggplot2)
library(ggpubr)
library(scales)
library(GSVA)
library(GSEABase)
library(mediation)
library(purrr)
library(openxlsx)

#import TEMPUS data####
#expression data
setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Expression Data")
TEMPUS <- as.data.frame(fread("TEMPUS RNAseq expression.csv")) %>%
  column_to_rownames(var = "V1")
colnames(TEMPUS) <- TEMPUS[1, ]
TEMPUS <- TEMPUS[-1, ]

#tme signature data
setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Heatmap - TME Signatures")
TME_sig <- as.data.frame(fread("GBM TME Classifier Results.csv")) %>% dplyr::select(-V1)
colnames(TME_sig) <- c("patient_id", "TME_signature")
TME_sig <- TME_sig[-1,]

#tumor subtypes
subtype_cols <- c("OPC_like_Cancer", "AC_like_Cancer", "NPC_like_Cancer", "MES_like_Cancer")
setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME")
tumor_subtype <- as.data.frame(fread("GBM Clinical-TME Cell States.csv")) %>%
  mutate(subtype = subtype_cols[
    max.col( as.matrix( across(all_of(subtype_cols)) ),
             ties.method = "first")]) %>%
  dplyr::select(patient_id, subtype)

#metadata
meta <- tumor_subtype %>%
  inner_join(TME_sig, by = "patient_id") %>%      # 1) merge
  distinct(patient_id, .keep_all = TRUE) %>%       # 2) keep one row per patient_id
  column_to_rownames(var = "patient_id")  

# make sure samples line up
common_samples <- intersect(colnames(TEMPUS), rownames(meta))
counts <- TEMPUS[, common_samples]
meta   <- meta[ common_samples, , drop=FALSE ]
counts_mat <- as.matrix(counts)

meta <- meta %>%
  mutate(immune_active = ifelse(TME_signature == "Immunoactive", 1, 0))

stopifnot(all(colnames(counts_mat) %in% rownames(meta)))
s
#determining STING and IFNa/y scores####
#Composite STING Score
all_c2 <- msigdbr(species = "Homo sapiens", category = "C2")
cgas_sting_genes <- all_c2 %>%
  filter(gs_name == "KEGG_MEDICUS_REFERENCE_CGAS_STING_SIGNALING_PATHWAY") %>%
  pull(gene_symbol)

reactome_sting_genes <- all_c2 %>%
  filter(gs_name == "REACTOME_STING_MEDIATED_INDUCTION_OF_HOST_IMMUNE_RESPONSES") %>%
  pull(gene_symbol)

param_cgas    <- ssgseaParam(exprData = counts_mat,
                             geneSets = list(CGAS_STING = cgas_sting_genes),
                             normalize = TRUE)             # ssGSEA normalization :contentReference[oaicite:0]{index=0}
param_react   <- ssgseaParam(exprData = counts_mat,
                             geneSets = list(REACTOME_STING = reactome_sting_genes),
                             normalize = TRUE)

res_cgas  <- gsva(param_cgas)                                         # note the new API: gsva(param) :contentReference[oaicite:1]{index=1}
res_reactome <- gsva(param_react)

#IFNa/y gene set
hallmark_h <- msigdbr(species = "Homo sapiens", category = "H")
ifna_genes  <- hallmark_h %>% filter(gs_name == "HALLMARK_INTERFERON_ALPHA_RESPONSE") %>% pull(gene_symbol)
ifng_genes  <- hallmark_h %>% filter(gs_name == "HALLMARK_INTERFERON_GAMMA_RESPONSE") %>% pull(gene_symbol)

param_ifna     <- ssgseaParam(exprData = counts_mat,
                              geneSets = list(IFNA                = ifna_genes),
                              normalize = TRUE)
param_ifng     <- ssgseaParam(exprData = counts_mat,
                              geneSets = list(IFNG                = ifng_genes),
                              normalize = TRUE)

res_ifna     <- gsva(param_ifna)
res_ifng     <- gsva(param_ifng)

#put into meta
meta <- meta %>%
  mutate(
    KEGG_STING_score          = as.numeric(res_cgas["CGAS_STING",          rownames(meta)]),
    REACTOME_STING_score      = as.numeric(res_reactome["REACTOME_STING",   rownames(meta)]),
    IFNA_score                = as.numeric(res_ifna["IFNA",                   rownames(meta)]),
    IFNG_score                = as.numeric(res_ifng["IFNG",                   rownames(meta)])
  )

#Sobel mediating analysis####
#this tests whether immune_active → STING_score → IFN_score (i.e. whether the difference in IFN between “immunosuppressed” and “immunoactive” samples is mediated by STING.)
# 1. Define your Sobel‐test function (same as before)
run_sobel_flipped <- function(df, subtype, treat, mediator, outcome) {
  # mediator model: STING_score ~ immune_active
  med_mod <- lm(as.formula(paste(mediator, "~", treat)), data = df)
  
  # outcome model: IFN_score ~ immune_active + STING_score
  out_mod <- lm(as.formula(paste(outcome, "~", treat, "+", mediator)), data = df)
  
  # extract paths & SEs
  a  <- coef(summary(med_mod)) [treat,     "Estimate"]
  sa <- coef(summary(med_mod)) [treat,     "Std. Error"]
  b  <- coef(summary(out_mod))[mediator,  "Estimate"]
  sb <- coef(summary(out_mod))[mediator,  "Std. Error"]
  
  # Sobel z & two‐tailed p
  z  <- a * b / sqrt(b^2 * sa^2 + a^2 * sb^2)
  p  <- 2 * (1 - pnorm(abs(z)))
  
  tibble(
    subtype   = subtype,
    treat     = treat,
    mediator  = mediator,
    outcome   = outcome,
    a         = a,
    b         = b,
    sobel_z   = z,
    p.value   = p
  )
}

# 2. Define your parameters
subtypes   <- unique(meta$subtype)
sting_scores <- c("KEGG_STING_score", "REACTOME_STING_score")
ifn_scores   <- c("IFNA_score",       "IFNG_score")
treat_var    <- "immune_active"

# 3. Loop over every subtype × STING × IFN combo (4 × 2 × 2 = 16)
results_list <- list()
ctr <- 1
for (sub in subtypes) {
  df_sub <- subset(meta, subtype == sub)
  for (med in sting_scores) {
    for (out in ifn_scores) {
      results_list[[ctr]] <- run_sobel_flipped(
        df       = df_sub,
        subtype  = sub,
        treat    = treat_var,
        mediator = med,
        outcome  = out
      )
      ctr <- ctr + 1
    }
  }
}

# 4. Combine into one tibble and view
sobel_results <- bind_rows(results_list)
print(sobel_results)
#prepare morpheus-ready heatmap####
filtered <- sobel_results %>%
  filter(mediator != "KEGG_STING_score")

filtered <- filtered %>%
  mutate(indirect = a * b)

# 3. Build wide matrix of p‐values
pvals_wide <- filtered %>%
  dplyr::mutate(log10p = -log10(p.value)) %>%
  dplyr::select(subtype, outcome, log10p) %>%
  pivot_wider(
    names_from  = outcome,
    values_from = log10p
  )


# 4. Build wide matrix of indirect effects
effects_wide <- filtered %>%
  dplyr::select(subtype, outcome, indirect) %>%
  pivot_wider(
    names_from  = outcome,
    values_from = indirect
  )


# 5. Write to Excel with two sheets
setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Sobel Mediating Analysis")
wb <- createWorkbook()
addWorksheet(wb, "p_values")
writeData(wb, "p_values", pvals_wide)
addWorksheet(wb, "indirect_effects")
writeData(wb, "indirect_effects", effects_wide)
saveWorkbook(wb, "sobel_summary.xlsx", overwrite = TRUE)

# 6. Confirm
pvals_wide
effects_wide