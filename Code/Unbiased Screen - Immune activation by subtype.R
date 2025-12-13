#libraries####
library(data.table)
library(dplyr)
library(tibble)
library(limma)
library(msigdbr)
library(clusterProfiler)
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
expr_mat <- as.matrix(TEMPUS)

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

#integrate data
meta <- tumor_subtype %>%
  inner_join(TME_sig, by = "patient_id") %>%      # 1) merge
  distinct(patient_id, .keep_all = TRUE) %>%       # 2) keep one row per patient_id
  column_to_rownames(var = "patient_id")  

# make sure samples line up
common_samples <- intersect(colnames(TEMPUS), rownames(meta))
counts <- TEMPUS[, common_samples]
meta   <- meta[ common_samples, , drop=FALSE ]

#unbiased screen of genes comparing immune activation between subtypes####
#–– 2. make sure your expr_mat and meta line up
expr <- expr_mat[, rownames(meta), drop=FALSE]
meta2 <- meta[colnames(expr), , drop=FALSE]

#–– 3. create a combined factor: "Subtype_TME"
meta2$group <- factor(paste0(meta2$subtype, "_", meta2$TME_signature))
design <- model.matrix(~0 + group, data = meta2)
colnames(design) <- levels(meta2$group)

#–– 5. fit linear model
fit <- lmFit(expr, design)

#–– 6. define contrasts: for each subtype, Immuneactive – Immunosuppressed
cont.matrix <- makeContrasts(
  AC  = AC_like_Cancer_Immunoactive   - AC_like_Cancer_Immunosuppressed,
  MES = MES_like_Cancer_Immunoactive  - MES_like_Cancer_Immunosuppressed,
  NPC = NPC_like_Cancer_Immunoactive  - NPC_like_Cancer_Immunosuppressed,
  OPC = OPC_like_Cancer_Immunoactive  - OPC_like_Cancer_Immunosuppressed,
  levels = design
)

fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

#–– 7. extract topTables for each subtype
res_list <- lapply(colnames(cont.matrix), function(coef) {
  tt <- topTable(fit2, coef = coef, number = Inf, sort.by = "none")
  tt$gene    <- rownames(tt)
  tt$subtype <- coef
  tt
})
res_df <- do.call(rbind, res_list)

#–– 8. reshape so each gene’s adj.P.Val is side-by-side
wide_p <- res_df %>%
  dplyr::select(gene, subtype, adj.P.Val) %>%
  pivot_wider(names_from = subtype, values_from = adj.P.Val)

#–– 9. pick out genes significant in MES only (for example)
alpha <- 0.05

sig_MES_only <- wide_p %>%
  filter(
    MES <  alpha,
    AC  >  alpha,
    NPC >  alpha,
    OPC >  alpha
  ) %>%
  pull(gene)

length(sig_MES_only)  # how many…

#–– 10. you can repeat for any subtype (“AC only”, “NPC only”, etc.) by swapping the filter
sig_AC_only <- wide_p %>%
  filter(
    AC  < alpha,
    MES > alpha,
    NPC > alpha,
    OPC > alpha
  ) %>%
  pull(gene)

length(sig_AC_only)  # how many…

#–– 11. finally, inspect your candidates
sig_NPC_only <- wide_p %>%
  filter(
    AC  > alpha,
    MES > alpha,
    NPC < alpha,
    OPC > alpha
  ) %>%
  pull(gene)

length(sig_NPC_only)  # how many…

#–– 11. finally, inspect your candidates
sig_OPC_only <- wide_p %>%
  filter(
    AC  > alpha,
    MES > alpha,
    NPC > alpha,
    OPC < alpha
  ) %>%
  pull(gene)

length(sig_OPC_only)  # how many…

#reactome pathway enrichment####
# 2. pull down all Reactome gene‐sets from MSigDB
reactome2gene <- msigdbr(
  species     = "Homo sapiens",
  category    = "C2",
  subcategory = "CP:REACTOME"
) %>%
  dplyr::select(gs_name, gene_symbol) %>%
  rename(term = gs_name, gene = gene_symbol)

# 2. your subtype‐specific gene lists
sig_genes_list <- list(
  AC  = sig_AC_only,
  MES = sig_MES_only,
  NPC = sig_NPC_only,
  OPC = sig_OPC_only
)

# 3. run ORA without q‐value filtering, and no p‐value adjustment if you prefer
ora_results <- lapply(sig_genes_list, function(genes) {
  enricher(
    gene          = genes,
    TERM2GENE     = reactome2gene,
    pAdjustMethod = "none",    # skip p.adjust
    pvalueCutoff  = 1,         # keep all by raw p
    qvalueCutoff  = 1          # keep all by q
  )
})

# 4. pull out only those with raw p < 0.10
signif_ora <- lapply(ora_results, function(e) {
  if (is.null(e) || nrow(e@result)==0) return(NULL)
  e@result %>%
    filter(pvalue < 0.10) %>%
    arrange(pvalue)
})

dim(signif_ora$AC)
dim(signif_ora$MES)
dim(signif_ora$OPC)
dim(signif_ora$NPC)

#Sobel Mediation Screen####
# 0. Prep: convert your counts → matrix
sub_map <- c(
  MES = "MES_like_Cancer",
  OPC = "OPC_like_Cancer",
  NPC = "NPC_like_Cancer",
  AC  = "AC_like_Cancer"
)

outcome_map <- list(
  MES_like_Cancer = c("HALLMARK_INTERFERON_ALPHA_RESPONSE",
                      "HALLMARK_INTERFERON_GAMMA_RESPONSE"),
  OPC_like_Cancer = c("HALLMARK_ESTROGEN_RESPONSE_LATE",
                      "HALLMARK_ESTROGEN_RESPONSE_EARLY"),
  NPC_like_Cancer = c("HALLMARK_COMPLEMENT",
                      "HALLMARK_IL6_JAK_STAT3_SIGNALING"),
  AC_like_Cancer  = c("HALLMARK_MYC_TARGETS_V1",
                      "HALLMARK_MTORC1_SIGNALING")
)

#–– 2. Pull MSigDB once
all_c2     <- msigdbr(species = "Homo sapiens",
                      category = "C2",
                      subcategory = "CP:REACTOME")
hallmark_h <- msigdbr(species = "Homo sapiens", category = "H")

#–– 3. Define the flipped Sobel function
run_sobel_flipped <- function(df, subtype, treat, mediator, outcome) {
  med_mod <- lm(as.formula(paste(mediator, "~", treat)), data = df)
  out_mod <- lm(as.formula(paste(outcome,   "~", treat, "+", mediator)), data = df)
  
  a  <- coef(summary(med_mod))[treat,    "Estimate"]
  sa <- coef(summary(med_mod))[treat,    "Std. Error"]
  b  <- coef(summary(out_mod))[mediator, "Estimate"]
  sb <- coef(summary(out_mod))[mediator, "Std. Error"]
  
  z  <- a * b / sqrt(b^2 * sa^2 + a^2 * sb^2)
  p  <- 2 * (1 - pnorm(abs(z)))
  
  tibble(
    subtype, treat, mediator, outcome,
    a, b,
    sobel_z = z,
    p.value = p
  )
}

#–– 4. Prepare counts matrix & meta
counts_mat <- as.matrix(counts)   # your original counts data.frame/matrix
meta       <- meta                # your original metadata data.frame
rownames(meta) <- colnames(counts_mat)  # ensure they match

#–– 5. Main loop to run Sobel mediation per subtype
mediation_results <- map_dfr(names(sub_map), function(key) {
  sub_name <- sub_map[[key]]
  
  # Subset to this subtype
  df_sub <- meta %>%
    filter(subtype == sub_name) %>%
    mutate(sample = row_number())  # temp row index if needed
  
  samples <- rownames(df_sub)
  
  # 5a. ORA pathways for this subtype (assuming you have signif_ora list)
  pathways <- signif_ora[[key]]$ID
  
  # 5b. Build gene‐list vectors
  mediatorLists <- set_names(pathways) %>%
    map(~ filter(all_c2, gs_name == .x)$gene_symbol)
  outcomeLists  <- set_names(outcome_map[[sub_name]]) %>%
    map(~ filter(hallmark_h, gs_name == .x)$gene_symbol)
  geneSets <- c(mediatorLists, outcomeLists)
  
  # 5c. Compute ssGSEA
  expr_sub <- counts_mat[, samples, drop = FALSE]
  my_param <- ssgseaParam(
    exprData  = expr_sub,
    geneSets  = geneSets,
    normalize = TRUE
  )
  scores_mat <- gsva(my_param)
  
  # 5d. Merge scores & add binary treatment column
  df_sub2 <- df_sub %>%
    mutate(sample = rownames(df_sub)) %>%
    left_join(
      as.data.frame(t(scores_mat)) %>%
        rownames_to_column("sample"),
      by = "sample"
    ) %>%
    mutate(
      immune_active = if_else(
        TME_signature == "Immunoactive",
        1L, 0L
      )
    )
  
  # 5e. Run Sobel for every mediator → outcome pair
  expand_grid(
    mediator = names(mediatorLists),
    outcome  = names(outcomeLists)
  ) %>%
    pmap_dfr(~ run_sobel_flipped(
      df       = df_sub2,
      subtype  = sub_name,
      treat    = "immune_active",
      mediator = ..1,
      outcome  = ..2
    ))
})

# 1. Tag significant tests
sig_thresh <- 0.05
med2 <- mediation_results %>%
  mutate(signif = p.value < sig_thresh)

# 2. Count how many subtypes are significant per mediator–outcome
unique_counts <- med2 %>%
  group_by(mediator, outcome) %>%
  summarize(n_sig = sum(signif), .groups="drop") %>%
  filter(n_sig == 1)

# 3. Pull only the rows that are both significant and in that “n_sig == 1” set
unique_sig <- med2 %>%
  inner_join(unique_counts, by = c("mediator","outcome")) %>%
  filter(signif) %>%
  dplyr::select(subtype, mediator, outcome, sobel_z, p.value)

# View
unique_sig


