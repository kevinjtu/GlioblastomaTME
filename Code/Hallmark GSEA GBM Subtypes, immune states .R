#libraries####
library(data.table)
library(dplyr)
library(tibble)
library(limma)
library(msigdbr)
library(clusterProfiler)
library(ggplot2)
library(ggpubr)

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

#Differential gene analysis####
# join metadata
meta <- tumor_subtype %>%
  inner_join(TME_sig, by = "patient_id") %>%      # 1) merge
  distinct(patient_id, .keep_all = TRUE) %>%       # 2) keep one row per patient_id
  column_to_rownames(var = "patient_id")  

# make sure samples line up
common_samples <- intersect(colnames(TEMPUS), rownames(meta))
counts <- TEMPUS[, common_samples]
meta   <- meta[ common_samples, , drop=FALSE ]

#deseq2 analysis
subtypes <- unique(meta$subtype)
deg_results <- list()
sig_genes   <- list()

# Rename TPM matrix so it doesn't shadow a function
expr_mat <- as.matrix(TEMPUS)

for(sub in subtypes) {
  # subset to this subtype
  keep      <- meta$subtype == sub
  expr_sub  <- expr_mat[, keep, drop = FALSE]
  md_sub    <- meta[keep, , drop = FALSE]
  
  # log2 TPM + pseudocount
  logTPM <- log2(expr_sub + 1)
  
  # filter low-expression (TPM ≥1 in ≥3 samples)
  keep_g <- rowSums(expr_sub >= 1) >= 3
  logTPM <- logTPM[keep_g, ]
  
  # set up design
  md_sub$TME_signature <- factor(
    md_sub$TME_signature,
    levels = c("Immunosuppressed","Immunoactive")
  )
  design <- model.matrix(~0 + TME_signature, data = md_sub)
  colnames(design) <- levels(md_sub$TME_signature)
  
  # contrast
  contrast <- makeContrasts(
    Immunoactive_vs_Immunosuppressed = Immunoactive - Immunosuppressed,
    levels = design
  )
  
  # fit limma
  fit  <- lmFit(logTPM, design)
  fitC <- contrasts.fit(fit, contrast)
  fitE <- eBayes(fitC, trend = TRUE)
  
  # extract & sort
  topTab <- topTable(fitE,
                     coef   = "Immunoactive_vs_Immunosuppressed",
                     number = Inf) %>%
    rownames_to_column("gene") %>%
    arrange(adj.P.Val)
  
  deg_results[[sub]] <- topTab
  
  # significant genes
  sig_genes[[sub]] <- list(
    up   = topTab %>% filter(adj.P.Val < 0.05, logFC > 1)  %>% pull(gene),
    down = topTab %>% filter(adj.P.Val < 0.05, logFC < -1) %>% pull(gene)
  )
}

#GSEA analysis####
#Build GSEA gene list
genes_for_GSEA <- unique(unlist(lapply(sig_genes, `[[`, "up")))

hallmark_df <- msigdbr(
  species  = "Homo sapiens",
  category = "H"         # “H” = Hallmark
) %>%
  dplyr::select(gs_name, gene_symbol) %>%
  dplyr::rename(term = gs_name, gene = gene_symbol)

#Run preranked GSEA for each subtype
gsea_results <- list()

set.seed(123)  # for reproducibility of the jitter

gsea_results <- list()

for(sub in names(deg_results)) {
  df <- deg_results[[sub]]
  
  # 1) build ranked vector and break ties by jitter
  rank_vec <- df$logFC
  names(rank_vec) <- df$gene
  # jitter only enough to separate exact ties
  rank_vec <- jitter(rank_vec, factor = 1e-8)
  rank_vec <- sort(rank_vec, decreasing = TRUE)
  
  # 2) run GSEA with eps = 0
  gsea_results[[sub]] <- GSEA(
    geneList     = rank_vec,
    TERM2GENE    = hallmark_df,
    pvalueCutoff = 0.05,
    eps          = 0,        # allow very small p-values
    verbose      = FALSE,
    seed         = TRUE      # fix random seed inside fgsea
  )
}


dotplot(gsea_results[["MES_like_Cancer"]], showCategory = 10) +
  ggtitle("Hallmark GSEA: MES-like Immunoactive vs Suppressed")

dotplot(gsea_results[["OPC_like_Cancer"]], showCategory = 10) +
  ggtitle("Hallmark GSEA: OPC-like Immunoactive vs Suppressed")

dotplot(gsea_results[["NPC_like_Cancer"]], showCategory = 10) +
  ggtitle("Hallmark GSEA: NPC-like Immunoactive vs Suppressed")

dotplot(gsea_results[["AC_like_Cancer"]], showCategory = 10) +
  ggtitle("Hallmark GSEA: AC-like Immunoactive vs Suppressed")

#STING Expression between subtypes####
expr_mat <- expr_mat[, rownames(meta), drop = FALSE]

sting_df <- data.frame(
  sample        = colnames(expr_mat),
  STING1        = expr_mat["STING1", ],              # raw TPM (or log2(TPM+1))
  subtype       = meta[colnames(expr_mat), "subtype"],
  TME_signature = meta[colnames(expr_mat), "TME_signature"]
)

# define comparisons: each vs MES
comparisons <- list(
  c("MES_like_Cancer", "OPC_like_Cancer"),
  c("MES_like_Cancer", "NPC_like_Cancer"),
  c("MES_like_Cancer", "AC_like_Cancer")
)

sting_df$TME_signature <- factor(
  sting_df$TME_signature,
  levels = c("Immunosuppressed", "Immunoactive")
)

ggplot(sting_df, aes(x = subtype, y = STING1, fill = TME_signature)) +
  # bars = group means
  stat_summary(
    fun   = mean,
    geom  = "bar",
    color = "black",
    position = position_dodge(width = 0.8),
    width    = 0.7
  ) +
  
  # error bars = mean ± SE
  stat_summary(
    fun.data  = mean_se,
    geom      = "errorbar",
    position  = position_dodge(width = 0.8),
    width     = 0.2
  ) +
  
  # individual points
  geom_jitter(
    aes(color = TME_signature),
    position = position_jitterdodge(
      jitter.width = 0.15,
      dodge.width  = 0.8
    ),
    alpha = 0.6,
    size  = 1.5,
    show.legend = FALSE
  ) +
  
  # pairwise t‐tests within each subtype
  stat_compare_means(
    aes(group = TME_signature),
    method     = "t.test",
    label      = "p.signif",
    comparisons = list(
      c("Immunosuppressed", "Immunoactive")
    ),
    position    = position_dodge(width = 0.8),
    tip.length  = 0.02
  ) +
  
  # overall ANOVA across all groups (optional)
  stat_compare_means(
    method = "anova",
    label  = "p.format",
    label.y = max(sting_df$STING1) * 1.1
  ) +
  
  scale_fill_brewer(palette = "Set1") +
  theme_bw(base_size = 14) +
  labs(
    x     = "GBM Subtype",
    y     = "STING1 expression (TPM)",
    fill  = "TME Signature",
    title = "STING1 Expression by Subtype and Immune State"
  )
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
library(dplyr)
library(tidyr)

wide_p <- res_df %>%
  select(gene, subtype, adj.P.Val) %>%
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