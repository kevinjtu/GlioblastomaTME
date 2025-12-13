#libraries####
library(data.table)
library(dplyr)
library(forcats)
library(compositions)
library(textshape)
library(survival)
library(broom)

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

#Cluster kaplan meier curves with survival####
setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Heatmap - TME Signatures")
clusters <- as.data.frame(fread("GBM TME Classifier Results.csv")) %>% dplyr::select(-V1)
colnames(clusters) <- c("sample_id", "TME_signature")
clusters <- clusters[-1,]

subtype_cols <- c("OPC_like_Cancer", "AC_like_Cancer", "NPC_like_Cancer", "MES_like_Cancer")

data <- data %>%
  mutate(
    subtype = subtype_cols[
      max.col( as.matrix( across(all_of(subtype_cols)) ),
               ties.method = "first")])

df_clusters <- data %>%
  filter(study != "IVYGAP", recurrence == "primary", IDH_status == "wt") %>%
  dplyr::select(sample_id, age, MGMT_status, IDH_status, study, 
                radiation, TMZ_treatment, resection, KarnPerfScore,
                OS_event, OS_years, DFS_event, DFS_years, 
                PFS_event, PFS_years, subtype) %>%
  left_join(clusters, by = "sample_id") %>%
  dplyr::mutate(TME_signature = factor(TME_signature), subtype = factor(subtype)) %>%
  dplyr::mutate(
    MGMT_status    = factor(MGMT_status, levels = c("hypomethylated","hypermethylated")),
    TME_signature  = factor(TME_signature, levels = c("Immunoactive","Immunosuppressed")),
    radiation    = factor(radiation, levels = c("0","1")),
    TMZ_treatment  = factor(TMZ_treatment, levels = c("0","1"))
  ) %>%
  dplyr::filter(study == "TCGA")

# Loop through each subtype and generate KM plots
# make sure factors are set
df_clusters$TME_signature <- factor(df_clusters$TME_signature,
                                    levels = c("Immunosuppressed", "Immunoactive"))
df_clusters$subtype       <- factor(df_clusters$subtype)

# 1. Overall cohort
# define survival object
surv_all <- with(df_clusters, Surv(time = OS_years, event = OS_event))

# fit KM by TME_signature
fit_all <- survfit(surv_all ~ TME_signature, data = df_clusters)

# plot
ggsurvplot(
  fit_all,
  data       = df_clusters,
  pval       = TRUE,                    # show log-rank p
  legend.title = "TME Signature",
  xlab       = "Time (years)",
  ylab       = "Overall Survival Probability",
  title      = "All Samples: OS by TME Signature",
  palette    = c("blue", "red")
)

# 2. Stratified by subtype
# loop over each subtype and plot
for(sub in levels(df_clusters$subtype)) {
  df_sub <- subset(df_clusters, subtype == sub)
  
  # only plot if there are at least two TME groups present
  if(length(unique(df_sub$TME_signature)) > 1 && 
     sum(!is.na(df_sub$OS_event) & !is.na(df_sub$OS_years)) > 0) {
    
    # define survival object
    surv_sub <- with(df_sub, Surv(time = OS_years, event = OS_event))
    
    # fit
    fit_sub <- survfit(surv_sub ~ TME_signature, data = df_sub)
    
    # plot
    print(
      ggsurvplot(
        fit_sub,
        data         = df_sub,
        pval         = TRUE,
        legend.title = "TME Signature",
        xlab         = "Time (years)",
        ylab         = "Overall Survival Probability",
        title        = paste0(sub, ": OS by TME Signature"),
        palette      = c("blue", "red")
      )$plot
    )
  }
}

#cox model for the entire cohort ####
#— 1. Read in clusters and rename 
setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Heatmap - TME Signatures")
clusters <- fread("GBM TME Classifier Results.csv") %>% 
  dplyr::select(-V1) %>% 
  setNames(c("sample_id","TME_signature")) %>% 
  slice(-1) %>% 
  mutate(TME_signature = factor(TME_signature,
                                levels = c("Immunosuppressed","Immunoactive")))

#— 2. Define subtype call (unchanged)
subtype_cols <- c("OPC_like_Cancer","AC_like_Cancer","NPC_like_Cancer","MES_like_Cancer")
data <- data %>% 
  mutate(subtype = subtype_cols[
    max.col(as.matrix(across(all_of(subtype_cols))), ties.method="first")
  ])

#— 3. Pull together all primary, IDH-wt samples across all studies
df_all <- data %>%
  filter(study != "IVYGAP",
         recurrence == "primary",
         IDH_status == "wt") %>%
  dplyr::select(sample_id, age, MGMT_status, IDH_status, study,
         radiation, TMZ_treatment, resection, KarnPerfScore,
         OS_event, OS_years, subtype) %>%
  left_join(clusters, by="sample_id") %>%
  mutate(
    study         = factor(study),
    TME_signature = factor(TME_signature,
                           levels = c("Immunosuppressed","Immunoactive")),
    MGMT_status   = factor(MGMT_status,
                           levels = c("hypomethylated","hypermethylated")),
    radiation     = factor(radiation, levels=c("0","1")),
    TMZ_treatment = factor(TMZ_treatment, levels=c("0","1"))
  )

#— 4. Fit Cox PH model for OS, stratifying by study:
surv_obj <- with(df_all, Surv(OS_years, OS_event))
cox_strata <- coxph(surv_obj ~ TME_signature + strata(study), data = df_all)

#— 5. Look at the results
summary(cox_strata)


#cox model for each subtype
# make sure subtype is a factor
df_all <- df_all %>% mutate(subtype = factor(subtype))

# get vector of subtype names
subtypes <- levels(df_all$subtype)

# loop and fit one Cox model per subtype
for(sub in subtypes) {
  df_sub <- df_all %>% filter(subtype == sub)
  
  # skip if only one TME group present
  if(length(unique(df_sub$TME_signature)) < 2) {
    message("Skipping ", sub, ": only one TME_signature present.")
    next
  }
  
  # define survival object
  surv_obj <- with(df_sub, Surv(OS_years, OS_event))
  
  # fit CoxPH: TME_signature + strata(study)
  cox_fit <- coxph(surv_obj ~ TME_signature + KarnPerfScore + strata(study), data = df_sub)
  
  # pull out the TME_signature coefficient
  smry <- summary(cox_fit)
  coef_tbl <- smry$coefficients
  # row name should be "TME_signatureImmunoactive" if your reference is Immunosuppressed
  row_name <- grep("^TME_signature", rownames(coef_tbl), value = TRUE)
  
  hr  <- exp(coef_tbl[row_name, "coef"])
  lo  <- exp(coef_tbl[row_name, "coef"] - 1.96 * coef_tbl[row_name, "se(coef)"])
  hi  <- exp(coef_tbl[row_name, "coef"] + 1.96 * coef_tbl[row_name, "se(coef)"])
  p   <- coef_tbl[row_name, "Pr(>|z|)"]
  
  cat("\nSubtype =", sub, "\n")
  cat(sprintf("  HR (Immunoactive vs. Immunosuppressed) = %.2f (95%% CI %.2f–%.2f), p = %.3g\n",
              hr, lo, hi, p))
}