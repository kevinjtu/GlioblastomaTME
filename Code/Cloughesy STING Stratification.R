#libraries####
library(InstaPrism)
library(data.table)
library(Biobase)
library(dplyr)
library(textshape)
library(tibble)
library(readxl)
library(biomaRt)
library("org.Hs.eg.db")
library(openxlsx)
library(AnnotationDbi)
library(limma)
library(msigdbr)
library(clusterProfiler)
library(GSVA)
library(stringr)
library(survival)
library(survminer)

setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Expression Data")

#STING Pathway Enrichment####
# Get Reactome pathways from msigdbr
sting_pathway <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME") %>%
  filter(gs_name == "REACTOME_STING_MEDIATED_INDUCTION_OF_HOST_IMMUNE_RESPONSES")

sting_genes <- unique(sting_pathway$gene_symbol)

# Define niche gene sets
niche_gsea <- list(
  niche_all = c(  # For overall enrichment
    "HILPDA", "LOX", "VEGFA", "CAV1", "CXCL8", "VIM", "TNC",
    "BCAN", "OLIG1", "KDR", "FLT1", "TMEM173", "PECAM1"
  ),
  niche_up = c(   # For directionality
    "HILPDA", "LOX", "VEGFA", "CAV1", "CXCL8", "VIM", "TNC"
  ),
  niche_down = c( # For directionality
    "BCAN", "OLIG1", "KDR", "FLT1", "TMEM173", "PECAM1"
  )
)

# Combine all gene sets into one list
gene_set_list <- c(
  list(REACTOME_STING = sting_genes),
  niche_gsea
)

# load the dataset, neo1
setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Expression Data")
neoadjuvant1 <- read_excel("NeoadjuvantStage1_counts.xlsx") %>%
  distinct(Genes, .keep_all = TRUE) %>%
  column_to_rownames(var = "Genes")
  
#get into correct data type
bulk <- as.data.frame(neoadjuvant1)
bulk_ma <- as.matrix(bulk)
if ("STING1" %in% rownames(bulk_ma)) {
  rownames(bulk_ma)[rownames(bulk_ma) == "STING1"] <- "TMEM173"}

params <- ssgseaParam(
  exprData = bulk_ma,
  geneSets = gene_set_list,
  norm     = TRUE       # normalize the scores
  )

ssgsea_eset_neo1 <- as.data.frame(t(gsva(
  params,
  verbose = FALSE
)))
ssgsea_eset_neo1$patient_id <- rownames(ssgsea_eset_neo1)

# load the dataset, neo2
setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Expression Data")
neoadjuvant2 <- read_excel("NeoadjuvantStage2_counts.xlsx") %>%
  distinct(Gene, .keep_all = TRUE) %>%
  column_to_rownames(var = "Gene")

#get into correct data type
bulk <- as.data.frame(neoadjuvant2)
bulk_ma <- as.matrix(bulk)
if ("STING1" %in% rownames(bulk_ma)) {
  rownames(bulk_ma)[rownames(bulk_ma) == "STING1"] <- "TMEM173"}

params <- ssgseaParam(
  exprData = bulk_ma,
  geneSets = gene_set_list,
  norm     = TRUE       # normalize the scores
)

ssgsea_eset_neo2 <- as.data.frame(t(gsva(
  params,
  verbose = FALSE
)))
ssgsea_eset_neo2$patient_id <- rownames(ssgsea_eset_neo2)

sting_scores <- rbind(ssgsea_eset_neo1, ssgsea_eset_neo2) %>%
  #dplyr::select(patient_id, REACTOME_STING) %>%
  mutate(patient_id = str_remove_all(patient_id, "Pt|_[A-Z]")) %>%
  mutate(patient_id = str_extract(as.character(patient_id), "^\\d+"))

#Mes-like classification####
subtype_cols <- c("OPC_like", "AC_like", "NPC_like", "MES_like")
setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Instaprism Results")
neo1_insta <- as.data.frame(read_excel("NeoadjuvantStage1 RNAseq InstaPrism results.xlsx", sheet = 4)) %>%
  dplyr::rename(patient_id = `...1`)
neo2_insta <- as.data.frame(read_excel("NeoadjuvantStage2 RNAseq InstaPrism results.xlsx", sheet = 4)) %>%
  dplyr::rename(patient_id = `...1`)
combined_dataset <- rbind(neo1_insta, neo2_insta)

new_names <- names(combined_dataset)
new_names <- str_replace_all(new_names, "\\-", "")       # Remove hyphens
new_names <- str_replace_all(new_names, "\\/", "_")       # Remove slashes
new_names <- str_replace_all(new_names, " ", "_")         # Replace spaces with underscores
new_names <- str_replace_all(new_names, "\\+", "")         # Remove plus signs
combined_dataset <- setNames(combined_dataset, new_names)

opc_like_columns <- grep("^OPClike", names(combined_dataset), value = TRUE)
ac_like_columns  <- grep("^AClike", names(combined_dataset), value = TRUE)
npc_like_columns <- grep("^NPClike", names(combined_dataset), value = TRUE)
mes_like_columns <- grep("^MESlike", names(combined_dataset), value = TRUE)

combined_dataset$OPC_like <- rowSums(combined_dataset[opc_like_columns], na.rm = TRUE)
combined_dataset$AC_like <- rowSums(combined_dataset[ac_like_columns], na.rm = TRUE)
combined_dataset$NPC_like <- rowSums(combined_dataset[npc_like_columns], na.rm = TRUE)
combined_dataset$MES_like <- rowSums(combined_dataset[mes_like_columns], na.rm = TRUE)
combined_dataset <- combined_dataset[, !(names(combined_dataset) %in% c(opc_like_columns, ac_like_columns, npc_like_columns, mes_like_columns))]

combined_dataset <- combined_dataset %>%
  mutate(subtype = subtype_cols[
    max.col( as.matrix( across(all_of(subtype_cols)) ),
             ties.method = "first")]) %>%
  mutate(MES_like_adj = MES_like / (MES_like + AC_like + NPC_like + OPC_like))

neo_subtype <- combined_dataset %>%
  dplyr::select(patient_id, subtype, MES_like, MES_like_adj) %>%
  mutate(patient_id = str_remove_all(patient_id, "Pt|_[A-Z]")) %>%
  mutate(patient_id = str_extract(as.character(patient_id), "^\\d+"))


#clinical data####
setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Clinical Data")
clinical <- as.data.frame(read_excel("Cloughesy Clinical.xlsx")) %>%
  dplyr::rename(patient_id = `Subject ID`, PFS_years = "Progression Free Survival (Days)",
         PFS_event = "PFS Censored Status (0 = PF)", OS_event = "OS Censored Status (0 = Alive, 1 = Dead)",
         OS_years = "Overall Survival (Days)", age = Age, mgmt_status = MGMT, 
         idh_status = IDH, KarnPerfScore = "KPS at registration", batch = Stage) %>%
  dplyr::mutate(OS_years = as.numeric(OS_years) / 365.25) %>%
  dplyr::mutate(PFS_years = as.numeric(PFS_years) / 365.25) %>%
  #mutate(patient_id = as.character(patient_id)) %>%
  merge(neo_subtype, by = "patient_id") %>%
  merge(sting_scores, by = "patient_id") %>%
  mutate(
    mgmt_status = case_when(
      tolower(mgmt_status) == "methylated"     ~ "Methylated",
      tolower(mgmt_status) == "unmethylated"   ~ "Unmethylated",
      TRUE                                     ~ "Unknown"),
    idh_status = case_when(
      tolower(idh_status) == "positive"        ~ "Positive",
      tolower(idh_status) == "negative"        ~ "Negative",
      TRUE                                     ~ "Unknown"))
  
#KM Curve####
# Make sure Treatment is a factor
clinical$Treatment <- as.factor(clinical$Treatment)
clinical_sub <- clinical %>%
  #filter(Treatment == "A") %>%
  #filter(batch == "1") %>%
  mutate(niche_all_tertile = ntile(niche_all, 2))

# Create survival object
surv_obj <- Surv(time = clinical_sub$OS_years, event = clinical_sub$OS_event)

# Fit survival model
fit <- survfit(surv_obj ~ niche_all_tertile, data = clinical_sub)

# Plot Kaplan-Meier curve with log-rank p-value
ggsurvplot(
  fit,
  data = clinical_sub,
  pval = TRUE,                       # Show log-rank p-value
  #conf.int = TRUE,                  # Show confidence interval
  risk.table = TRUE,                # Add risk table
  xlab = "Time (Years)", 
  ylab = "Overall Survival Probability",
  legend.title = "Treatment",
)

#cox model
cox <- coxph(Surv(OS_years, OS_event) ~ niche_all + KarnPerfScore + age + Treatment, data = clinical_sub)
summary(cox)

#Cox model####
clinical_idh_neg <- subset(clinical, idh_status == "Negative")
#clinical_idh_neg <- subset(clinical_idh_neg, `Treatment Group` == "A")

# Optional: ensure data types are correct
#clinical_idh_neg$Treatment <- as.factor(clinical_idh_neg$Treatment)
clinical_idh_neg$subtype <- as.factor(clinical_idh_neg$subtype)
clinical_idh_neg$mgmt_status <- as.factor(clinical_idh_neg$subtype)

# Split into MES_like and non-MES_like groups
mes_data <- subset(clinical_idh_neg, subtype == "MES_like")
nonmes_data <- subset(clinical_idh_neg, subtype != "MES_like")

# Cox model for MES_like tumors
cox_mes <- coxph(Surv(OS_years, OS_event) ~ REACTOME_STING + age + KarnPerfScore + `Treatment Group`, data = mes_data)
summary(cox_mes)

# Cox model for non-MES_like tumors
cox_nonmes <- coxph(Surv(OS_years, OS_event) ~ REACTOME_STING + age + KarnPerfScore + `Treatment Group` + batch, data = nonmes_data)
summary(cox_nonmes)






clinical$MES_like_scaled <- scale(log10(clinical$MES_like_adj + 1e-8))
clinical$subtype_group <- ifelse(clinical$subtype == "MES_like", "MES", "non-MES")

clinical_idh_neg <- subset(clinical, idh_status == "Negative")
clinical_idh_neg <- subset(clinical_idh_neg, `Treatment Group` == "A")

#for all data, strata by subtype group
summary(coxph(Surv(OS_years, OS_event) ~ REACTOME_STING + subtype_group + KarnPerfScore + age + mgmt_status, data = clinical_idh_neg))

#only non-mes
clinical_idh_neg <- subset(clinical_idh_neg, subtype_group == "non-MES")
summary(coxph(Surv(OS_years, OS_event) ~ REACTOME_STING + KarnPerfScore + age + mgmt_status, data = clinical_idh_neg))


# Re-run the model with scaled MES_like
cox_model_scaled <- coxph(
  Surv(OS_years, OS_event) ~ REACTOME_STING * subtype_group + KarnPerfScore + age + mgmt_status,
  data = clinical_idh_neg
)

summary(cox_model_scaled)