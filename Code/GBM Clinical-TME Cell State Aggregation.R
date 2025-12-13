#libraries####
library(data.table)
library(Biobase)
library(dplyr)
library(textshape)
library(tibble)
library(readxl)
library(stringr)
library(biomaRt)
library("org.Hs.eg.db")

sting_aliases <- c("STING1", "ERIS", "FLJ38577", "MITA", "MPYS", "NET23", "TMEM173")
arid1a_aliases <- c("ARID1A", "BAF250", "B120", "BAF250a", "SMARCF1", "C1orf4", "P270")

fpkm_to_tpm <- function(fpkm_mat) {
  scaling_factors <- colSums(fpkm_mat, na.rm = TRUE) / 1e6
  tpm <- sweep(fpkm_mat, 2, scaling_factors, "/")
  return(tpm)
}
#notes: IDH_mutation status refers to idh1 mutation status
#notes: check IDH1 mutation status in tempus and tcga - Kevin 5/13/2025

#TCGA cleaning####
setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Clinical Data")
tcga_patient <- fread("tcga_clinical_patient.txt")
tcga_patient <- tcga_patient[-(1:3), ]
colnames(tcga_patient) <- as.character(tcga_patient[1, ])
tcga_patient <- tcga_patient[-1, ] 
tcga_sample <- fread("tcga_clinical_sample.txt")
tcga_sample <- tcga_sample[-(1:3), ]
colnames(tcga_sample) <- as.character(tcga_sample[1, ])
tcga_sample <- tcga_sample[-1, ] 
tcga_subtype <- read_xlsx("TCGA Subtype.xlsx")
tcga_subtype <- tcga_subtype[,(1:2)]
tcga_subtype <- tcga_subtype %>% 
  mutate(sampleId = gsub("\\.", "-", sampleId)) %>%
  dplyr::rename(sample_id = sampleId) %>%
  dplyr::rename(Subtype = Group)

tcga_meth <- fread("TCGA_MGMT_Methylation_HM27.txt")
tcga_meth <- tcga_meth %>%
  mutate(MGMT_status = ifelse(`MGMT: Methylation (HM27)` > 0.30, "hypermethylated", "hypomethylated")) %>%
  dplyr::select("Sample ID", MGMT_status) %>%
  dplyr::rename(SAMPLE_ID = "Sample ID")

tcga_idh <- fread("idhmut_gbm_tcga_clinical_data.tsv", check.names = TRUE)
tcga_mut <- fread("gbm_mutation_tcga.tsv", check.names = TRUE)
tcga_idh <- tcga_idh %>% mutate(IDH_status = "IDHmut") %>% dplyr::select(Sample.ID, IDH_status)
tcga_mut <- tcga_mut %>%
  merge(tcga_idh, by = "Sample.ID", all.x = TRUE) %>%
  mutate(
    IDH_status = ifelse(is.na(IDH_status), "IDHwt", IDH_status)
  ) %>%
  dplyr::rename(SAMPLE_ID = "Sample.ID") %>%
  dplyr::select(SAMPLE_ID, IDH_status)

tcga <- tcga_sample %>%
  merge(tcga_patient, by = "PATIENT_ID", all.x = TRUE) %>%
  mutate(OS_MONTHS = ifelse(OS_MONTHS == "[Not Available]", NA, OS_MONTHS),
         DFS_MONTHS = ifelse(DFS_MONTHS == "[Not Available]", NA, DFS_MONTHS)) %>%
  mutate(
    OS_MONTHS = as.numeric(OS_MONTHS),
    DFS_MONTHS = as.numeric(DFS_MONTHS)
  ) %>%
  mutate(OS_years = round(OS_MONTHS / 12, 4),
         DFS_years = round(DFS_MONTHS / 12, 4)) %>%
  merge(tcga_meth, by = "SAMPLE_ID", all.x = TRUE) %>%
  merge(tcga_mut, by = "SAMPLE_ID", all.x = TRUE) %>%
  dplyr::select(PATIENT_ID, SAMPLE_ID, AGE, SEX, RACE, OS_STATUS, OS_years, 
                DFS_STATUS, DFS_years, SAMPLE_TYPE, KARNOFSKY_PERFORMANCE_SCORE, 
                MGMT_status, IDH_status, PHARMACEUTICAL_TX_ADJUVANT, RADIATION_TREATMENT_ADJUVANT) %>%
  dplyr::rename(patient_id = PATIENT_ID, sample_id = SAMPLE_ID, age = AGE, sex = SEX, race = RACE,
                OS_event = OS_STATUS, DFS_event = DFS_STATUS, recurrence = SAMPLE_TYPE,
                KarnPerfScore = KARNOFSKY_PERFORMANCE_SCORE, chemo = PHARMACEUTICAL_TX_ADJUVANT, radiation = RADIATION_TREATMENT_ADJUVANT)%>%
  mutate(study = "TCGA", location = NA, PFS_event = NA, PFS_years = NA,
         resection = NA, chemo = NA, radiation, STING_expr = NA, 
         ARID1A_expr = NA, grade = NA, platform = "rnaseq") %>%
  dplyr::select(patient_id, sample_id, age, sex, race, grade,  
                OS_event, OS_years, DFS_event, DFS_years, PFS_event, PFS_years,
                recurrence, location, resection, chemo, radiation, KarnPerfScore, 
                MGMT_status, IDH_status, study, platform)

setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Instaprism Results")
tcga_immune <- read_xlsx("TCGA RNAseq InstaPrism results.xlsx", sheet = 4)
colnames(tcga_immune)[1] <- "sample_id"

setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Sample GSEA")
tcga_gsea <- read_xlsx("TCGA_GSEA.xlsx", sheet = 1)

setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Expression Data")
tcga_expression <- fread("TCGA RNAseq expression.txt")
tcga_expression <- tcga_expression %>%
  filter(!is.na(Hugo_Symbol)) %>%
  distinct(Hugo_Symbol, .keep_all = TRUE) %>%
  dplyr::select(-Entrez_Gene_Id)
tcga_expression <- column_to_rownames(tcga_expression, var = "Hugo_Symbol")
tcga_tpm <- apply(tcga_expression, 2, function(x) {
  if (sum(x, na.rm = TRUE) == 0) {
    return(rep(0, length(x)))
  } else {
    return((x / sum(x, na.rm = TRUE)) * 1e6)
  }
})
tcga_tpm <- as.data.frame(t(tcga_tpm))

gene_names_upper <- toupper(colnames(tcga_tpm))
sting_match <- gene_names_upper %in% toupper(sting_aliases)
arid1a_match <- gene_names_upper %in% toupper(arid1a_aliases)
tcga_genes <- tcga_tpm[, sting_match | arid1a_match, drop = FALSE]
tcga_genes$sample_id <- as.character(rownames(tcga_tpm))

tcga_merged <- tcga %>%
  merge(tcga_subtype, by = "sample_id", all.x = TRUE) %>%
  merge(tcga_genes, by = "sample_id", all.x = TRUE) %>%
  merge(tcga_gsea, by = "sample_id", all.x = TRUE) %>%
  merge(tcga_immune, by = "sample_id", all.x = TRUE) %>%
  filter(!is.na(Oligodendrocyte))

tcga_merged <- tcga_merged %>%
  dplyr::rename(STING1_tpm = TMEM173, ARID1A_tpm = ARID1A) %>%
  dplyr::select(-Subtype) %>%
  dplyr::select(
    patient_id,
    sample_id,
    everything()
  )

#TEMPUS cleaning####
setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Clinical Data/TEMPUS")
tempus_patient <- fread("patient_summary[56].csv")
tempus_patient <- tempus_patient %>%
  dplyr::select(patient_id, gender, birth_time_from_index, race_concept_canonical_name) %>%
  dplyr::rename(age = birth_time_from_index, sex = gender, race = race_concept_canonical_name) %>%
  mutate(age = as.numeric(age) / -365.25) %>%
  mutate(sample_id = NA,OS_event = NA, OS_years = NA, DFS_event = NA, DFS_years = NA, 
         recurrence = "primary", PFS_event = NA, PFS_years = NA, IDH_status = NA,
         grade = NA, location = NA, resection = NA, chemo = NA, radiation = NA, KarnPerfScore = NA,
         MGMT_status = NA, chemo = NA, STING_expr = NA, 
         ARID1A_expr = NA, study = "TEMPUS", platform = "rnaseq") %>%
  dplyr::select(patient_id, sample_id, age, sex, race, grade,  
                OS_event, OS_years, DFS_event, DFS_years, PFS_event, PFS_years,
                recurrence, location, resection, chemo, radiation, KarnPerfScore, 
                MGMT_status, IDH_status, study, platform)

setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Instaprism Results")
tempus_immune <- readxl::read_xlsx("TEMPUS RNAseq InstaPrism results.xlsx", sheet = 4)

colnames(tempus_immune)[1] <- "patient_id"
tempus_immune$patient_id <- sub("^X", "", tempus_immune$patient_id)
tempus_patient$patient_id <- as.character(tempus_patient$patient_id)
tempus_immune$patient_id <- as.character(tempus_immune$patient_id)

setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Sample GSEA")
tempus_gsea <- readxl::read_xlsx("TEMPUS_GSEA.xlsx", sheet = 1)

setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Expression Data")
tempus_expression <- fread("TEMPUS RNAseq expression.csv")
tempus_expression <- as.data.frame(tempus_expression)
rownames(tempus_expression) <- tempus_expression[[1]]
tempus_expression <- tempus_expression[, -1]
tempus_expression <- t(tempus_expression)
tempus_expression <- as.data.frame(tempus_expression)

gene_names_upper <- toupper(colnames(tempus_expression))
sting_match <- gene_names_upper %in% toupper(sting_aliases)
arid1a_match <- gene_names_upper %in% toupper(arid1a_aliases)
tempus_genes <- tempus_expression[, sting_match | arid1a_match, drop = FALSE]
tempus_genes$patient_id <- as.character(tempus_expression$name)
tempus_genes <- tempus_genes[, c("patient_id", setdiff(colnames(tempus_genes), "patient_id"))]

tempus_merged <- tempus_patient %>%
  #merge(tempus_subtype, by = "sample_id", all.x = TRUE) %>%
  merge(tempus_genes, by = "patient_id", all.x = TRUE) %>%
  merge(tempus_gsea, by = "patient_id", all.x = TRUE) %>%
  merge(tempus_immune, by = "patient_id", all.x = TRUE) %>%
  filter(!is.na(Oligodendrocyte))

tempus_merged <- tempus_merged %>%
  dplyr::rename(STING1_tpm = STING1, ARID1A_tpm = ARID1A)

#Wu2020 cleaning####
setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Clinical Data")
wu_patient <- read_excel("Wu2020 Clinical.xlsx")
wu_patient <- wu_patient %>%
  dplyr::select(`Sample ID`, `WHO grade`, `survival status`, `OS (months)`, progression, 
                `PFS  (months)`, `IDH1 mutaton status`, `IDH2 mutation status`) %>%
  dplyr::rename(patient_id = `Sample ID`, grade = `WHO grade`, OS_event = `survival status`, OS_years = `OS (months)`,,
                PFS_event = progression, PFS_years = `PFS  (months)`, IDH_status = `IDH1 mutaton status`) %>%
  mutate(OS_years = as.numeric(OS_years) / 12) %>%
  mutate(PFS_years = as.numeric(PFS_years) / 12) %>%
  mutate(sample_id = patient_id, age = NA, race = NA, sex = NA, DFS_event = NA, DFS_years = NA, recurrence = "primary", 
         location = NA, STING_expr = NA, ARID1A_expr = NA, study = "Wu2020", resection = NA,
         KarnPerfScore = NA, MGMT_status = NA, chemo = NA, radiation = NA, platform = "rnaseq") %>%
  dplyr::select(patient_id, sample_id, age, sex, race, grade,  
                OS_event, OS_years, DFS_event, DFS_years, PFS_event, PFS_years,
                recurrence, location, resection, chemo, radiation, KarnPerfScore, 
                MGMT_status, IDH_status, study, platform)

setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Expression Data")
wu_expression <- fread("Wu2020 TPM Expression.txt")
wu_expression[, ensg_id := sub("\\.\\d+$", "", V1)]  # Remove version numbers
hugo_symbols <- mapIds(org.Hs.eg.db,
                       keys = wu_expression$ensg_id,
                       column = "SYMBOL",
                       keytype = "ENSEMBL",
                       multiVals = "first")
wu_expression[, V1 := hugo_symbols][, ensg_id := NULL]  # Update V1 and clean up
wu_expression <- wu_expression[!is.na(V1)]
wu_expression <- as.data.frame(t(wu_expression))
colnames(wu_expression) <- wu_expression[1, ]
wu_expression <- wu_expression[-1, ]

gene_names_upper <- toupper(colnames(wu_expression))
sting_match <- gene_names_upper %in% toupper(sting_aliases)
arid1a_match <- gene_names_upper %in% toupper(arid1a_aliases)
wu_genes <- wu_expression[, sting_match | arid1a_match, drop = FALSE]
wu_genes$patient_id <- as.character(rownames(wu_expression))

setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Instaprism Results")
wu_immune <- read_xlsx("Wu2020 RNAseq InstaPrism results.xlsx", sheet = 4)
colnames(wu_immune)[1] <- "patient_id"

setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Sample GSEA")
wu_gsea <- read_xlsx("Wu2020_GSEA.xlsx", sheet = 1)

wu_merged <- wu_patient %>%
  #merge(tcga_subtype, by = "sample_id", all.x = TRUE) %>%
  merge(wu_genes, by = "patient_id", all.x = TRUE) %>%
  merge(wu_gsea, by = "patient_id", all.x = TRUE) %>%
  merge(wu_immune, by = "patient_id", all.x = TRUE) %>%
  filter(!is.na(Oligodendrocyte))

wu_merged <- wu_merged %>%
  dplyr::rename(STING1_tpm = STING1, ARID1A_tpm = ARID1A)


#CGGA Batch 2 cleaning####
setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Clinical Data")
cgga_patient <- fread("CGGA Clinical.txt") #CGGA.mRNAseq_693_clinical.20200506.txt
cgga_patient <- cgga_patient %>%
  filter(grepl("GBM", Histology, ignore.case = FALSE)) %>%
  dplyr::select(CGGA_ID,PRS_type, Grade, Gender, Age, OS, `Censor (alive=0; dead=1)`, IDH_mutation_status, 
                `Radio_status (treated=1;un-treated=0)`, `Chemo_status (TMZ treated=1;un-treated=0)`, MGMTp_methylation_status) %>%
  dplyr::rename(patient_id = CGGA_ID, grade = Grade, sex = Gender, age = Age, grade = Grade, OS_event = `Censor (alive=0; dead=1)`, OS_years = OS,
                IDH_status = IDH_mutation_status, 
                chemo = `Chemo_status (TMZ treated=1;un-treated=0)`, 
                radiation =  `Radio_status (treated=1;un-treated=0)`,
                MGMT_status = MGMTp_methylation_status) %>%
  mutate(OS_years = as.numeric(OS_years) / 365.25) %>%
  mutate(sample_id = patient_id, race = "Chinese", DFS_event = NA, DFS_years = NA, 
         PFS_event = NA, PFS_years = NA, recurrence = PRS_type,
         location = NA, STING_expr = NA, ARID1A_expr = NA, resection = NA, 
         KarnPerfScore = NA, study = "CGGA Batch 2", platform = "rnaseq") %>%
  dplyr::select(patient_id, sample_id, age, sex, race, grade,  
                OS_event, OS_years, DFS_event, DFS_years, PFS_event, PFS_years,
                recurrence, location, resection, chemo, radiation, KarnPerfScore, 
                MGMT_status, IDH_status, study, platform)

setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Expression Data")
cgga_expression <- fread("CGGA FPKM Expression.txt")
gene_names <- cgga_expression$Gene_Name
cgga_fpkm <- as.matrix(cgga_expression[, -1, with = FALSE])
cgga_tpm <- fpkm_to_tpm(cgga_fpkm)
cgga_tpm <- as.data.table(cgga_tpm)
cgga_tpm[, Gene_Name := gene_names]
setcolorder(cgga_tpm, "Gene_Name")
cgga_tpm <- t(cgga_tpm)
colnames(cgga_tpm) <- cgga_tpm[1, ]
cgga_tpm <- cgga_tpm[-1, ]

gene_names_upper <- toupper(colnames(cgga_tpm))
sting_match <- gene_names_upper %in% toupper(sting_aliases)
arid1a_match <- gene_names_upper %in% toupper(arid1a_aliases)
cgga_genes <- cgga_tpm[, sting_match | arid1a_match, drop = FALSE]
cgga_genes <- as.data.frame(cgga_genes)
cgga_genes$patient_id <- as.character(rownames(cgga_tpm))

setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Instaprism Results")
cgga_immune <- read_xlsx("CGGA RNAseq InstaPrism results.xlsx", sheet = 4)
colnames(cgga_immune)[1] <- "patient_id"

setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Sample GSEA")
cgga_gsea <- read_xlsx("cgga_batch1_GSEA.xlsx", sheet = 1)

cgga_merged <- cgga_patient %>%
  #merge(tcga_subtype, by = "sample_id", all.x = TRUE) %>%
  merge(cgga_genes, by = "patient_id", all.x = TRUE) %>%
  merge(cgga_gsea, by = "patient_id", all.x = TRUE) %>%
  merge(cgga_immune, by = "patient_id", all.x = TRUE) %>%
  filter(!is.na(Oligodendrocyte))

cgga_merged <- cgga_merged %>%
  dplyr::rename(STING1_tpm = TMEM173, ARID1A_tpm = ARID1A)

#CGGA Batch 1 cleaning####
setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Clinical Data")
cgga_patient_b <- fread("CGGA batch1 clinical.txt") #CGGA.mRNAseq_325_clinical.20200506.txt
cgga_patient_b <- cgga_patient_b %>%
  filter(grepl("GBM", Histology, ignore.case = FALSE)) %>%
  dplyr::select(CGGA_ID,PRS_type, Grade, Gender, Age, OS, `Censor (alive=0; dead=1)`, IDH_mutation_status, 
                `Radio_status (treated=1;un-treated=0)`, `Chemo_status (TMZ treated=1;un-treated=0)`, MGMTp_methylation_status) %>%
  dplyr::rename(patient_id = CGGA_ID, grade = Grade, sex = Gender, age = Age, grade = Grade, OS_event = `Censor (alive=0; dead=1)`, OS_years = OS,
                IDH_status = IDH_mutation_status, 
                chemo = `Chemo_status (TMZ treated=1;un-treated=0)`, 
                radiation =  `Radio_status (treated=1;un-treated=0)`,
                MGMT_status = MGMTp_methylation_status) %>%
  mutate(OS_years = as.numeric(OS_years) / 365.25) %>%
  mutate(sample_id = patient_id, race = "Chinese", DFS_event = NA, DFS_years = NA, 
         PFS_event = NA, PFS_years = NA, recurrence = PRS_type,
         location = NA, STING_expr = NA, ARID1A_expr = NA, resection = NA, 
         KarnPerfScore = NA, study = "CGGA Batch 1", platform = "rnaseq") %>%
  dplyr::select(patient_id, sample_id, age, sex, race, grade,  
                OS_event, OS_years, DFS_event, DFS_years, PFS_event, PFS_years,
                recurrence, location, resection, chemo, radiation, KarnPerfScore, 
                MGMT_status, IDH_status, study, platform)

setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Expression Data")
cgga_expression <- fread("CGGA batch1 fpkm.txt")
gene_names <- cgga_expression$Gene_Name
cgga_fpkm <- as.matrix(cgga_expression[, -1, with = FALSE])
cgga_tpm <- fpkm_to_tpm(cgga_fpkm)
cgga_tpm <- as.data.table(cgga_tpm)
cgga_tpm[, Gene_Name := gene_names]
setcolorder(cgga_tpm, "Gene_Name")
cgga_tpm <- t(cgga_tpm)
colnames(cgga_tpm) <- cgga_tpm[1, ]
cgga_tpm <- cgga_tpm[-1, ]

gene_names_upper <- toupper(colnames(cgga_tpm))
sting_match <- gene_names_upper %in% toupper(sting_aliases)
arid1a_match <- gene_names_upper %in% toupper(arid1a_aliases)
cgga_b_genes <- cgga_tpm[, sting_match | arid1a_match, drop = FALSE]
cgga_b_genes <- as.data.frame(cgga_b_genes)
cgga_b_genes$patient_id <- as.character(rownames(cgga_tpm))

setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Instaprism Results")
cgga_b_immune <- read_xlsx("CGGA Batch 1 RNAseq InstaPrism results.xlsx", sheet = 4)
colnames(cgga_b_immune)[1] <- "patient_id"

setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Sample GSEA")
cgga1_gsea <- read_xlsx("cgga_batch2_GSEA.xlsx", sheet = 1)

cgga_b_merged <- cgga_patient_b %>%
  #merge(tcga_subtype, by = "sample_id", all.x = TRUE) %>%
  merge(cgga_b_genes, by = "patient_id", all.x = TRUE) %>%
  merge(cgga1_gsea, by = "patient_id", all.x = TRUE) %>%
  merge(cgga_b_immune, by = "patient_id", all.x = TRUE) %>%
  filter(!is.na(Oligodendrocyte))

cgga_b_merged <- cgga_b_merged %>%
  dplyr::rename(STING1_tpm = TMEM173, ARID1A_tpm = ARID1A)

#G-SAM cleaning####
setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Clinical Data")
gsam_patient <- read_excel("GSAM Clinical.xlsx") #obtained from Table S1 from orignal paper
gsam_patient <- gsam_patient %>%
  dplyr::mutate(sid = str_replace(sid, ".n", "-new")) %>%
  dplyr::select(sid, pid, resection, svvl.stat, svvl.r1.days, 
                svvl.r2.days, pfs.days, `ssGSEA subtype`, MGMT, TMZ, RT) %>%
  dplyr::rename(sample_id = sid, patient_id = pid, MGMT_status = MGMT, chemo = TMZ, radiation = RT,
                OS_event = svvl.stat, OS_years = svvl.r1.days,
                recurrence = resection, PFS_years = pfs.days, subtype = `ssGSEA subtype`) %>%
  mutate(OS_years = as.numeric(OS_years) / 365.25) %>%
  mutate(PFS_years = as.numeric(PFS_years) / 365.25) %>%
  mutate(race = NA, age = NA, sex = NA, grade = NA, DFS_event = NA, DFS_years = NA, PFS_event = "recurrence",
         location = NA, STING_expr = NA, ARID1A_expr = NA, IDH_status = "WT", 
         KarnPerfScore = NA, resection = NA, study = "GSAM", platform = "rnaseq") %>%
  dplyr::select(patient_id, sample_id, age, sex, race, grade,  
                OS_event, OS_years, DFS_event, DFS_years, PFS_event, PFS_years,
                recurrence, location, resection, chemo, radiation, KarnPerfScore, 
                MGMT_status, IDH_status, study, platform)

setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Instaprism Results")
gsam_immune <- read_xlsx("GSAM RNAseq InstaPrism results.xlsx", sheet = 4)
colnames(gsam_immune)[1] <- "sample_id"

setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Sample GSEA")
gsam_gsea <- read_xlsx("gsam_GSEA.xlsx", sheet = 1)

setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Expression Data")
gsam_tpm <- fread("gsam_TPM.csv")
gsam_tpm <- t(gsam_tpm)
colnames(gsam_tpm) <- gsam_tpm[1, ]
gsam_tpm <- gsam_tpm[-1, ]

gene_names_upper <- toupper(colnames(gsam_tpm))
sting_match <- gene_names_upper %in% toupper(sting_aliases)
arid1a_match <- gene_names_upper %in% toupper(arid1a_aliases)
gsam_genes <- gsam_tpm[, sting_match | arid1a_match, drop = FALSE]
gsam_genes <- as.data.frame(gsam_genes)
gsam_genes$sample_id <- as.character(rownames(gsam_tpm))

gsam_merged <- gsam_patient %>%
  #merge(tcga_subtype, by = "sample_id", all.x = TRUE) %>%
  merge(gsam_genes, by = "sample_id", all.x = TRUE) %>%
  merge(gsam_gsea, by = "sample_id", all.x = TRUE) %>%
  merge(gsam_immune, by = "sample_id", all.x = TRUE) %>%
  dplyr::filter(!is.na(Oligodendrocyte))
gsam_merged <- gsam_merged[, c("patient_id", "sample_id", setdiff(names(gsam_merged), c("patient_id", "sample_id")))]

gsam_merged <- gsam_merged %>%
  dplyr::rename(STING1_tpm = STING1, ARID1A_tpm = ARID1A)

#IVYGAP cleaning####
#input anatomic tumor locations
setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/IVYGAP")
ivygap_matching <- read_excel("columns-samples.xlsx")
ivygap_matching <- ivygap_matching %>%
  dplyr::rename("sample_id" = "rna_well_id")%>%
  mutate(IDH_status = if_else(tumor_name %in% c("W10-1-1", "W4-1-1", "W31-1-1"), "MUT", "WT"))

setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Clinical Data")
ivygap_patient <- fread("IVYGAP Clinical.csv") #obtained from rnaseq sample data from website
ivygap_patient <- ivygap_patient %>%
  merge(ivygap_matching, by = "sample_id", all.x = TRUE) %>%
  dplyr::select(donor_id, sample_id, anatomic_location, survival_days, surgery,
                molecular_subtype, survival_days, age_in_years, anatomic_location, 
                IDH_status, mgmt_methylation, initial_kps, extent_of_resection) %>%
  dplyr::rename(patient_id = donor_id, subtype = molecular_subtype, age = age_in_years,
                recurrence = surgery, location = anatomic_location, OS_years = survival_days,
                MGMT_status = mgmt_methylation, KarnPerfScore = initial_kps, resection = extent_of_resection) %>%
  dplyr::mutate(OS_years = as.numeric(OS_years) / 365.25) %>%
  mutate(race = NA, sex = NA, grade = NA, DFS_event = NA, DFS_years = NA, OS_event = "died", PFS_event = NA,
         PFS_years = NA, STING_expr = NA, chemo = NA, 
         radiation = NA, ARID1A_expr = NA, study = "IVYGAP", platform = "rnaseq") %>%
  dplyr::select(patient_id, sample_id, age, sex, race, grade,  
                OS_event, OS_years, DFS_event, DFS_years, PFS_event, PFS_years,
                recurrence, location, resection, chemo, radiation, KarnPerfScore, 
                MGMT_status, IDH_status, study, platform) %>%
  mutate(OS_event = if_else(is.na(OS_years), NA_character_, OS_event))

setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Instaprism Results")
ivygap_immune <- read_xlsx("IVYGAP RNAseq InstaPrism results.xlsx", sheet = 4)
colnames(ivygap_immune)[1] <- "sample_id"

setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Sample GSEA")
ivygap_gsea <- read_xlsx("IVYGAP_GSEA.xlsx", sheet = 1)

ivygap_patient$sample_id <- as.character(ivygap_patient$sample_id)
ivygap_immune$sample_id <- as.character(ivygap_immune$sample_id)
ivygap_immune$sample_id <- gsub("^X", "", ivygap_immune$sample_id)
ivygap_gsea$sample_id <- gsub("^X", "", ivygap_gsea$sample_id)

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
ivygap_tpm <- as.data.frame(t(tpm_df))

gene_names_upper <- toupper(colnames(ivygap_tpm))
sting_match <- gene_names_upper %in% toupper(sting_aliases)
arid1a_match <- gene_names_upper %in% toupper(arid1a_aliases)
ivygap_genes <- ivygap_tpm[, sting_match | arid1a_match, drop = FALSE]
ivygap_genes$sample_id <- as.character(rownames(ivygap_tpm))
ivygap_genes$sample_id <- gsub("^X", "", ivygap_genes$sample_id)

ivygap_merged <- ivygap_patient %>%
  #merge(tcga_subtype, by = "sample_id", all.x = TRUE) %>%
  merge(ivygap_genes, by = "sample_id", all.x = TRUE) %>%
  merge(ivygap_gsea, by = "sample_id", all.x = TRUE) %>%
  merge(ivygap_immune, by = "sample_id", all.x = TRUE) %>%
  dplyr::filter(!is.na(Oligodendrocyte))
ivygap_merged <- as.data.frame(ivygap_merged)
ivygap_merged <- ivygap_merged[, c("patient_id", "sample_id", setdiff(names(ivygap_merged), c("patient_id", "sample_id")))]

ivygap_merged <- ivygap_merged %>%
  dplyr::rename(STING1_tpm = TMEM173, ARID1A_tpm = ARID1A) 

#GLASS cleaning####
#input anatomic tumor locations
setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Clinical Data")
glass_surgeries <- fread("GLASS Surgeries.csv")
glass_cases <- fread("GLASS Cases.csv")
glass_patient <- merge(glass_surgeries, glass_cases, by = "case_barcode", all.x = TRUE)
glass_patient <- glass_patient %>%
  filter(histology == "Glioblastoma") %>%
  filter(!grepl("tcga", sample_barcode, ignore.case = TRUE) &
        !grepl("tcga", case_barcode, ignore.case = TRUE)) %>%
  dplyr::select(case_barcode, sample_barcode, surgery_number,
                grade, idh_status, case_sex, case_age_diagnosis_years,
                case_vital_status, case_overall_survival_mo, mgmt_methylation, surgery_extent_of_resection,
                treatment_tmz, treatment_radiotherapy) %>%
  dplyr::rename(patient_id = case_barcode, sample_id = sample_barcode,
                IDH_status = idh_status, sex = case_sex, age = case_age_diagnosis_years,
                OS_event = case_vital_status, OS_years = case_overall_survival_mo, 
                recurrence = surgery_number, chemo = treatment_tmz, radiation = treatment_radiotherapy,
                MGMT_status = mgmt_methylation, resection = surgery_extent_of_resection) %>%
  dplyr::mutate(OS_years = as.numeric(OS_years) / 12) %>%
  mutate(race = NA, DFS_event = NA, DFS_years = NA, PFS_event = NA,
         PFS_years = NA, STING_expr = NA, ARID1A_expr = NA, KarnPerfScore = NA, 
         location = NA, study = "GLASS", platform = "rnaseq") %>%
  dplyr::select(patient_id, sample_id, age, sex, race, grade,  
                OS_event, OS_years, DFS_event, DFS_years, PFS_event, PFS_years,
                recurrence, location, resection, chemo, radiation, KarnPerfScore, 
                MGMT_status, IDH_status, study, platform)

  setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Instaprism Results")
glass_immune <- read_xlsx("GLASS RNAseq InstaPrism results.xlsx", sheet = 4)
colnames(glass_immune)[1] <- "sample_id"
glass_immune <- glass_immune %>% dplyr::mutate(sample_id = str_replace_all(sample_id, "\\.", "-"))


setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Sample GSEA")
glass_gsea <- read_xlsx("glass_GSEA.xlsx", sheet = 1)
glass_gsea <- glass_gsea %>% dplyr::mutate(sample_id = str_replace_all(sample_id, "\\.", "-"))

setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Clinical Data")
glass_subtype <- read_xlsx("GLASS Subtype.xlsx")
glass_subtype <- glass_subtype %>%
  filter(p_value < 0.05) %>%
  group_by(aliquot_barcode) %>%
  summarise(signatures = list(signature_name)) %>%
  dplyr::rename(subtype = signatures, sample_id = aliquot_barcode) %>%
  dplyr::mutate(across(where(is.list), ~ sapply(., function(x) paste(x, collapse = "; "))))
glass_subtype <- as.data.frame(glass_subtype)
glass_test <- merge(glass_subtype, glass_gsea, by = "sample_id", all.x = TRUE)
glass_test <- merge(glass_test, glass_immune, by = "sample_id", all.x = TRUE)

glass_test <- glass_test %>%
  dplyr::mutate(sample_id = str_remove(sample_id, "-RNA-.*")) %>%
  dplyr::filter(!grepl("tcga", sample_id, ignore.case = TRUE)) %>%
  filter(!str_detect(sample_id, "-\\d{2}R$") | str_detect(sample_id, "-01R$")) %>%  # Keep only -01R samples
  mutate(sample_id = str_sub(sample_id, 1, -5))

setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Expression Data")
df <- fread("GLASS Expression.tsv")
df <- df %>%
  filter(!is.na(Gene_symbol)) %>%
  distinct(Gene_symbol, .keep_all = TRUE)
df <- column_to_rownames(df, var = "Gene_symbol")
glass_tpm <- as.data.frame(t(df))

gene_names_upper <- toupper(colnames(glass_tpm))
sting_match <- gene_names_upper %in% toupper(sting_aliases)
arid1a_match <- gene_names_upper %in% toupper(arid1a_aliases)
glass_genes <- glass_tpm[, sting_match | arid1a_match, drop = FALSE]
glass_genes$sample_id <- as.character(rownames(glass_tpm))
glass_genes <- glass_genes %>% dplyr::mutate(sample_id = str_replace_all(sample_id, "\\.", "-"))
glass_genes <- glass_genes %>%
  dplyr::mutate(sample_id = str_remove(sample_id, "-RNA-.*")) %>%
  dplyr::filter(!grepl("tcga", sample_id, ignore.case = TRUE)) %>%
  filter(!str_detect(sample_id, "-\\d{2}R$") | str_detect(sample_id, "-01R$")) %>%  # Keep only -01R samples
  mutate(sample_id = str_sub(sample_id, 1, -5))

glass_merged <- glass_patient %>%
  #merge(tcga_subtype, by = "sample_id", all.x = TRUE) %>%
  merge(glass_genes, by = "sample_id", all.x = TRUE) %>%
  merge(glass_test, by = "sample_id", all.x = TRUE) %>%
  dplyr::filter(!is.na(Oligodendrocyte))
glass_merged <- as.data.frame(glass_merged)
glass_merged <- glass_merged[, c("patient_id", "sample_id", setdiff(names(glass_merged), c("patient_id", "sample_id")))]

glass_merged <- glass_merged %>%
  dplyr::rename(STING1_tpm = TMEM173, ARID1A_tpm = ARID1A) %>%
  dplyr::select(-subtype)

#TCGA UG133A####
setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/GlioVis Datasets")
TCGA_GBM <- readRDS("TCGA_GBM.Rds")
tcga_karn <- fread("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Clinical Data/tcga_Karnofsky_Performance_Score.txt")
tcga_microarray_clinical <- TCGA_GBM[["pData"]] %>%
  mutate(Sample = gsub("\\.", "-", Sample)) %>%
  merge(tcga_karn, by.x = "Sample", by.y = "Patient ID", all.x = TRUE) %>%
  dplyr::filter(Histology == "GBM") %>%
  dplyr::rename(patient_id = Sample, grade = Grade, age = Age, sex = Gender, 
                recurrence = Recurrence, OS_years = survival, OS_event = status,
                IDH_status = IDH1_status, chemo = Therapy_Class, 
                KarnPerfScore = "Karnofsky Performance Score") %>%
  mutate(OS_years = as.numeric(OS_years) / 12) %>%
  mutate(sample_id = NA, DFS_event = NA, DFS_years = NA, PFS_event = NA, 
         PFS_years = NA, resection = NA, radiation = NA, location = NA,
         chemo = NA, race = NA, ARID1A_tpm = NA, 
         STING1_tpm = NA, study = "TCGA_microarray", platform = "microarray") %>%
  dplyr::select(patient_id, sample_id, age, sex, race, grade,  
                  OS_event, OS_years, DFS_event, DFS_years, PFS_event, PFS_years,
                  recurrence, location, resection, chemo, radiation, KarnPerfScore, 
                  MGMT_status, IDH_status, study, platform, ARID1A_tpm, STING1_tpm)
      
setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Instaprism Results")
microarray_immune <- read_xlsx("TCGA Microarray InstaPrism results.xlsx", sheet = 4)
colnames(microarray_immune)[1] <- "patient_id"
microarray_immune <- mutate(microarray_immune, patient_id = gsub("\\.", "-", patient_id))

setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Sample GSEA")
microarray_gsea <- read_xlsx("tcga_microarray_GSEA.xlsx", sheet = 1) %>%
  dplyr::rename(patient_id = sample_id) %>%
  mutate(patient_id = gsub("\\.", "-", patient_id))

tcga_microarray_merged <- tcga_microarray_clinical %>%
  merge(microarray_gsea, by = "patient_id", all.x = TRUE) %>%
  merge(microarray_immune, by = "patient_id", all.x = TRUE) %>%
  dplyr::filter(!is.na(Oligodendrocyte))


#combine datasets####
#cleaning data types
clean_os_event <- function(x) {
  x_chr <- tolower(as.character(x))
  case_when(
    str_detect(x_chr, "deceased|dead|died")  ~ 1L,
    str_detect(x_chr, "alive|living|censored") ~ 0L,
    x_chr == "1" ~ 1L,
    x_chr == "0" ~ 0L,
    TRUE       ~ NA_integer_
  )
}
clean_dfs_event <- function(x) {
  x_chr <- tolower(as.character(x))
  case_when(
    str_detect(x_chr, "recurred|relapse|progression|event") ~ 1L,
    str_detect(x_chr, "no event|disease-free|censored")    ~ 0L,
    x_chr == "1" ~ 1L,
    x_chr == "0" ~ 0L,
    TRUE         ~ NA_integer_
  )
}
clean_pfs_event <- function(x) {
  x_chr <- tolower(as.character(x))
  case_when(
    str_detect(x_chr, "progressed|recurrence|event")  ~ 1L,
    str_detect(x_chr, "no event|stable|censored")     ~ 0L,
    x_chr == "1" ~ 1L,
    x_chr == "0" ~ 0L,
    TRUE         ~ NA_integer_
  )
}
clean_recurrence <- function(df, cohort_name) {
  df %>%
    # coerce to plain lower‐case text
    mutate(
      .rec = tolower(as.character(recurrence))
    ) %>%
    mutate(
      recurrence = case_when(
        ## TCGA: anything with "primary" → "primary", else NA
        cohort_name == "TCGA" & str_detect(.rec, "primary") ~ "primary",
        cohort_name == "TCGA" & str_detect(.rec, "recurrence") ~ "recurrent",
        
        cohort_name == "TCGA_microarray" & str_detect(.rec, "primary") ~ "primary",
        cohort_name == "TCGA_microarray" & str_detect(.rec, "recurrent") ~ "recurrent",
        cohort_name == "TCGA_microarray" & str_detect(.rec, "secondary") ~ "secondary",
        
        ## Tempus & WU (and any other single‐label "primary" cohorts)
        cohort_name %in% c("TEMPUS", "WU") & .rec == "primary" ~ "primary",
        
        ## IVYGAP: "primary"/"recurrent"
        cohort_name == "IVYGAP" & .rec %in% c("primary", "recurrent") ~ .rec,
        
        ## CGGA (no secondary)
        cohort_name == "CGGA" & .rec %in% c("primary", "recurrent") ~ .rec,
        
        ## CGGA_B (has secondary)
        cohort_name == "CGGA_B" & .rec %in% c("primary", "recurrent", "secondary") ~ .rec,
        
        ## GSAM: r1 → primary, r2 → recurrent
        cohort_name == "GSAM" & .rec == "r1" ~ "primary",
        cohort_name == "GSAM" & .rec == "r2" ~ "recurrent",
        
        ## GLASS: numeric 1 → primary; 2/3/4 → recurrent
        cohort_name == "GLASS" & .rec == "1"                      ~ "primary",
        cohort_name == "GLASS" & .rec %in% c("2", "3", "4")       ~ "recurrent",
        
        ## all else → NA
        TRUE ~ NA_character_
      )
    ) %>%
    dplyr::select(-.rec)
}
clean_radiation <- function(x) {
  x_chr <- tolower(as.character(x))
  dplyr::case_when(
    x_chr %in% c("yes", "y", "1")          ~ 1L,
    x_chr %in% c("no",  "n", "0")          ~ 0L,
    x_chr %in% c("[not available]", "", NA_character_) ~ NA_integer_,
    TRUE                                   ~ NA_integer_
  )
}
clean_chemo <- function(x) {
  x_chr <- tolower(as.character(x))
  dplyr::case_when(
    x_chr %in% c("tmz/radiotherapy")      ~ "comb",  # specific combo first
    x_chr %in% c("rt", "radiotherapy")    ~ "rad",
    x_chr %in% c("yes", "y", "1", "ct")   ~ "tmz",
    x_chr %in% c("no", "n", "0")          ~ "none",
    x_chr %in% c("[not available]", "", NA) ~ NA_character_,
    TRUE                                   ~ NA_character_
  )
}
clean_karnperfscore <- function(x) {
  x_chr  <- as.character(x)
  digits <- gsub("[^0-9]", "", x_chr)
  as.integer(ifelse(digits == "", NA, digits))
}
clean_resection <- function(x) {
  x_chr <- tolower(as.character(x))
  dplyr::case_when(
    x_chr %in% c("biopsy")                   ~ "Biopsy",
    x_chr %in% c("sub-total", "subtotal")    ~ "Subtotal",
    x_chr %in% c("complete", "total")        ~ "Total",
    TRUE                                     ~ NA_character_
  )
}
clean_mgmt_status <- function(x) {
  x_chr <- tolower(as.character(x))
  dplyr::case_when(
    x_chr %in% c(
      "hypermethylated", "methylated", "Meth", "M", "Yes", "Methylated"
    ) ~ "hypermethylated",
    x_chr %in% c(
      "hypomethylated", "No", "U", "un-methylated", "Unmeth", "unmethylated", "Unmethylated"
    ) ~ "hypomethylated",
    TRUE ~ NA_character_  
  )
}
clean_all <- function(df, cohort_name) {
  df %>%
    mutate(
      patient_id   = as.character(patient_id),
      # strip non-numeric text (units, commas, etc.) then parse
      age          = as.double(gsub("[^0-9.]", "", age)),
      ARID1A_tpm   = as.double(gsub("[^0-9.]", "", ARID1A_tpm)),
      STING1_tpm   = as.double(gsub("[^0-9.]", "", STING1_tpm)),
      OS_event     = clean_os_event(OS_event),
      DFS_event    = clean_dfs_event(DFS_event),
      PFS_event    = clean_pfs_event(PFS_event),
      radiation    = clean_radiation(radiation),
      chemo        = clean_chemo(chemo),
      MGMT_status  = clean_mgmt_status(MGMT_status),
      KarnPerfScore  = clean_karnperfscore(KarnPerfScore),
      resection      = clean_resection(resection)
    ) %>%
    clean_recurrence(cohort_name)
}

tcga_merged_clean   <- clean_all(tcga_merged,   "TCGA")
tempus_merged_clean <- clean_all(tempus_merged, "TEMPUS")
wu_merged_clean     <- clean_all(wu_merged,     "WU")
cgga_merged_clean   <- clean_all(cgga_merged,   "CGGA")
cgga_b_merged_clean <- clean_all(cgga_b_merged, "CGGA_B")
gsam_merged_clean   <- clean_all(gsam_merged,   "GSAM")
ivygap_merged_clean <- clean_all(ivygap_merged, "IVYGAP")
glass_merged_clean  <- clean_all(glass_merged,  "GLASS")
tcga_microarray_merged_clean  <- clean_all(tcga_microarray_merged,  "TCGA_microarray")

combined_dataset <- bind_rows(
  tcga_merged_clean,
  tempus_merged_clean,
  wu_merged_clean,
  cgga_merged_clean,
  cgga_b_merged_clean,
  gsam_merged_clean,
  ivygap_merged_clean,
  glass_merged_clean,
  tcga_microarray_merged_clean
)

combined_dataset <- combined_dataset %>%
  dplyr::mutate(
    treatment = dplyr::case_when(
      chemo == "none" & radiation == 1 ~ "rad",
      chemo == "none" & radiation == 0 ~ "none",
      chemo == "rad" ~ "rad",
      chemo == "comb" ~ "tmz/rad",
      chemo == "tmz" & radiation == 0 ~ "tmz",
      chemo == "tmz" & radiation == 1 ~ "tmz/rad",
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::select(-c(chemo, radiation))

combined_dataset <- combined_dataset %>%
  dplyr::select(
    patient_id, sample_id, age, sex, race, grade, OS_event, OS_years, DFS_event,
    DFS_years, PFS_event, PFS_years, recurrence, location, resection, KarnPerfScore,
    MGMT_status, IDH_status, treatment, everything()
  )

table(tcga_merged$treatment)
table(tempus_merged$treatment)
table(wu_merged$treatment)
table(cgga_merged$treatment)
table(cgga_b_merged$treatment)
table(gsam_merged$treatment)
table(ivygap_merged$treatment)
table(glass_merged$treatment)
table(tcga_microarray_merged$treatment)

#sample_id
combined_dataset <- combined_dataset %>%
  mutate(
    sample_id = if_else(
      is.na(sample_id) | sample_id == "",
      patient_id,
      sample_id))

#Sex
combined_dataset <- combined_dataset %>%
  mutate(
    sex = case_when(
      tolower(sex) == "male"   ~ "Male",
      tolower(sex) == "female" ~ "Female",
      tolower(sex) %in% c("unknown", "[not available]", "", NA_character_) ~ NA_character_,
      TRUE                      ~ NA_character_
    )
  )

#race
combined_dataset <- combined_dataset %>%
  mutate(
    race = case_when(
      tolower(race) == "white"                     ~ "White",
      tolower(race) == "black or african american" ~ "Black or African American",
      tolower(race) == "asian"                     ~ "Asian",
      tolower(race) %in% c("unknown", "[not available]", "other race", "", "other") ~ NA_character_,
      TRUE                                         ~ NA_character_
    )
  )

#grade
combined_dataset <- combined_dataset %>%
  mutate(
    grade = case_when(
      tolower(grade) %in% c("iv", "who iv") ~ "IV",
      TRUE ~ NA_character_
    )
  )

#idh mutatiion
combined_dataset <- combined_dataset %>%
  mutate(
    IDH_status = tolower(IDH_status),
    IDH_status = case_when(
      IDH_status %in% c("wt", "wildtype", "idhwt", "wild-type") ~ "wt",
      IDH_status %in% c("r132h", "mutant", "mut", "idhmut", "Mutant") ~ "mut",
      TRUE ~ NA_character_
    )
  )

#cell types
combined_dataset <- as.data.frame(combined_dataset)

new_names <- names(combined_dataset)
new_names <- str_replace_all(new_names, "\\-", "")       # Remove hyphens
new_names <- str_replace_all(new_names, "\\/", "_")       # Remove slashes
new_names <- str_replace_all(new_names, " ", "_")         # Replace spaces with underscores
new_names <- str_replace_all(new_names, "\\+", "")         # Remove plus signs
combined_dataset <- setNames(combined_dataset, new_names)
print(names(combined_dataset))

opc_like_columns <- grep("^OPClike", names(combined_dataset), value = TRUE)
ac_like_columns  <- grep("^AClike", names(combined_dataset), value = TRUE)
npc_like_columns <- grep("^NPClike", names(combined_dataset), value = TRUE)
mes_like_columns <- grep("^MESlike", names(combined_dataset), value = TRUE)

combined_dataset$OPC_like <- rowSums(combined_dataset[opc_like_columns], na.rm = TRUE)
combined_dataset$AC_like <- rowSums(combined_dataset[ac_like_columns], na.rm = TRUE)
combined_dataset$NPC_like <- rowSums(combined_dataset[npc_like_columns], na.rm = TRUE)
combined_dataset$MES_like <- rowSums(combined_dataset[mes_like_columns], na.rm = TRUE)
combined_dataset <- combined_dataset[, !(names(combined_dataset) %in% c(opc_like_columns, ac_like_columns, npc_like_columns, mes_like_columns))]

#rename cel types
# Rename and combine columns in the dataset
# Define mapping: original names to new names
rename_map <- c(
  "Oligodendrocyte" = "Oligodendrocyte",
  "Pericyte" = "Pericyte",
  "Neuron" = "Neuron",
  "Astrocyte" = "Astrocyte",
  "Endo_arterial" = "Arterial_Vessel",
  "Tiplike" = "Tiplike_Vessel",
  "Endo_capilar" = "Capilary_Vessel",
  "TAMBDM_hypoxia_MES" = "BDM_Hypoxia_and_MES_like",
  "cDC2" = "Conventional_DC_2",
  "TAMMG_aging_sig" = "Microglia_Aging_Signature",
  "TAMBDM_antiinfl" = "Anti_inflammatory_BDM",
  "Mono_antiinfl" = "Anti_inflammatory_Monocyte",
  "Mono_hypoxia" = "Hypoxia_associated_Monocytes",
  "Mono_naive" = "Naive_Monocyte",
  "TAMMG_proinfl_I" = "Pro_inflammatory_Microglia",
  "TAMBDM_MHC" = "BDM_MHC_Expression",
  "TAMMG_proinfl_II" = "Pro_inflammatory_Microglia",
  "TAMMG_prolif" = "Proliferative_Microglia",
  "TAMBDM_INF" = "BDM_IFN_Signature",
  "Mast" = "Mast_Cell",
  "DC3" = "Dendritic_like_Cell",
  "DC1" = "Dendritic_like_Cell",
  "DC2" = "Dendritic_like_Cell",
  "cDC1" = "Conventional_DC_1",
  "pDC" = "Plasmacytoid_DC",
  "Stress_sig" = "T_Cells_Stress_Signature",
  "Prolif_T" = "Proliferating_T_Cells",
  "CD8_EM" = "CD8_Effector_Memory_T_Cells",
  "CD8_NK_sig" = "NK_like_CD8_T_Cells",
  "Reg_T" = "Regulatory_T_Cells",
  "CD8_cytotoxic" = "Cytotoxic_CD8_T_Cell",
  "B_cell" = "B_Cell",
  "CD4_rest" = "Resting_CD4_T_Cells",
  "CD4_INF" = "CD4_T_Cell_IFN_Signature",
  "Plasma_B" = "Plasma_Cell",
  "NK" = "Natural_Killer_Cells",
  "OPC" = "OPC",
  "RG" = "Radial_Glia",
  "OPC_like" = "OPC_like_Cancer",
  "AC_like" = "AC_like_Cancer",
  "NPC_like" = "NPC_like_Cancer",
  "MES_like" = "MES_like_Cancer"
)

metadata_cols      <- combined_dataset[, 1:28]
celltype_matrix    <- combined_dataset[, names(rename_map)]
colnames(celltype_matrix) <- rename_map[colnames(celltype_matrix)]

new_names <- colnames(celltype_matrix)
unique_names <- unique(new_names)

celltype_summed <- as.data.frame(
  sapply(unique_names, function(nm) {
    # select the columns corresponding to this new name
    cols <- celltype_matrix[, new_names == nm, drop = FALSE]
    # rowSums handles single‐column and multi‐column cases
    rowSums(cols)
  }),
  row.names = NULL
)

colnames(celltype_summed) <- unique_names
combined_summed <- cbind(metadata_cols, celltype_summed)

combined_summed <- combined_summed %>% 
  mutate(MGMT_status = clean_mgmt_status(MGMT_status))

#remove duplicates ids between tcgas
# First get TCGA sample IDs
tcga_samples <- combined_summed %>%
  filter(study == "TCGA") %>%
  pull(patient_id)

# Now filter out duplicate IDs from TCGA_microarray
combined_summed <- combined_summed %>%
  filter(!(study == "TCGA_microarray" & sample_id %in% tcga_samples))

# Save the combined dataset
setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME")
write.csv(combined_summed, "GBM Clinical-TME Cell States.csv", row.names = FALSE)

#Huang2021 cleaning####
#Gill2014 cleaning####
#Claughtery cleaning####
#input anatomic tumor locations
setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Clinical Data")
cloughesy_patient <- read_excel("Cloughesy Clinical.xlsx")
cloughesy_patient <- glass_patient %>%
  dplyr::select("Subject ID", "Treatment Group", "Progression Free Survival (Days)",
                "PFS Censored Status (0 = PF)", "Overall Survival (Days)", 
                "OS Censored Status (0 = Alive, 1 = Dead)", Sex, Race,
                Age, MGMT, IDH, "KPS at registration")%>%
  dplyr::rename(patient_id = "Subject ID",
                IDH_status = IDH, sex = Sex, age = Age,
                OS_event = "OS Censored Status (0 = Alive, 1 = Dead)", 
                OS_years = "Overall Survival (Days)", 
                chemo = "Treatment Group", 
                PFS_event = "PFS Censored Status (0 = PF)",
                PFS_years = "Progression Free Survival (Days)",
                KarnPerfScore = "KPS at registration",
                MGMT_status = MGMT, 
                IDH_status = IDH, 
                race = Race) %>%
  dplyr::mutate(OS_years = as.numeric(OS_years) / 365.25) %>%
  dplyr::mutate(PFS_years = as.numeric(OS_years) / 365.25) %>%
  mutate(DFS_event = NA, DFS_years = NA, STING_expr = NA, ARID1A_expr = NA, 
         location = NA, radiation = NA, study = "Cloughesy", subtype = NA) %>%
  dplyr::select(patient_id, sample_id, age, sex, race, grade,  
                OS_event, OS_years, DFS_event, DFS_years, PFS_event, PFS_years,
                recurrence, location, resection, chemo, radiation, KarnPerfScore, 
                MGMT_status, IDH_status, study, subtype)

setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Instaprism Results")
cloughesy_immune <- read_xlsx("Cloughesy RNAseq InstaPrism results.xlsx", sheet = 4)
colnames(cloughesy_immune)[1] <- "patient_id"

setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Expression Data")
df <- fread("Cloughesy TPM Expression.tsv")
df <- df %>%
  filter(!is.na(Gene_symbol)) %>%
  distinct(Gene_symbol, .keep_all = TRUE)
df <- column_to_rownames(df, var = "Gene_symbol")
cloughesy_tpm <- as.data.frame(t(df))

gene_names_upper <- toupper(colnames(glass_tpm))
sting_match <- gene_names_upper %in% toupper(sting_aliases)
arid1a_match <- gene_names_upper %in% toupper(arid1a_aliases)
glass_genes <- glass_tpm[, sting_match | arid1a_match, drop = FALSE]
glass_genes$sample_id <- as.character(rownames(glass_tpm))
glass_genes <- glass_genes %>% dplyr::mutate(sample_id = str_replace_all(sample_id, "\\.", "-"))
glass_genes <- glass_genes %>%
  dplyr::mutate(sample_id = str_remove(sample_id, "-RNA-.*")) %>%
  dplyr::filter(!grepl("tcga", sample_id, ignore.case = TRUE)) %>%
  filter(!str_detect(sample_id, "-\\d{2}R$") | str_detect(sample_id, "-01R$")) %>%  # Keep only -01R samples
  mutate(sample_id = str_sub(sample_id, 1, -5))

glass_merged <- glass_patient %>%
  #merge(tcga_subtype, by = "sample_id", all.x = TRUE) %>%
  merge(glass_genes, by = "sample_id", all.x = TRUE) %>%
  merge(glass_test, by = "sample_id", all.x = TRUE) %>%
  dplyr::filter(!is.na(Oligodendrocyte))
glass_merged <- as.data.frame(glass_merged)
glass_merged <- glass_merged[, c("patient_id", "sample_id", setdiff(names(glass_merged), c("patient_id", "sample_id")))]

glass_merged <- glass_merged %>%
  dplyr::rename(STING1_tpm = TMEM173, ARID1A_tpm = ARID1A) %>%
  dplyr::select(-subtype)


#Rembrandt####
setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/GlioVis Datasets")
GBM <- readRDS("Rembrandt.Rds")
microarray_clinical <- GBM[["pData"]] %>%
  dplyr::filter(Histology == "GBM") %>%
  dplyr::rename(patient_id = Sample, grade = Grade, 
                recurrence = Recurrence, OS_years = survival, OS_event = status) %>%
  mutate(OS_years = as.numeric(OS_years) / 12) %>%
  mutate(sample_id = NA, DFS_event = NA, DFS_years = NA, PFS_event = NA, 
         PFS_years = NA, resection = NA, radiation = NA, location = NA,
         chemo = NA, race = NA, ARID1A_tpm = NA, STING1_tpm = NA, KarnPerfScore = NA,
         study = "Rembrandt", IDH_status = NA, MGMT_status = NA, age = NA,
         race = NA, sex = NA, platform = "microarray") %>%
  dplyr::select(patient_id, sample_id, age, sex, race, grade,  
                OS_event, OS_years, DFS_event, DFS_years, PFS_event, PFS_years,
                recurrence, location, resection, chemo, radiation, KarnPerfScore, 
                MGMT_status, IDH_status, study, platform, ARID1A_tpm, STING1_tpm)

setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Instaprism Results")
microarray_immune <- read_xlsx("Rembrandt Microarray InstaPrism results.xlsx", sheet = 4)
colnames(microarray_immune)[1] <- "patient_id"

setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Sample GSEA")
microarray_gsea <- read_xlsx("Rembrandt_microarray_GSEA.xlsx", sheet = 1) %>%
  rename(patient_id = sample_id)

Rembrandt_merged <- microarray_clinical %>%
  merge(microarray_gsea, by = "patient_id", all.x = TRUE) %>%
  merge(microarray_immune, by = "patient_id", all.x = TRUE) %>%
  dplyr::filter(!is.na(Oligodendrocyte))


#Gravendeel ####
setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/GlioVis Datasets")
GBM <- readRDS("Gravendeel.Rds")
microarray_clinical <- GBM[["pData"]] %>%
  dplyr::filter(Histology == "GBM") %>%
  dplyr::rename(patient_id = Sample, grade = Grade, 
                recurrence = Recurrence, OS_years = survival, OS_event = status) %>%
  mutate(OS_years = as.numeric(OS_years) / 12) %>%
  mutate(sample_id = NA, DFS_event = NA, DFS_years = NA, PFS_event = NA, 
         PFS_years = NA, resection = NA, radiation = NA, location = NA,
         chemo = NA, race = NA, ARID1A_tpm = NA, STING1_tpm = NA, KarnPerfScore = NA,
         study = "Gravendeel", IDH_status = NA, MGMT_status = NA, age = NA,
         race = NA, sex = NA, platform = "microarray") %>%
  dplyr::select(patient_id, sample_id, age, sex, race, grade,  
                OS_event, OS_years, DFS_event, DFS_years, PFS_event, PFS_years,
                recurrence, location, resection, chemo, radiation, KarnPerfScore, 
                MGMT_status, IDH_status, study, platform, ARID1A_tpm, STING1_tpm)

setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Instaprism Results")
microarray_immune <- read_xlsx("Gravendeel Microarray InstaPrism results.xlsx", sheet = 4)
colnames(microarray_immune)[1] <- "patient_id"

setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Sample GSEA")
microarray_gsea <- read_xlsx("Gravendeel_microarray_GSEA.xlsx", sheet = 1) %>%
  rename(patient_id = sample_id)

Gravendeel_merged <- microarray_clinical %>%
  merge(microarray_gsea, by = "patient_id", all.x = TRUE) %>%
  merge(microarray_immune, by = "patient_id", all.x = TRUE) %>%
  dplyr::filter(!is.na(Oligodendrocyte))

#LeeY ####
setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/GlioVis Datasets")
GBM <- readRDS("LeeY.Rds")
microarray_clinical <- GBM[["pData"]] %>%
  dplyr::filter(Histology == "GBM") %>%
  dplyr::rename(patient_id = Sample, grade = Grade, 
                recurrence = Recurrence, OS_years = survival, OS_event = status,
                study = Dataset_code, sex = Gender, age = Age) %>%
  mutate(OS_years = as.numeric(OS_years) / 12) %>%
  mutate(sample_id = NA, DFS_event = NA, DFS_years = NA, PFS_event = NA, 
         PFS_years = NA, resection = NA, radiation = NA, location = NA,
         chemo = NA, race = NA, ARID1A_tpm = NA, STING1_tpm = NA, KarnPerfScore = NA,
         IDH_status = NA, MGMT_status = NA,race = NA, platform = "microarray") %>%
  dplyr::select(patient_id, sample_id, age, sex, race, grade,  
                OS_event, OS_years, DFS_event, DFS_years, PFS_event, PFS_years,
                recurrence, location, resection, chemo, radiation, KarnPerfScore, 
                MGMT_status, IDH_status, study, platform, ARID1A_tpm, STING1_tpm)
microarray_clinical$study <- paste0("LeeY", microarray_clinical$study)


setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Instaprism Results")
microarray_immune <- read_xlsx("LeeY Microarray InstaPrism results.xlsx", sheet = 4)
colnames(microarray_immune)[1] <- "patient_id"

setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Sample GSEA")
microarray_gsea <- read_xlsx("LeeY_microarray_GSEA.xlsx", sheet = 1) %>%
  rename(patient_id = sample_id)

LeeY_merged <- microarray_clinical %>%
  merge(microarray_gsea, by = "patient_id", all.x = TRUE) %>%
  merge(microarray_immune, by = "patient_id", all.x = TRUE) %>%
  dplyr::filter(!is.na(Oligodendrocyte))


#Phillips ####
setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/GlioVis Datasets")
GBM <- readRDS("Phillips.Rds")
microarray_clinical <- GBM[["pData"]] %>%
  dplyr::select(-Subtype.ssgsea) %>%
  dplyr::filter(Histology == "GBM") %>%
  dplyr::rename(patient_id = Sample, grade = Grade, 
                recurrence = Recurrence, OS_years = survival, OS_event = status,
                sample_id = Sample_ID, age = Age, sex = Gender) %>%
  mutate(OS_years = as.numeric(OS_years) / 12) %>%
  mutate(DFS_event = NA, DFS_years = NA, PFS_event = NA, 
         PFS_years = NA, resection = NA, radiation = NA, location = NA,
         chemo = NA, race = NA, ARID1A_tpm = NA, STING1_tpm = NA, KarnPerfScore = NA,
         study = "Phillips", IDH_status = NA, MGMT_status = NA, race = NA, platform = "microarray") %>%
  dplyr::select(patient_id, sample_id, age, sex, race, grade,  
                OS_event, OS_years, DFS_event, DFS_years, PFS_event, PFS_years,
                recurrence, location, resection, chemo, radiation, KarnPerfScore, 
                MGMT_status, IDH_status, study, platform, ARID1A_tpm, STING1_tpm)

setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Instaprism Results")
microarray_immune <- read_xlsx("Phillips Microarray InstaPrism results.xlsx", sheet = 4)
colnames(microarray_immune)[1] <- "patient_id"

setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Sample GSEA")
microarray_gsea <- read_xlsx("Phillips_microarray_GSEA.xlsx", sheet = 1) %>%
  rename(patient_id = sample_id)

Phillips_merged <- microarray_clinical %>%
  merge(microarray_gsea, by = "patient_id", all.x = TRUE) %>%
  merge(microarray_immune, by = "patient_id", all.x = TRUE) %>%
  dplyr::filter(!is.na(Oligodendrocyte))



#Freije ####
setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/GlioVis Datasets")
GBM <- readRDS("Freije.Rds")
microarray_clinical <- GBM[["pData"]] %>%
  dplyr::filter(Histology == "GBM") %>%
  dplyr::rename(patient_id = Sample, grade = Grade, 
                recurrence = Recurrence, OS_years = survival, OS_event = status,
                age = Age, sex = Gender) %>%
  mutate(OS_years = as.numeric(OS_years) / 12) %>%
  mutate(sample_id = NA, DFS_event = NA, DFS_years = NA, PFS_event = NA, 
         PFS_years = NA, resection = NA, radiation = NA, location = NA,
         chemo = NA, race = NA, ARID1A_tpm = NA, STING1_tpm = NA, KarnPerfScore = NA,
         study = "Freije", IDH_status = NA, MGMT_status = NA, race = NA, platform = "microarray") %>%
  dplyr::select(patient_id, sample_id, age, sex, race, grade,  
                OS_event, OS_years, DFS_event, DFS_years, PFS_event, PFS_years,
                recurrence, location, resection, chemo, radiation, KarnPerfScore, 
                MGMT_status, IDH_status, study, platform, ARID1A_tpm, STING1_tpm)

setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Instaprism Results")
microarray_immune <- read_xlsx("Freije Microarray InstaPrism results.xlsx", sheet = 4)
colnames(microarray_immune)[1] <- "patient_id"

setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Sample GSEA")
microarray_gsea <- read_xlsx("Freije_microarray_GSEA.xlsx", sheet = 1) %>%
  rename(patient_id = sample_id)

Freije_merged <- microarray_clinical %>%
  merge(microarray_gsea, by = "patient_id", all.x = TRUE) %>%
  merge(microarray_immune, by = "patient_id", all.x = TRUE) %>%
  dplyr::filter(!is.na(Oligodendrocyte))



#Murat ####
setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/GlioVis Datasets")
GBM <- readRDS("Murat.Rds")
microarray_clinical <- GBM[["pData"]] %>%
  dplyr::filter(Histology == "GBM") %>%
  dplyr::rename(patient_id = Sample, grade = Grade, 
                recurrence = Recurrence, OS_years = survival, OS_event = status,
                age = Age, sex = Gender, chemo = Therapy) %>%
  mutate(OS_years = as.numeric(OS_years) / 12) %>%
  mutate(sample_id = NA, DFS_event = NA, DFS_years = NA, PFS_event = NA, 
         PFS_years = NA, resection = NA, radiation = NA, location = NA,
         race = NA, ARID1A_tpm = NA, STING1_tpm = NA, KarnPerfScore = NA,
         study = "Murat", IDH_status = NA, platform = "microarray") %>%
  dplyr::select(patient_id, sample_id, age, sex, race, grade,  
                OS_event, OS_years, DFS_event, DFS_years, PFS_event, PFS_years,
                recurrence, location, resection, chemo, radiation, KarnPerfScore, 
                MGMT_status, IDH_status, study, platform, ARID1A_tpm, STING1_tpm)

setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Instaprism Results")
microarray_immune <- read_xlsx("Murat Microarray InstaPrism results.xlsx", sheet = 4)
colnames(microarray_immune)[1] <- "patient_id"

setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Sample GSEA")
microarray_gsea <- read_xlsx("Murat_microarray_GSEA.xlsx", sheet = 1) %>%
  rename(patient_id = sample_id)

Murat_merged <- microarray_clinical %>%
  merge(microarray_gsea, by = "patient_id", all.x = TRUE) %>%
  merge(microarray_immune, by = "patient_id", all.x = TRUE) %>%
  dplyr::filter(!is.na(Oligodendrocyte))



#Joo ####
setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/GlioVis Datasets")
GBM <- readRDS("Joo.Rds")
microarray_clinical <- GBM[["pData"]] %>%
  dplyr::filter(Histology == "GBM") %>%
  dplyr::rename(sample_id = Sample, grade = Grade, 
                recurrence = Recurrence, OS_years = survival, OS_event = status,
                sex = Gender, age = Age, patient_id = Patient_ID) %>%
  mutate(OS_years = as.numeric(OS_years) / 12) %>%
  mutate(DFS_event = NA, DFS_years = NA, PFS_event = NA, 
         PFS_years = NA, resection = NA, radiation = NA, location = NA,
         chemo = NA, race = NA, ARID1A_tpm = NA, STING1_tpm = NA, KarnPerfScore = NA,
         study = "Joo", IDH_status = NA, MGMT_status = NA, platform = "microarray") %>%
  dplyr::select(patient_id, sample_id, age, sex, race, grade,  
                OS_event, OS_years, DFS_event, DFS_years, PFS_event, PFS_years,
                recurrence, location, resection, chemo, radiation, KarnPerfScore, 
                MGMT_status, IDH_status, study, platform, ARID1A_tpm, STING1_tpm)

setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Instaprism Results")
microarray_immune <- read_xlsx("Joo Microarray InstaPrism results.xlsx", sheet = 4)
colnames(microarray_immune)[1] <- "sample_id"

setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Sample GSEA")
microarray_gsea <- read_xlsx("Joo_microarray_GSEA.xlsx", sheet = 1)

Joo_merged <- microarray_clinical %>%
  merge(microarray_gsea, by = "sample_id", all.x = TRUE) %>%
  merge(microarray_immune, by = "sample_id", all.x = TRUE) %>%
  dplyr::filter(!is.na(Oligodendrocyte))



#Ducray ####
setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/GlioVis Datasets")
GBM <- readRDS("Ducray.Rds")
microarray_clinical <- GBM[["pData"]] %>%
  dplyr::filter(Histology == "GBM") %>%
  dplyr::rename(patient_id = Sample, grade = Grade, 
                recurrence = Recurrence, OS_years = survival, OS_event = status,
                age = Age, chemo = Treatment) %>%
  mutate(OS_years = as.numeric(OS_years) / 12) %>%
  mutate(sample_id = NA, DFS_event = NA, DFS_years = NA, PFS_event = NA, 
         PFS_years = NA, resection = NA, radiation = NA, location = NA,
         race = NA, ARID1A_tpm = NA, STING1_tpm = NA, KarnPerfScore = NA,
         study = "Ducray", IDH_status = NA, sex = NA, platform = "microarray") %>%
  dplyr::select(patient_id, sample_id, age, sex, race, grade,  
                OS_event, OS_years, DFS_event, DFS_years, PFS_event, PFS_years,
                recurrence, location, resection, chemo, radiation, KarnPerfScore, 
                MGMT_status, IDH_status, study, platform, ARID1A_tpm, STING1_tpm)

setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Instaprism Results")
microarray_immune <- read_xlsx("Ducray Microarray InstaPrism results.xlsx", sheet = 4)
colnames(microarray_immune)[1] <- "patient_id"

setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Sample GSEA")
microarray_gsea <- read_xlsx("Ducray_microarray_GSEA.xlsx", sheet = 1) %>%
  rename(patient_id = sample_id)

Ducray_merged <- microarray_clinical %>%
  merge(microarray_gsea, by = "patient_id", all.x = TRUE) %>%
  merge(microarray_immune, by = "patient_id", all.x = TRUE) %>%
  dplyr::filter(!is.na(Oligodendrocyte))


