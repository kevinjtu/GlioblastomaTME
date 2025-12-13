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


setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Expression Data")

#Define gene sets####
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

#NOVA geneset
nova_gsea <- list(
  nova = c(  # For overall enrichment
    "CCL19","DCN","FN1","IGFBP4","IGFBP5","IGFBP7","IGHM","IGKC","IL7R",
    "LTB","LUM","MGP","MS4A1","SELL","TCF7", "APOE","AQP4","CCL5","CD8A","CDK6",
    "CXCL10","CXCL11","CXCR6","FCGR3A","FEZ1","IFIT2","IFIT3","NCAM1","NES","OLIG1",
    "OLIG2","PTPRZ1","S100B","SERPINA3","SLC1A3","SOX2","SOX9","SPARCL1","STAT1"
  )
)



# Combine all gene sets into one list
gene_set_list <- c(
  list(REACTOME_STING = sting_genes),
  niche_gsea,
  nova_gsea
)

#TEMPUS RNAseq####
# load the dataset
setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Expression Data")
df <- as.data.frame(fread("TEMPUS RNAseq expression.csv", header = TRUE))
df <- tibble::column_to_rownames(df, var = "name")

#get into correct data type
bulk <- as.data.frame(df)
bulk_ma <- as.matrix(bulk)
if ("TMEM173" %in% rownames(bulk_ma)) {
  rownames(bulk_ma)[rownames(bulk_ma) == "TMEM173"] <- "STING1"
}

params <- ssgseaParam(
  exprData = bulk_ma,
  geneSets = gene_set_list,
  norm     = TRUE       # normalize the scores
)

ssgsea_eset <- as.data.frame(t(gsva(
  params,
  verbose = FALSE
)))
ssgsea_eset$patient_id <- rownames(ssgsea_eset)

setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Sample GSEA")
wb <- createWorkbook()
addWorksheet(wb, "GSEA")
writeData(wb, "GSEA", ssgsea_eset)
saveWorkbook(wb, "TEMPUS_GSEA.xlsx", overwrite = TRUE)


#TGCA RNAseq####
# load the dataset
setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Expression Data")
df <- fread("TCGA RNAseq expression.txt")
df <- df %>%
  filter(!is.na(Hugo_Symbol)) %>%
  distinct(Hugo_Symbol, .keep_all = TRUE) %>%
  dplyr::select(-Entrez_Gene_Id)
df <- tibble::column_to_rownames(df, var = "Hugo_Symbol")

df_tpm <- apply(df, 2, function(x) {
  if (sum(x, na.rm = TRUE) == 0) {
    return(rep(0, length(x)))
  } else {
    return((x / sum(x, na.rm = TRUE)) * 1e6)
  }
})

#get into correct data type
bulk <- as.data.frame(df)
bulk_ma <- as.matrix(bulk)
if ("TMEM173" %in% rownames(bulk_ma)) {
  rownames(bulk_ma)[rownames(bulk_ma) == "TMEM173"] <- "STING1"
}

params <- ssgseaParam(
  exprData = bulk_ma,
  geneSets = gene_set_list,
  norm     = TRUE       # normalize the scores
)

ssgsea_eset <- as.data.frame(t(gsva(
  params,
  verbose = FALSE
)))
ssgsea_eset$sample_id <- rownames(ssgsea_eset)

setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Sample GSEA")
wb <- createWorkbook()
addWorksheet(wb, "GSEA")
writeData(wb, "GSEA", ssgsea_eset)
saveWorkbook(wb, "TCGA_GSEA.xlsx", overwrite = TRUE)

#IVY GAP RNAseq####
# load the dataset
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
df <- tibble::column_to_rownames(df, var = "Hugo_Symbol")

df <- data.frame(lapply(df, function(x) as.numeric(as.character(x))),
                 row.names = rownames(df))

tpm_df <- as.data.frame(apply(df, 2, function(col) {
  (col / sum(col, na.rm = TRUE)) * 1e6
}))
colSums(tpm_df)

#get into correct data type
bulk <- as.data.frame(tpm_df)
bulk_ma <- as.matrix(bulk)
if ("TMEM173" %in% rownames(bulk_ma)) {
  rownames(bulk_ma)[rownames(bulk_ma) == "TMEM173"] <- "STING1"
}

params <- ssgseaParam(
  exprData = bulk_ma,
  geneSets = gene_set_list,
  norm     = TRUE       # normalize the scores
)

ssgsea_eset <- as.data.frame(t(gsva(
  params,
  verbose = FALSE
)))
ssgsea_eset$sample_id <- rownames(ssgsea_eset)

setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Sample GSEA")
wb <- createWorkbook()
addWorksheet(wb, "GSEA")
writeData(wb, "GSEA", ssgsea_eset)
saveWorkbook(wb, "IVYGAP_GSEA.xlsx", overwrite = TRUE)

#Wu 2020 RNAseq####
setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Expression Data")
df <- fread("Wu2020 TPM Expression.txt")
# Convert ENSG IDs in V1 to HUGO symbols
df[, ensg_id := sub("\\.\\d+$", "", V1)]  # Remove version numbers
hugo_symbols <- mapIds(org.Hs.eg.db,
                       keys = df$ensg_id,
                       column = "SYMBOL",
                       keytype = "ENSEMBL",
                       multiVals = "first")
df[, V1 := hugo_symbols][, ensg_id := NULL]  # Update V1 and clean up

df <- df %>%
  filter(!is.na(V1)) %>%
  distinct(V1, .keep_all = TRUE)
df <- tibble::column_to_rownames(df, var = "V1")

#get into correct data type
bulk <- as.data.frame(df)
bulk_ma <- as.matrix(bulk)
if ("TMEM173" %in% rownames(bulk_ma)) {
  rownames(bulk_ma)[rownames(bulk_ma) == "TMEM173"] <- "STING1"
}

params <- ssgseaParam(
  exprData = bulk_ma,
  geneSets = gene_set_list,
  norm     = TRUE       # normalize the scores
)

ssgsea_eset <- as.data.frame(t(gsva(
  params,
  verbose = FALSE
)))
ssgsea_eset$patient_id <- rownames(ssgsea_eset)

setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Sample GSEA")
wb <- createWorkbook()
addWorksheet(wb, "GSEA")
writeData(wb, "GSEA", ssgsea_eset)
saveWorkbook(wb, "Wu2020_GSEA.xlsx", overwrite = TRUE)

#CGGA - batch 2 RNAseq####
setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Expression Data")
df <- fread("CGGA FPKM Expression.txt")

df <- df %>%
  filter(!is.na(Gene_Name)) %>%
  distinct(Gene_Name, .keep_all = TRUE)
df <- tibble::column_to_rownames(df, var = "Gene_Name")

#get into correct data type
bulk <- as.data.frame(df)
bulk_ma <- as.matrix(bulk)
if ("TMEM173" %in% rownames(bulk_ma)) {
  rownames(bulk_ma)[rownames(bulk_ma) == "TMEM173"] <- "STING1"
}

params <- ssgseaParam(
  exprData = bulk_ma,
  geneSets = gene_set_list,
  norm     = TRUE       # normalize the scores
)

ssgsea_eset <- as.data.frame(t(gsva(
  params,
  verbose = FALSE
)))
ssgsea_eset$patient_id <- rownames(ssgsea_eset)

setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Sample GSEA")
wb <- createWorkbook()
addWorksheet(wb, "GSEA")
writeData(wb, "GSEA", ssgsea_eset)
saveWorkbook(wb, "cgga_batch1_GSEA.xlsx", overwrite = TRUE)


#CGGA - batch 1 RNAseq####
setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Expression Data")
df <- fread("CGGA batch1 fpkm.txt")

df <- df %>%
  filter(!is.na(Gene_Name)) %>%
  distinct(Gene_Name, .keep_all = TRUE)
df <- tibble::column_to_rownames(df, var = "Gene_Name")

#get into correct data type
bulk <- as.data.frame(df)
bulk_ma <- as.matrix(bulk)
if ("TMEM173" %in% rownames(bulk_ma)) {
  rownames(bulk_ma)[rownames(bulk_ma) == "TMEM173"] <- "STING1"
}

params <- ssgseaParam(
  exprData = bulk_ma,
  geneSets = gene_set_list,
  norm     = TRUE       # normalize the scores
)

ssgsea_eset <- as.data.frame(t(gsva(
  params,
  verbose = FALSE
)))
ssgsea_eset$patient_id <- rownames(ssgsea_eset)

setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Sample GSEA")
wb <- createWorkbook()
addWorksheet(wb, "GSEA")
writeData(wb, "GSEA", ssgsea_eset)
saveWorkbook(wb, "cgga_batch2_GSEA.xlsx", overwrite = TRUE)


#GLASS RNAseq####
setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Expression Data")
df <- fread("GLASS Expression.tsv")
df <- df %>%
  filter(!is.na(Gene_symbol)) %>%
  distinct(Gene_symbol, .keep_all = TRUE)
df <- tibble::column_to_rownames(df, var = "Gene_symbol")

#get into correct data type
bulk <- as.data.frame(df)
bulk_ma <- as.matrix(bulk)
if ("TMEM173" %in% rownames(bulk_ma)) {
  rownames(bulk_ma)[rownames(bulk_ma) == "TMEM173"] <- "STING1"
}

params <- ssgseaParam(
  exprData = bulk_ma,
  geneSets = gene_set_list,
  norm     = TRUE       # normalize the scores
)

ssgsea_eset <- as.data.frame(t(gsva(
  params,
  verbose = FALSE
)))
ssgsea_eset$sample_id <- rownames(ssgsea_eset)

setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Sample GSEA")
wb <- createWorkbook()
addWorksheet(wb, "GSEA")
writeData(wb, "GSEA", ssgsea_eset)
saveWorkbook(wb, "glass_GSEA.xlsx", overwrite = TRUE)

#G-SAM RNAseq####
setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Expression Data")
df <- fread("gsam raw counts.csv")
df <- df[, V1 := tstrsplit(V1, "\\|")[[2]]]
df <- df %>%
  filter(!is.na(V1)) %>%
  distinct(V1, .keep_all = TRUE)
df <- tibble::column_to_rownames(df, var = "V1")

#get into correct data type
bulk <- as.data.frame(df)
bulk_ma <- as.matrix(bulk)
if ("TMEM173" %in% rownames(bulk_ma)) {
  rownames(bulk_ma)[rownames(bulk_ma) == "TMEM173"] <- "STING1"
}

params <- ssgseaParam(
  exprData = bulk_ma,
  geneSets = gene_set_list,
  norm     = TRUE       # normalize the scores
)

ssgsea_eset <- as.data.frame(t(gsva(
  params,
  verbose = FALSE
)))
ssgsea_eset$sample_id <- rownames(ssgsea_eset)

setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Sample GSEA")
wb <- createWorkbook()
addWorksheet(wb, "GSEA")
writeData(wb, "GSEA", ssgsea_eset)
saveWorkbook(wb, "gsam_GSEA.xlsx", overwrite = TRUE)

#TCGA UG133A####
setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/GlioVis Datasets")
TCGA_GBM <- readRDS("TCGA_GBM.Rds")
df <- TCGA_GBM[["expr"]] %>%
  dplyr::filter(Histology == "GBM") %>%
  dplyr::select(-c(Sample, Histology, Grade, Recurrence, Subtype, CIMP_status, survival, status)) %>%
  t()

#get into correct data type
bulk <- as.data.frame(df)
bulk_ma <- as.matrix(bulk)

if ("TMEM173" %in% rownames(bulk_ma)) {
  rownames(bulk_ma)[rownames(bulk_ma) == "TMEM173"] <- "STING1"
}

params <- ssgseaParam(
  exprData = bulk_ma,
  geneSets = gene_set_list,
  norm     = TRUE       # normalize the scores
)

ssgsea_eset <- as.data.frame(t(gsva(
  params,
  verbose = FALSE
)))
ssgsea_eset$sample_id <- rownames(ssgsea_eset)

setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Sample GSEA")
wb <- createWorkbook()
addWorksheet(wb, "GSEA")
writeData(wb, "GSEA", ssgsea_eset)
saveWorkbook(wb, "tcga_microarray_GSEA.xlsx", overwrite = TRUE)

#Rembrandt####
setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/GlioVis Datasets")
RDS <- readRDS("Rembrandt.Rds")
df <- RDS[["expr"]] %>%
  dplyr::filter(Histology == "GBM") %>%
  dplyr::select(-c(Sample, Histology, Grade, Recurrence, Subtype, CIMP_status, survival, status)) %>%
  t()

#get into correct data type
bulk <- as.data.frame(df)
bulk_ma <- as.matrix(bulk)

if ("TMEM173" %in% rownames(bulk_ma)) {
  rownames(bulk_ma)[rownames(bulk_ma) == "TMEM173"] <- "STING1"
}

params <- ssgseaParam(
  exprData = bulk_ma,
  geneSets = gene_set_list,
  norm     = TRUE       # normalize the scores
)

ssgsea_eset <- as.data.frame(t(gsva(
  params,
  verbose = FALSE
)))
ssgsea_eset$sample_id <- rownames(ssgsea_eset)

setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Sample GSEA")
wb <- createWorkbook()
addWorksheet(wb, "GSEA")
writeData(wb, "GSEA", ssgsea_eset)
saveWorkbook(wb, "Rembrandt_microarray_GSEA.xlsx", overwrite = TRUE)

#Gravendeel ####
setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/GlioVis Datasets")
RDS <- readRDS("Gravendeel.Rds")
df <- RDS[["expr"]] %>%
  dplyr::filter(Histology == "GBM") %>%
  dplyr::select(-c(Sample, Histology, Grade, Recurrence, Subtype, CIMP_status, survival, status)) %>%
  t()

#get into correct data type
bulk <- as.data.frame(df)
bulk_ma <- as.matrix(bulk)

if ("TMEM173" %in% rownames(bulk_ma)) {
  rownames(bulk_ma)[rownames(bulk_ma) == "TMEM173"] <- "STING1"
}

params <- ssgseaParam(
  exprData = bulk_ma,
  geneSets = gene_set_list,
  norm     = TRUE       # normalize the scores
)

ssgsea_eset <- as.data.frame(t(gsva(
  params,
  verbose = FALSE
)))
ssgsea_eset$sample_id <- rownames(ssgsea_eset)

setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Sample GSEA")
wb <- createWorkbook()
addWorksheet(wb, "GSEA")
writeData(wb, "GSEA", ssgsea_eset)
saveWorkbook(wb, "Gravendeel_microarray_GSEA.xlsx", overwrite = TRUE)


#LeeY ####
setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/GlioVis Datasets")
RDS <- readRDS("LeeY.Rds")
df <- RDS[["expr"]] 
colnames(df) <- make.unique(colnames(df))
df <- df %>%
  dplyr::filter(Histology == "GBM") %>%
  dplyr::select(-c(Sample, Histology, Grade, Recurrence, Subtype, CIMP_status, CIMP_status.1, survival, status)) %>%
  t()

#get into correct data type
bulk <- as.data.frame(df)
bulk_ma <- as.matrix(bulk)

if ("TMEM173" %in% rownames(bulk_ma)) {
  rownames(bulk_ma)[rownames(bulk_ma) == "TMEM173"] <- "STING1"
}

params <- ssgseaParam(
  exprData = bulk_ma,
  geneSets = gene_set_list,
  norm     = TRUE       # normalize the scores
)

ssgsea_eset <- as.data.frame(t(gsva(
  params,
  verbose = FALSE
)))
ssgsea_eset$sample_id <- rownames(ssgsea_eset)

setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Sample GSEA")
wb <- createWorkbook()
addWorksheet(wb, "GSEA")
writeData(wb, "GSEA", ssgsea_eset)
saveWorkbook(wb, "LeeY_microarray_GSEA.xlsx", overwrite = TRUE)


#Phillips ####
setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/GlioVis Datasets")
RDS <- readRDS("Phillips.Rds")
df <- RDS[["expr"]] %>%
  dplyr::filter(Histology == "GBM") %>%
  dplyr::select(-c(Sample, Histology, Grade, Recurrence, Subtype, CIMP_status, survival, status)) %>%
  t()

#get into correct data type
bulk <- as.data.frame(df)
bulk_ma <- as.matrix(bulk)

if ("TMEM173" %in% rownames(bulk_ma)) {
  rownames(bulk_ma)[rownames(bulk_ma) == "TMEM173"] <- "STING1"
}

params <- ssgseaParam(
  exprData = bulk_ma,
  geneSets = gene_set_list,
  norm     = TRUE       # normalize the scores
)

ssgsea_eset <- as.data.frame(t(gsva(
  params,
  verbose = FALSE
)))
ssgsea_eset$sample_id <- rownames(ssgsea_eset)

setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Sample GSEA")
wb <- createWorkbook()
addWorksheet(wb, "GSEA")
writeData(wb, "GSEA", ssgsea_eset)
saveWorkbook(wb, "Phillips_microarray_GSEA.xlsx", overwrite = TRUE)



#Freije ####
setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/GlioVis Datasets")
RDS <- readRDS("Freije.Rds")
df <- RDS[["expr"]] %>%
  dplyr::filter(Histology == "GBM") %>%
  dplyr::select(-c(Sample, Histology, Grade, Recurrence, Subtype, CIMP_status, survival, status)) %>%
  t()

#get into correct data type
bulk <- as.data.frame(df)
bulk_ma <- as.matrix(bulk)

if ("TMEM173" %in% rownames(bulk_ma)) {
  rownames(bulk_ma)[rownames(bulk_ma) == "TMEM173"] <- "STING1"
}

params <- ssgseaParam(
  exprData = bulk_ma,
  geneSets = gene_set_list,
  norm     = TRUE       # normalize the scores
)

ssgsea_eset <- as.data.frame(t(gsva(
  params,
  verbose = FALSE
)))
ssgsea_eset$sample_id <- rownames(ssgsea_eset)

setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Sample GSEA")
wb <- createWorkbook()
addWorksheet(wb, "GSEA")
writeData(wb, "GSEA", ssgsea_eset)
saveWorkbook(wb, "Freije_microarray_GSEA.xlsx", overwrite = TRUE)

#Murat ####
setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/GlioVis Datasets")
RDS <- readRDS("Murat.Rds")
df <- RDS[["expr"]] %>%
  dplyr::filter(Histology == "GBM") %>%
  dplyr::select(-c(Sample, Histology, Grade, Recurrence, Subtype, CIMP_status, survival, status)) %>%
  t()

#get into correct data type
bulk <- as.data.frame(df)
bulk_ma <- as.matrix(bulk)

if ("TMEM173" %in% rownames(bulk_ma)) {
  rownames(bulk_ma)[rownames(bulk_ma) == "TMEM173"] <- "STING1"
}

params <- ssgseaParam(
  exprData = bulk_ma,
  geneSets = gene_set_list,
  norm     = TRUE       # normalize the scores
)

ssgsea_eset <- as.data.frame(t(gsva(
  params,
  verbose = FALSE
)))
ssgsea_eset$sample_id <- rownames(ssgsea_eset)

setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Sample GSEA")
wb <- createWorkbook()
addWorksheet(wb, "GSEA")
writeData(wb, "GSEA", ssgsea_eset)
saveWorkbook(wb, "Murat_microarray_GSEA.xlsx", overwrite = TRUE)



#Joo ####
setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/GlioVis Datasets")
RDS <- readRDS("Joo.Rds")
df <- RDS[["expr"]] %>%
  dplyr::filter(Histology == "GBM") %>%
  dplyr::select(-c(Sample, Histology, Grade, Recurrence, Subtype, CIMP_status, survival, status)) %>%
  t()

#get into correct data type
bulk <- as.data.frame(df)
bulk_ma <- as.matrix(bulk)

if ("TMEM173" %in% rownames(bulk_ma)) {
  rownames(bulk_ma)[rownames(bulk_ma) == "TMEM173"] <- "STING1"
}

params <- ssgseaParam(
  exprData = bulk_ma,
  geneSets = gene_set_list,
  norm     = TRUE       # normalize the scores
)

ssgsea_eset <- as.data.frame(t(gsva(
  params,
  verbose = FALSE
)))
ssgsea_eset$sample_id <- rownames(ssgsea_eset)

setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Sample GSEA")
wb <- createWorkbook()
addWorksheet(wb, "GSEA")
writeData(wb, "GSEA", ssgsea_eset)
saveWorkbook(wb, "Joo_microarray_GSEA.xlsx", overwrite = TRUE)



#Ducray ####
setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/GlioVis Datasets")
RDS <- readRDS("Ducray.Rds")
df <- RDS[["expr"]] %>%
  dplyr::filter(Histology == "GBM") %>%
  dplyr::select(-c(Sample, Histology, Grade, Recurrence, Subtype, CIMP_status, survival, status)) %>%
  t()

#get into correct data type
bulk <- as.data.frame(df)
bulk_ma <- as.matrix(bulk)

if ("TMEM173" %in% rownames(bulk_ma)) {
  rownames(bulk_ma)[rownames(bulk_ma) == "TMEM173"] <- "STING1"
}

params <- ssgseaParam(
  exprData = bulk_ma,
  geneSets = gene_set_list,
  norm     = TRUE       # normalize the scores
)

ssgsea_eset <- as.data.frame(t(gsva(
  params,
  verbose = FALSE
)))
ssgsea_eset$sample_id <- rownames(ssgsea_eset)

setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Sample GSEA")
wb <- createWorkbook()
addWorksheet(wb, "GSEA")
writeData(wb, "GSEA", ssgsea_eset)
saveWorkbook(wb, "Ducray_microarray_GSEA.xlsx", overwrite = TRUE)



