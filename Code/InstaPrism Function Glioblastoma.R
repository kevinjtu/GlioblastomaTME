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
library(preprocessCore)

setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Expression Data")

#TEMPUS RNAseq####
# load the dataset
df <- read.csv("TEMPUS RNAseq expression.csv")
df <- column_to_rownames(df, var = "name")

#get into correct data type
bulk <- as.data.frame(df)
bulk_ma <- as.matrix(bulk)
bulk_ar <- array(bulk_ma, dim = c(nrow(bulk_ma), ncol(bulk_ma)))
rownames(bulk_ar) <- rownames(bulk_ma)
colnames(bulk_ar) <- colnames(bulk_ma)
class(bulk_ar)

# take GBM_refPhi for reference
setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME")
GBM_refPhi <- readRDS("GBM_refPhi.RDS")
length(intersect(rownames(bulk_ar), rownames(GBM_refPhi@phi.cs)))

InstaPrism.res <- InstaPrism(input_type = 'refPhi_cs', filter=FALSE, bulk_Expr = bulk_ar, refPhi_cs = GBM_refPhi, n.core = 16)

setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Instaprism Results")
wb <- createWorkbook()
addWorksheet(wb, "Cell Types")
writeData(wb, "Cell Types", InstaPrism.res@Post.ini.ct@theta, rowNames = TRUE)
addWorksheet(wb, "Cell States")
writeData(wb, "Cell States", InstaPrism.res@Post.ini.cs@theta, rowNames = TRUE)
addWorksheet(wb, "Cell Types Transposed")
writeData(wb, "Cell Types Transposed", t(InstaPrism.res@Post.ini.ct@theta), rowNames = TRUE)
addWorksheet(wb, "Cell States Transposed")
writeData(wb, "Cell States Transposed", t(InstaPrism.res@Post.ini.cs@theta), rowNames = TRUE)

file_name <- "TEMPUS RNAseq InstaPrism results.xlsx"
saveWorkbook(wb, file = file_name, overwrite = TRUE)

#TGCA RNAseq####
# load the dataset
df <- fread("TCGA RNAseq expression.txt")
df <- df %>%
  filter(!is.na(Hugo_Symbol)) %>%
  distinct(Hugo_Symbol, .keep_all = TRUE) %>%
  dplyr::select(-Entrez_Gene_Id)
df <- column_to_rownames(df, var = "Hugo_Symbol")

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
bulk_ar <- array(bulk_ma, dim = c(nrow(bulk_ma), ncol(bulk_ma)))
rownames(bulk_ar) <- rownames(bulk_ma)
colnames(bulk_ar) <- colnames(bulk_ma)
class(bulk_ar)

# take GBM_refPhi for reference
setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME")
GBM_refPhi <- readRDS("GBM_refPhi.RDS")
length(intersect(rownames(bulk_ar), rownames(GBM_refPhi@phi.cs)))

InstaPrism.res <- InstaPrism(input_type = 'refPhi_cs', filter=FALSE, bulk_Expr = bulk_ar, refPhi_cs = GBM_refPhi, n.core = 16)

setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Instaprism Results")
wb <- createWorkbook()
addWorksheet(wb, "Cell Types")
writeData(wb, "Cell Types", InstaPrism.res@Post.ini.ct@theta, rowNames = TRUE)
addWorksheet(wb, "Cell States")
writeData(wb, "Cell States", InstaPrism.res@Post.ini.cs@theta, rowNames = TRUE)
addWorksheet(wb, "Cell Types Transposed")
writeData(wb, "Cell Types Transposed", t(InstaPrism.res@Post.ini.ct@theta), rowNames = TRUE)
addWorksheet(wb, "Cell States Transposed")
writeData(wb, "Cell States Transposed", t(InstaPrism.res@Post.ini.cs@theta), rowNames = TRUE)

file_name <- "TCGA RNAseq InstaPrism results.xlsx"
saveWorkbook(wb, file = file_name, overwrite = TRUE)

#IVY GAP RNAseq####
# load the dataset
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

#get into correct data type
bulk <- as.data.frame(tpm_df)
bulk_ma <- as.matrix(bulk)
bulk_ar <- array(bulk_ma, dim = c(nrow(bulk_ma), ncol(bulk_ma)))
rownames(bulk_ar) <- rownames(bulk_ma)
colnames(bulk_ar) <- colnames(bulk_ma)
class(bulk_ar)

# take GBM_refPhi for reference
setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME")
GBM_refPhi <- readRDS("GBM_refPhi.RDS")
length(intersect(rownames(bulk_ar), rownames(GBM_refPhi@phi.cs)))

InstaPrism.res <- InstaPrism(input_type = 'refPhi_cs', filter=FALSE, bulk_Expr = bulk_ar, refPhi_cs = GBM_refPhi, n.core = 16)

setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Instaprism Results")
wb <- createWorkbook()
addWorksheet(wb, "Cell Types")
writeData(wb, "Cell Types", InstaPrism.res@Post.ini.ct@theta, rowNames = TRUE)
addWorksheet(wb, "Cell States")
writeData(wb, "Cell States", InstaPrism.res@Post.ini.cs@theta, rowNames = TRUE)
addWorksheet(wb, "Cell Types Transposed")
writeData(wb, "Cell Types Transposed", t(InstaPrism.res@Post.ini.ct@theta), rowNames = TRUE)
addWorksheet(wb, "Cell States Transposed")
writeData(wb, "Cell States Transposed", t(InstaPrism.res@Post.ini.cs@theta), rowNames = TRUE)

file_name <- "IVYGAP RNAseq InstaPrism results.xlsx"
saveWorkbook(wb, file = file_name, overwrite = TRUE)

#Wu 2020 RNAseq####
df <- fread("Wu2020 Expression.txt")
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
df <- column_to_rownames(df, var = "V1")

#get into correct data type
bulk <- as.data.frame(df)
bulk_ma <- as.matrix(bulk)
bulk_ar <- array(bulk_ma, dim = c(nrow(bulk_ma), ncol(bulk_ma)))
rownames(bulk_ar) <- rownames(bulk_ma)
colnames(bulk_ar) <- colnames(bulk_ma)
class(bulk_ar)

# take GBM_refPhi for reference
setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME")
GBM_refPhi <- readRDS("GBM_refPhi.RDS")
length(intersect(rownames(bulk_ar), rownames(GBM_refPhi@phi.cs)))

InstaPrism.res <- InstaPrism(input_type = 'refPhi_cs', filter=FALSE, bulk_Expr = bulk_ar, refPhi_cs = GBM_refPhi, n.core = 16)

setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Instaprism Results")
wb <- createWorkbook()
addWorksheet(wb, "Cell Types")
writeData(wb, "Cell Types", InstaPrism.res@Post.ini.ct@theta, rowNames = TRUE)
addWorksheet(wb, "Cell States")
writeData(wb, "Cell States", InstaPrism.res@Post.ini.cs@theta, rowNames = TRUE)
addWorksheet(wb, "Cell Types Transposed")
writeData(wb, "Cell Types Transposed", t(InstaPrism.res@Post.ini.ct@theta), rowNames = TRUE)
addWorksheet(wb, "Cell States Transposed")
writeData(wb, "Cell States Transposed", t(InstaPrism.res@Post.ini.cs@theta), rowNames = TRUE)

file_name <- "Wu2020 RNAseq InstaPrism results.xlsx"
saveWorkbook(wb, file = file_name, overwrite = TRUE)

#Huang 2021 RNAseq####
# Specify the directory containing the files
directory <- "~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Expression Data/Huang 2021 Expression"

# Get list of all GSM files in the directory
file_list <- list.files(
  path = directory,
  pattern = "^GSM\\d+_.*\\.txt\\.gz$",
  full.names = TRUE
)

# Check if files were found
if (length(file_list) == 0) {
  stop("No matching files found in the directory.")
}

# Function to extract SYC ID from filename
extract_syc_id <- function(filename) {
  # Split filename into parts and get last segment
  parts <- strsplit(basename(filename), "_")[[1]]
  last_part <- parts[length(parts)]
  # Remove .txt.gz extension
  sub("\\.txt\\.gz$", "", last_part)
}

# Initialize with first file
first_file <- file_list[1]
patient_id <- extract_syc_id(first_file)
first_data <- read.table(
  first_file,
  sep = "\t",
  header = FALSE,
  col.names = c("ENSG", patient_id)
)
combined_df <- first_data

# Process remaining files
for (file in file_list[-1]) {
  current_patient <- extract_syc_id(file)
  current_data <- read.table(
    file,
    sep = "\t",
    header = FALSE,
    col.names = c("ENSG", current_patient)
  )
  
  # Verify ENSG IDs match
  if (!identical(combined_df$ENSG, current_data$ENSG)) {
    stop(paste("ENSG mismatch in file:", file))
  }
  
  combined_df[[current_patient]] <- current_data[[current_patient]]
}

# Convert to matrix (if needed)
gene_ids <- combined_df$ENSG
count_matrix <- as.matrix(combined_df[, -1])
rownames(count_matrix) <- gene_ids

# View the resulting matrix (first 6 rows)
df <- combined_df
head(df)

# Clean ENSG IDs by removing version numbers
df$ENSG <- sub("\\.\\d+$", "", df$ENSG)

# Map to HUGO symbols
hugo_symbols <- mapIds(org.Hs.eg.db,
                       keys = df$ENSG,
                       column = "SYMBOL",
                       keytype = "ENSEMBL",
                       multiVals = "first")
df$ENSG <- hugo_symbols
colnames(df)[1] <- "SYMBOL"

df <- df %>%
  filter(!is.na(SYMBOL)) %>%
  distinct(SYMBOL, .keep_all = TRUE)
df <- column_to_rownames(df, var = "SYMBOL")


#get into correct data type
bulk <- as.data.frame(df)
bulk_ma <- as.matrix(bulk)
bulk_ar <- array(bulk_ma, dim = c(nrow(bulk_ma), ncol(bulk_ma)))
rownames(bulk_ar) <- rownames(bulk_ma)
colnames(bulk_ar) <- colnames(bulk_ma)
class(bulk_ar)

# take GBM_refPhi for reference
setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME")
GBM_refPhi <- readRDS("GBM_refPhi.RDS")
length(intersect(rownames(bulk_ar), rownames(GBM_refPhi@phi.cs)))

InstaPrism.res <- InstaPrism(input_type = 'refPhi_cs', filter=FALSE, bulk_Expr = bulk_ar, refPhi_cs = GBM_refPhi, n.core = 16)

setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Instaprism Results")
wb <- createWorkbook()
addWorksheet(wb, "Cell Types")
writeData(wb, "Cell Types", InstaPrism.res@Post.ini.ct@theta, rowNames = TRUE)
addWorksheet(wb, "Cell States")
writeData(wb, "Cell States", InstaPrism.res@Post.ini.cs@theta, rowNames = TRUE)
addWorksheet(wb, "Cell Types Transposed")
writeData(wb, "Cell Types Transposed", t(InstaPrism.res@Post.ini.ct@theta), rowNames = TRUE)
addWorksheet(wb, "Cell States Transposed")
writeData(wb, "Cell States Transposed", t(InstaPrism.res@Post.ini.cs@theta), rowNames = TRUE)

file_name <- "Huang2021 RNAseq InstaPrism results.xlsx"
saveWorkbook(wb, file = file_name, overwrite = TRUE)

#CGGA - batch 2 RNAseq####
df <- fread("CGGA Expression.txt")

df <- df %>%
  filter(!is.na(gene_name)) %>%
  distinct(gene_name, .keep_all = TRUE)
df <- column_to_rownames(df, var = "gene_name")

#get into correct data type
bulk <- as.data.frame(df)
bulk_ma <- as.matrix(bulk)
bulk_ar <- array(bulk_ma, dim = c(nrow(bulk_ma), ncol(bulk_ma)))
rownames(bulk_ar) <- rownames(bulk_ma)
colnames(bulk_ar) <- colnames(bulk_ma)
class(bulk_ar)

# take GBM_refPhi for reference
setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME")
GBM_refPhi <- readRDS("GBM_refPhi.RDS")
length(intersect(rownames(bulk_ar), rownames(GBM_refPhi@phi.cs)))

InstaPrism.res <- InstaPrism(input_type = 'refPhi_cs', filter=FALSE, bulk_Expr = bulk_ar, refPhi_cs = GBM_refPhi, n.core = 16)

setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Instaprism Results")
wb <- createWorkbook()
addWorksheet(wb, "Cell Types")
writeData(wb, "Cell Types", InstaPrism.res@Post.ini.ct@theta, rowNames = TRUE)
addWorksheet(wb, "Cell States")
writeData(wb, "Cell States", InstaPrism.res@Post.ini.cs@theta, rowNames = TRUE)
addWorksheet(wb, "Cell Types Transposed")
writeData(wb, "Cell Types Transposed", t(InstaPrism.res@Post.ini.ct@theta), rowNames = TRUE)
addWorksheet(wb, "Cell States Transposed")
writeData(wb, "Cell States Transposed", t(InstaPrism.res@Post.ini.cs@theta), rowNames = TRUE)

file_name <- "CGGA RNAseq InstaPrism results.xlsx"
saveWorkbook(wb, file = file_name, overwrite = TRUE)


#CGGA - batch 1 RNAseq####
df <- fread("CGGA batch1 raw counts.txt")

df <- df %>%
  filter(!is.na(gene_name)) %>%
  distinct(gene_name, .keep_all = TRUE)
df <- column_to_rownames(df, var = "gene_name")

#get into correct data type
bulk <- as.data.frame(df)
bulk_ma <- as.matrix(bulk)
bulk_ar <- array(bulk_ma, dim = c(nrow(bulk_ma), ncol(bulk_ma)))
rownames(bulk_ar) <- rownames(bulk_ma)
colnames(bulk_ar) <- colnames(bulk_ma)
class(bulk_ar)

# take GBM_refPhi for reference
setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME")
GBM_refPhi <- readRDS("GBM_refPhi.RDS")
length(intersect(rownames(bulk_ar), rownames(GBM_refPhi@phi.cs)))

InstaPrism.res <- InstaPrism(input_type = 'refPhi_cs', filter=FALSE, bulk_Expr = bulk_ar, refPhi_cs = GBM_refPhi, n.core = 16)

setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Instaprism Results")
wb <- createWorkbook()
addWorksheet(wb, "Cell Types")
writeData(wb, "Cell Types", InstaPrism.res@Post.ini.ct@theta, rowNames = TRUE)
addWorksheet(wb, "Cell States")
writeData(wb, "Cell States", InstaPrism.res@Post.ini.cs@theta, rowNames = TRUE)
addWorksheet(wb, "Cell Types Transposed")
writeData(wb, "Cell Types Transposed", t(InstaPrism.res@Post.ini.ct@theta), rowNames = TRUE)
addWorksheet(wb, "Cell States Transposed")
writeData(wb, "Cell States Transposed", t(InstaPrism.res@Post.ini.cs@theta), rowNames = TRUE)

file_name <- "CGGA Batch 1 RNAseq InstaPrism results.xlsx"
saveWorkbook(wb, file = file_name, overwrite = TRUE)


#Gill 2014 - Core RNAseq####
df <- fread("Gill2014 Core Biopsy FPKM.txt")
df <- df %>%
  filter(!is.na(V1)) %>%
  distinct(V1, .keep_all = TRUE)
df <- column_to_rownames(df, var = "V1")

df <- data.frame(lapply(df, function(x) as.numeric(as.character(x))),
                 row.names = rownames(df))

tpm_df <- as.data.frame(apply(df, 2, function(col) {
  (col / sum(col, na.rm = TRUE)) * 1e6
}))

#get into correct data type
bulk <- as.data.frame(tpm_df)
bulk_ma <- as.matrix(bulk)
bulk_ar <- array(bulk_ma, dim = c(nrow(bulk_ma), ncol(bulk_ma)))
rownames(bulk_ar) <- rownames(bulk_ma)
colnames(bulk_ar) <- colnames(bulk_ma)
class(bulk_ar)

# take GBM_refPhi for reference
setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME")
GBM_refPhi <- readRDS("GBM_refPhi.RDS")
length(intersect(rownames(bulk_ar), rownames(GBM_refPhi@phi.cs)))

InstaPrism.res <- InstaPrism(input_type = 'refPhi_cs', filter=FALSE, bulk_Expr = bulk_ar, refPhi_cs = GBM_refPhi, n.core = 16)

setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Instaprism Results")
wb <- createWorkbook()
addWorksheet(wb, "Cell Types")
writeData(wb, "Cell Types", InstaPrism.res@Post.ini.ct@theta, rowNames = TRUE)
addWorksheet(wb, "Cell States")
writeData(wb, "Cell States", InstaPrism.res@Post.ini.cs@theta, rowNames = TRUE)
addWorksheet(wb, "Cell Types Transposed")
writeData(wb, "Cell Types Transposed", t(InstaPrism.res@Post.ini.ct@theta), rowNames = TRUE)
addWorksheet(wb, "Cell States Transposed")
writeData(wb, "Cell States Transposed", t(InstaPrism.res@Post.ini.cs@theta), rowNames = TRUE)

file_name <- "Gill2014Core RNAseq InstaPrism results.xlsx"
saveWorkbook(wb, file = file_name, overwrite = TRUE)

#Gill 2014 - Margin RNAseq####
df <- fread("Gill2014 Margin Biopsy FPKM.txt")
df <- df %>%
  filter(!is.na(V1)) %>%
  distinct(V1, .keep_all = TRUE)
df <- column_to_rownames(df, var = "V1")

df <- data.frame(lapply(df, function(x) as.numeric(as.character(x))),
                 row.names = rownames(df))

tpm_df <- as.data.frame(apply(df, 2, function(col) {
  (col / sum(col, na.rm = TRUE)) * 1e6
}))

#get into correct data type
bulk <- as.data.frame(tpm_df)
bulk_ma <- as.matrix(bulk)
bulk_ar <- array(bulk_ma, dim = c(nrow(bulk_ma), ncol(bulk_ma)))
rownames(bulk_ar) <- rownames(bulk_ma)
colnames(bulk_ar) <- colnames(bulk_ma)
class(bulk_ar)

# take GBM_refPhi for reference
setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME")
GBM_refPhi <- readRDS("GBM_refPhi.RDS")
length(intersect(rownames(bulk_ar), rownames(GBM_refPhi@phi.cs)))

InstaPrism.res <- InstaPrism(input_type = 'refPhi_cs', filter=FALSE, bulk_Expr = bulk_ar, refPhi_cs = GBM_refPhi, n.core = 16)

setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Instaprism Results")
wb <- createWorkbook()
addWorksheet(wb, "Cell Types")
writeData(wb, "Cell Types", InstaPrism.res@Post.ini.ct@theta, rowNames = TRUE)
addWorksheet(wb, "Cell States")
writeData(wb, "Cell States", InstaPrism.res@Post.ini.cs@theta, rowNames = TRUE)
addWorksheet(wb, "Cell Types Transposed")
writeData(wb, "Cell Types Transposed", t(InstaPrism.res@Post.ini.ct@theta), rowNames = TRUE)
addWorksheet(wb, "Cell States Transposed")
writeData(wb, "Cell States Transposed", t(InstaPrism.res@Post.ini.cs@theta), rowNames = TRUE)

file_name <- "Gill2014Margin RNAseq InstaPrism results.xlsx"
saveWorkbook(wb, file = file_name, overwrite = TRUE)
#Gill 2014 - Noncancer RNAseq####
df <- fread("Gill2014 Nonneoplastic Biopsy FPKM.txt")
df <- df %>%
  filter(!is.na(V1)) %>%
  distinct(V1, .keep_all = TRUE)
df <- column_to_rownames(df, var = "V1")

df <- data.frame(lapply(df, function(x) as.numeric(as.character(x))),
                 row.names = rownames(df))

tpm_df <- as.data.frame(apply(df, 2, function(col) {
  (col / sum(col, na.rm = TRUE)) * 1e6
}))

#get into correct data type
bulk <- as.data.frame(tpm_df)
bulk_ma <- as.matrix(bulk)
bulk_ar <- array(bulk_ma, dim = c(nrow(bulk_ma), ncol(bulk_ma)))
rownames(bulk_ar) <- rownames(bulk_ma)
colnames(bulk_ar) <- colnames(bulk_ma)
class(bulk_ar)

# take GBM_refPhi for reference
setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME")
GBM_refPhi <- readRDS("GBM_refPhi.RDS")
length(intersect(rownames(bulk_ar), rownames(GBM_refPhi@phi.cs)))

InstaPrism.res <- InstaPrism(input_type = 'refPhi_cs', filter=FALSE, bulk_Expr = bulk_ar, refPhi_cs = GBM_refPhi, n.core = 16)

setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Instaprism Results")
wb <- createWorkbook()
addWorksheet(wb, "Cell Types")
writeData(wb, "Cell Types", InstaPrism.res@Post.ini.ct@theta, rowNames = TRUE)
addWorksheet(wb, "Cell States")
writeData(wb, "Cell States", InstaPrism.res@Post.ini.cs@theta, rowNames = TRUE)
addWorksheet(wb, "Cell Types Transposed")
writeData(wb, "Cell Types Transposed", t(InstaPrism.res@Post.ini.ct@theta), rowNames = TRUE)
addWorksheet(wb, "Cell States Transposed")
writeData(wb, "Cell States Transposed", t(InstaPrism.res@Post.ini.cs@theta), rowNames = TRUE)

file_name <- "Gill2014Noncancer RNAseq InstaPrism results.xlsx"
saveWorkbook(wb, file = file_name, overwrite = TRUE)
#GLASS RNAseq####
df <- fread("GLASS Expression.tsv")
df <- df %>%
  filter(!is.na(Gene_symbol)) %>%
  distinct(Gene_symbol, .keep_all = TRUE)
df <- column_to_rownames(df, var = "Gene_symbol")

#get into correct data type
bulk <- as.data.frame(df)
bulk_ma <- as.matrix(bulk)
bulk_ar <- array(bulk_ma, dim = c(nrow(bulk_ma), ncol(bulk_ma)))
rownames(bulk_ar) <- rownames(bulk_ma)
colnames(bulk_ar) <- colnames(bulk_ma)
class(bulk_ar)

# take GBM_refPhi for reference
setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME")
GBM_refPhi <- readRDS("GBM_refPhi.RDS")
length(intersect(rownames(bulk_ar), rownames(GBM_refPhi@phi.cs)))

InstaPrism.res <- InstaPrism(input_type = 'refPhi_cs', filter=FALSE, bulk_Expr = bulk_ar, refPhi_cs = GBM_refPhi, n.core = 16)

setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Instaprism Results")
wb <- createWorkbook()
addWorksheet(wb, "Cell Types")
writeData(wb, "Cell Types", InstaPrism.res@Post.ini.ct@theta, rowNames = TRUE)
addWorksheet(wb, "Cell States")
writeData(wb, "Cell States", InstaPrism.res@Post.ini.cs@theta, rowNames = TRUE)
addWorksheet(wb, "Cell Types Transposed")
writeData(wb, "Cell Types Transposed", t(InstaPrism.res@Post.ini.ct@theta), rowNames = TRUE)
addWorksheet(wb, "Cell States Transposed")
writeData(wb, "Cell States Transposed", t(InstaPrism.res@Post.ini.cs@theta), rowNames = TRUE)

file_name <- "GLASS RNAseq InstaPrism results.xlsx"
saveWorkbook(wb, file = file_name, overwrite = TRUE)

#G-SAM RNAseq####
df <- fread("gsam raw counts.csv")
df <- df[, V1 := tstrsplit(V1, "\\|")[[2]]]
df <- df %>%
  filter(!is.na(V1)) %>%
  distinct(V1, .keep_all = TRUE)
df <- column_to_rownames(df, var = "V1")

#get into correct data type
bulk <- as.data.frame(df)
bulk_ma <- as.matrix(bulk)
bulk_ar <- array(bulk_ma, dim = c(nrow(bulk_ma), ncol(bulk_ma)))
rownames(bulk_ar) <- rownames(bulk_ma)
colnames(bulk_ar) <- colnames(bulk_ma)
class(bulk_ar)

# take GBM_refPhi for reference
setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME")
GBM_refPhi <- readRDS("GBM_refPhi.RDS")
length(intersect(rownames(bulk_ar), rownames(GBM_refPhi@phi.cs)))

InstaPrism.res <- InstaPrism(input_type = 'refPhi_cs', filter=FALSE, bulk_Expr = bulk_ar, refPhi_cs = GBM_refPhi, n.core = 16)

setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Instaprism Results")
wb <- createWorkbook()
addWorksheet(wb, "Cell Types")
writeData(wb, "Cell Types", InstaPrism.res@Post.ini.ct@theta, rowNames = TRUE)
addWorksheet(wb, "Cell States")
writeData(wb, "Cell States", InstaPrism.res@Post.ini.cs@theta, rowNames = TRUE)
addWorksheet(wb, "Cell Types Transposed")
writeData(wb, "Cell Types Transposed", t(InstaPrism.res@Post.ini.ct@theta), rowNames = TRUE)
addWorksheet(wb, "Cell States Transposed")
writeData(wb, "Cell States Transposed", t(InstaPrism.res@Post.ini.cs@theta), rowNames = TRUE)

file_name <- "GSAM RNAseq InstaPrism results.xlsx"
saveWorkbook(wb, file = file_name, overwrite = TRUE)

#Neoadjuvant1 RNAseq####
setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Expression Data")
df <- read_excel("NeoadjuvantStage1_counts.xlsx") %>%
  distinct(Genes, .keep_all = TRUE) %>%
  column_to_rownames(var = "Genes")

#get into correct data type
bulk <- as.data.frame(df)
bulk_ma <- as.matrix(bulk)
bulk_ar <- array(bulk_ma, dim = c(nrow(bulk_ma), ncol(bulk_ma)))
rownames(bulk_ar) <- rownames(bulk_ma)
colnames(bulk_ar) <- colnames(bulk_ma)
class(bulk_ar)

# take GBM_refPhi for reference
setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME")
GBM_refPhi <- readRDS("GBM_refPhi.RDS")
length(intersect(rownames(bulk_ar), rownames(GBM_refPhi@phi.cs)))

InstaPrism.res <- InstaPrism(input_type = 'refPhi_cs', filter=FALSE, bulk_Expr = bulk_ar, refPhi_cs = GBM_refPhi, n.core = 16)

setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Instaprism Results")
wb <- createWorkbook()
addWorksheet(wb, "Cell Types")
writeData(wb, "Cell Types", InstaPrism.res@Post.ini.ct@theta, rowNames = TRUE)
addWorksheet(wb, "Cell States")
writeData(wb, "Cell States", InstaPrism.res@Post.ini.cs@theta, rowNames = TRUE)
addWorksheet(wb, "Cell Types Transposed")
writeData(wb, "Cell Types Transposed", t(InstaPrism.res@Post.ini.ct@theta), rowNames = TRUE)
addWorksheet(wb, "Cell States Transposed")
writeData(wb, "Cell States Transposed", t(InstaPrism.res@Post.ini.cs@theta), rowNames = TRUE)

file_name <- "NeoadjuvantStage1 RNAseq InstaPrism results.xlsx"
saveWorkbook(wb, file = file_name, overwrite = TRUE)

#Neoadjuvant2 RNAseq####
setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Expression Data")
df <- read_excel("NeoadjuvantStage2_counts.xlsx") %>%
  distinct(Gene, .keep_all = TRUE) %>%
  column_to_rownames(var = "Gene")

#get into correct data type
bulk <- as.data.frame(df)
bulk_ma <- as.matrix(bulk)
bulk_ar <- array(bulk_ma, dim = c(nrow(bulk_ma), ncol(bulk_ma)))
rownames(bulk_ar) <- rownames(bulk_ma)
colnames(bulk_ar) <- colnames(bulk_ma)
class(bulk_ar)

# take GBM_refPhi for reference
setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME")
GBM_refPhi <- readRDS("GBM_refPhi.RDS")
length(intersect(rownames(bulk_ar), rownames(GBM_refPhi@phi.cs)))

InstaPrism.res <- InstaPrism(input_type = 'refPhi_cs', filter=FALSE, bulk_Expr = bulk_ar, refPhi_cs = GBM_refPhi, n.core = 16)

setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Instaprism Results")
wb <- createWorkbook()
addWorksheet(wb, "Cell Types")
writeData(wb, "Cell Types", InstaPrism.res@Post.ini.ct@theta, rowNames = TRUE)
addWorksheet(wb, "Cell States")
writeData(wb, "Cell States", InstaPrism.res@Post.ini.cs@theta, rowNames = TRUE)
addWorksheet(wb, "Cell Types Transposed")
writeData(wb, "Cell Types Transposed", t(InstaPrism.res@Post.ini.ct@theta), rowNames = TRUE)
addWorksheet(wb, "Cell States Transposed")
writeData(wb, "Cell States Transposed", t(InstaPrism.res@Post.ini.cs@theta), rowNames = TRUE)

file_name <- "NeoadjuvantStage2 RNAseq InstaPrism results.xlsx"
saveWorkbook(wb, file = file_name, overwrite = TRUE)

#TCGA UG133A####
setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/GlioVis Datasets")
TCGA_GBM <- readRDS("TCGA_GBM.Rds")
df <- TCGA_GBM[["expr"]] %>%
  dplyr::filter(Histology == "GBM") %>%
  dplyr::select(-c(Sample, Histology, Grade, Recurrence, Subtype, CIMP_status, survival, status)) %>%
  t()
  
df_linear <- 2^df  # genes in rows, samples in columns
df_tpm <- sweep(df_linear, 2, colSums(df_linear), FUN = "/") * 1e6
df_tpm_matrix <- t(df_tpm)  # transpose: samples as rows
df_tpm_qn <- normalize.quantiles(df_tpm_matrix)
rownames(df_tpm_qn) <- rownames(df_tpm_matrix)
colnames(df_tpm_qn) <- colnames(df_tpm_matrix)
df_tpm_qn_final <- t(df_tpm_qn)

#get into correct data type
bulk <- as.data.frame(df_tpm_qn_final)
bulk_ma <- as.matrix(bulk)
bulk_ar <- array(bulk_ma, dim = c(nrow(bulk_ma), ncol(bulk_ma)))
rownames(bulk_ar) <- rownames(bulk_ma)
colnames(bulk_ar) <- colnames(bulk_ma)
class(bulk_ar)

# take GBM_refPhi for reference
setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME")
GBM_refPhi <- readRDS("GBM_refPhi.RDS")
length(intersect(rownames(bulk_ar), rownames(GBM_refPhi@phi.cs)))

InstaPrism.res <- InstaPrism(input_type = 'refPhi_cs', filter=FALSE, bulk_Expr = bulk_ar, refPhi_cs = GBM_refPhi, n.core = 16)

setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Instaprism Results")
wb <- createWorkbook()
addWorksheet(wb, "Cell Types")
writeData(wb, "Cell Types", InstaPrism.res@Post.ini.ct@theta, rowNames = TRUE)
addWorksheet(wb, "Cell States")
writeData(wb, "Cell States", InstaPrism.res@Post.ini.cs@theta, rowNames = TRUE)
addWorksheet(wb, "Cell Types Transposed")
writeData(wb, "Cell Types Transposed", t(InstaPrism.res@Post.ini.ct@theta), rowNames = TRUE)
addWorksheet(wb, "Cell States Transposed")
writeData(wb, "Cell States Transposed", t(InstaPrism.res@Post.ini.cs@theta), rowNames = TRUE)

file_name <- "TCGA Microarray InstaPrism results.xlsx"
saveWorkbook(wb, file = file_name, overwrite = TRUE)


#Rembrandt####
setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/GlioVis Datasets")
Rembrandt <- readRDS("Rembrandt.Rds")
df <- Rembrandt[["expr"]] %>%
  dplyr::filter(Histology == "GBM") %>%
  dplyr::select(-c(Sample, Histology, Grade, Recurrence, Subtype, CIMP_status, survival, status)) %>%
  t()

df_linear <- 2^df  # genes in rows, samples in columns
df_tpm <- sweep(df_linear, 2, colSums(df_linear), FUN = "/") * 1e6
df_tpm_matrix <- t(df_tpm)  # transpose: samples as rows
df_tpm_qn <- normalize.quantiles(df_tpm_matrix)
rownames(df_tpm_qn) <- rownames(df_tpm_matrix)
colnames(df_tpm_qn) <- colnames(df_tpm_matrix)
df_tpm_qn_final <- t(df_tpm_qn)

#get into correct data type
bulk <- as.data.frame(df_tpm_qn_final)
bulk_ma <- as.matrix(bulk)
bulk_ar <- array(bulk_ma, dim = c(nrow(bulk_ma), ncol(bulk_ma)))
rownames(bulk_ar) <- rownames(bulk_ma)
colnames(bulk_ar) <- colnames(bulk_ma)
class(bulk_ar)

# take GBM_refPhi for reference
setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME")
GBM_refPhi <- readRDS("GBM_refPhi.RDS")
length(intersect(rownames(bulk_ar), rownames(GBM_refPhi@phi.cs)))

InstaPrism.res <- InstaPrism(input_type = 'refPhi_cs', filter=FALSE, bulk_Expr = bulk_ar, refPhi_cs = GBM_refPhi, n.core = 16)

setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Instaprism Results")
wb <- createWorkbook()
addWorksheet(wb, "Cell Types")
writeData(wb, "Cell Types", InstaPrism.res@Post.ini.ct@theta, rowNames = TRUE)
addWorksheet(wb, "Cell States")
writeData(wb, "Cell States", InstaPrism.res@Post.ini.cs@theta, rowNames = TRUE)
addWorksheet(wb, "Cell Types Transposed")
writeData(wb, "Cell Types Transposed", t(InstaPrism.res@Post.ini.ct@theta), rowNames = TRUE)
addWorksheet(wb, "Cell States Transposed")
writeData(wb, "Cell States Transposed", t(InstaPrism.res@Post.ini.cs@theta), rowNames = TRUE)

file_name <- "Rembrandt Microarray InstaPrism results.xlsx"
saveWorkbook(wb, file = file_name, overwrite = TRUE)


#Gravendeel ####
setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/GlioVis Datasets")
Gravendeel <- readRDS("Gravendeel.Rds")
df <- Gravendeel[["expr"]] %>%
  dplyr::filter(Histology == "GBM") %>%
  dplyr::select(-c(Sample, Histology, Grade, Recurrence, Subtype, CIMP_status, survival, status)) %>%
  t()

df_linear <- 2^df  # genes in rows, samples in columns
df_tpm <- sweep(df_linear, 2, colSums(df_linear), FUN = "/") * 1e6
df_tpm_matrix <- t(df_tpm)  # transpose: samples as rows
df_tpm_qn <- normalize.quantiles(df_tpm_matrix)
rownames(df_tpm_qn) <- rownames(df_tpm_matrix)
colnames(df_tpm_qn) <- colnames(df_tpm_matrix)
df_tpm_qn_final <- t(df_tpm_qn)

#get into correct data type
bulk <- as.data.frame(df_tpm_qn_final)
bulk_ma <- as.matrix(bulk)
bulk_ar <- array(bulk_ma, dim = c(nrow(bulk_ma), ncol(bulk_ma)))
rownames(bulk_ar) <- rownames(bulk_ma)
colnames(bulk_ar) <- colnames(bulk_ma)
class(bulk_ar)

# take GBM_refPhi for reference
setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME")
GBM_refPhi <- readRDS("GBM_refPhi.RDS")
length(intersect(rownames(bulk_ar), rownames(GBM_refPhi@phi.cs)))

InstaPrism.res <- InstaPrism(input_type = 'refPhi_cs', filter=FALSE, bulk_Expr = bulk_ar, refPhi_cs = GBM_refPhi, n.core = 16)

setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Instaprism Results")
wb <- createWorkbook()
addWorksheet(wb, "Cell Types")
writeData(wb, "Cell Types", InstaPrism.res@Post.ini.ct@theta, rowNames = TRUE)
addWorksheet(wb, "Cell States")
writeData(wb, "Cell States", InstaPrism.res@Post.ini.cs@theta, rowNames = TRUE)
addWorksheet(wb, "Cell Types Transposed")
writeData(wb, "Cell Types Transposed", t(InstaPrism.res@Post.ini.ct@theta), rowNames = TRUE)
addWorksheet(wb, "Cell States Transposed")
writeData(wb, "Cell States Transposed", t(InstaPrism.res@Post.ini.cs@theta), rowNames = TRUE)

file_name <- "Gravendeel Microarray InstaPrism results.xlsx"
saveWorkbook(wb, file = file_name, overwrite = TRUE)


#LeeY ####
setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/GlioVis Datasets")
LeeY <- readRDS("LeeY.Rds")
df <- LeeY[["expr"]] 
colnames(df) <- make.unique(colnames(df))
df <- df %>%
  dplyr::filter(Histology == "GBM") %>%
  dplyr::select(-c(Sample, Histology, Grade, Recurrence, Subtype, CIMP_status, CIMP_status.1, survival, status)) %>%
  t()

df_linear <- 2^df  # genes in rows, samples in columns
df_tpm <- sweep(df_linear, 2, colSums(df_linear), FUN = "/") * 1e6
df_tpm_matrix <- t(df_tpm)  # transpose: samples as rows
df_tpm_qn <- normalize.quantiles(df_tpm_matrix)
rownames(df_tpm_qn) <- rownames(df_tpm_matrix)
colnames(df_tpm_qn) <- colnames(df_tpm_matrix)
df_tpm_qn_final <- t(df_tpm_qn)

#get into correct data type
bulk <- as.data.frame(df_tpm_qn_final)
bulk_ma <- as.matrix(bulk)
bulk_ar <- array(bulk_ma, dim = c(nrow(bulk_ma), ncol(bulk_ma)))
rownames(bulk_ar) <- rownames(bulk_ma)
colnames(bulk_ar) <- colnames(bulk_ma)
class(bulk_ar)

# take GBM_refPhi for reference
setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME")
GBM_refPhi <- readRDS("GBM_refPhi.RDS")
length(intersect(rownames(bulk_ar), rownames(GBM_refPhi@phi.cs)))

InstaPrism.res <- InstaPrism(input_type = 'refPhi_cs', filter=FALSE, bulk_Expr = bulk_ar, refPhi_cs = GBM_refPhi, n.core = 16)

setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Instaprism Results")
wb <- createWorkbook()
addWorksheet(wb, "Cell Types")
writeData(wb, "Cell Types", InstaPrism.res@Post.ini.ct@theta, rowNames = TRUE)
addWorksheet(wb, "Cell States")
writeData(wb, "Cell States", InstaPrism.res@Post.ini.cs@theta, rowNames = TRUE)
addWorksheet(wb, "Cell Types Transposed")
writeData(wb, "Cell Types Transposed", t(InstaPrism.res@Post.ini.ct@theta), rowNames = TRUE)
addWorksheet(wb, "Cell States Transposed")
writeData(wb, "Cell States Transposed", t(InstaPrism.res@Post.ini.cs@theta), rowNames = TRUE)

file_name <- "LeeY Microarray InstaPrism results.xlsx"
saveWorkbook(wb, file = file_name, overwrite = TRUE)


#Phillips ####
setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/GlioVis Datasets")
Phillips <- readRDS("Phillips.Rds")
df <- Phillips[["expr"]] %>%
  dplyr::filter(Histology == "GBM") %>%
  dplyr::select(-c(Sample, Histology, Grade, Recurrence, Subtype, CIMP_status, survival, status)) %>%
  t()

df_linear <- 2^df  # genes in rows, samples in columns
df_tpm <- sweep(df_linear, 2, colSums(df_linear), FUN = "/") * 1e6
df_tpm_matrix <- t(df_tpm)  # transpose: samples as rows
df_tpm_qn <- normalize.quantiles(df_tpm_matrix)
rownames(df_tpm_qn) <- rownames(df_tpm_matrix)
colnames(df_tpm_qn) <- colnames(df_tpm_matrix)
df_tpm_qn_final <- t(df_tpm_qn)

#get into correct data type
bulk <- as.data.frame(df_tpm_qn_final)
bulk_ma <- as.matrix(bulk)
bulk_ar <- array(bulk_ma, dim = c(nrow(bulk_ma), ncol(bulk_ma)))
rownames(bulk_ar) <- rownames(bulk_ma)
colnames(bulk_ar) <- colnames(bulk_ma)
class(bulk_ar)

# take GBM_refPhi for reference
setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME")
GBM_refPhi <- readRDS("GBM_refPhi.RDS")
length(intersect(rownames(bulk_ar), rownames(GBM_refPhi@phi.cs)))

InstaPrism.res <- InstaPrism(input_type = 'refPhi_cs', filter=FALSE, bulk_Expr = bulk_ar, refPhi_cs = GBM_refPhi, n.core = 16)

setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Instaprism Results")
wb <- createWorkbook()
addWorksheet(wb, "Cell Types")
writeData(wb, "Cell Types", InstaPrism.res@Post.ini.ct@theta, rowNames = TRUE)
addWorksheet(wb, "Cell States")
writeData(wb, "Cell States", InstaPrism.res@Post.ini.cs@theta, rowNames = TRUE)
addWorksheet(wb, "Cell Types Transposed")
writeData(wb, "Cell Types Transposed", t(InstaPrism.res@Post.ini.ct@theta), rowNames = TRUE)
addWorksheet(wb, "Cell States Transposed")
writeData(wb, "Cell States Transposed", t(InstaPrism.res@Post.ini.cs@theta), rowNames = TRUE)

file_name <- "Phillips Microarray InstaPrism results.xlsx"
saveWorkbook(wb, file = file_name, overwrite = TRUE)



#Freije ####
setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/GlioVis Datasets")
Freije <- readRDS("Freije.Rds")
df <- Freije[["expr"]] %>%
  dplyr::filter(Histology == "GBM") %>%
  dplyr::select(-c(Sample, Histology, Grade, Recurrence, Subtype, CIMP_status, survival, status)) %>%
  t()

df_linear <- 2^df  # genes in rows, samples in columns
df_tpm <- sweep(df_linear, 2, colSums(df_linear), FUN = "/") * 1e6
df_tpm_matrix <- t(df_tpm)  # transpose: samples as rows
df_tpm_qn <- normalize.quantiles(df_tpm_matrix)
rownames(df_tpm_qn) <- rownames(df_tpm_matrix)
colnames(df_tpm_qn) <- colnames(df_tpm_matrix)
df_tpm_qn_final <- t(df_tpm_qn)

#get into correct data type
bulk <- as.data.frame(df_tpm_qn_final)
bulk_ma <- as.matrix(bulk)
bulk_ar <- array(bulk_ma, dim = c(nrow(bulk_ma), ncol(bulk_ma)))
rownames(bulk_ar) <- rownames(bulk_ma)
colnames(bulk_ar) <- colnames(bulk_ma)
class(bulk_ar)

# take GBM_refPhi for reference
setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME")
GBM_refPhi <- readRDS("GBM_refPhi.RDS")
length(intersect(rownames(bulk_ar), rownames(GBM_refPhi@phi.cs)))

InstaPrism.res <- InstaPrism(input_type = 'refPhi_cs', filter=FALSE, bulk_Expr = bulk_ar, refPhi_cs = GBM_refPhi, n.core = 16)

setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Instaprism Results")
wb <- createWorkbook()
addWorksheet(wb, "Cell Types")
writeData(wb, "Cell Types", InstaPrism.res@Post.ini.ct@theta, rowNames = TRUE)
addWorksheet(wb, "Cell States")
writeData(wb, "Cell States", InstaPrism.res@Post.ini.cs@theta, rowNames = TRUE)
addWorksheet(wb, "Cell Types Transposed")
writeData(wb, "Cell Types Transposed", t(InstaPrism.res@Post.ini.ct@theta), rowNames = TRUE)
addWorksheet(wb, "Cell States Transposed")
writeData(wb, "Cell States Transposed", t(InstaPrism.res@Post.ini.cs@theta), rowNames = TRUE)

file_name <- "Freije Microarray InstaPrism results.xlsx"
saveWorkbook(wb, file = file_name, overwrite = TRUE)



#Murat ####
setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/GlioVis Datasets")
Murat <- readRDS("Murat.Rds")
df <- Murat[["expr"]] %>%
  dplyr::filter(Histology == "GBM") %>%
  dplyr::select(-c(Sample, Histology, Grade, Recurrence, Subtype, CIMP_status, survival, status)) %>%
  t()

df_linear <- 2^df  # genes in rows, samples in columns
df_tpm <- sweep(df_linear, 2, colSums(df_linear), FUN = "/") * 1e6
df_tpm_matrix <- t(df_tpm)  # transpose: samples as rows
df_tpm_qn <- normalize.quantiles(df_tpm_matrix)
rownames(df_tpm_qn) <- rownames(df_tpm_matrix)
colnames(df_tpm_qn) <- colnames(df_tpm_matrix)
df_tpm_qn_final <- t(df_tpm_qn)

#get into correct data type
bulk <- as.data.frame(df_tpm_qn_final)
bulk_ma <- as.matrix(bulk)
bulk_ar <- array(bulk_ma, dim = c(nrow(bulk_ma), ncol(bulk_ma)))
rownames(bulk_ar) <- rownames(bulk_ma)
colnames(bulk_ar) <- colnames(bulk_ma)
class(bulk_ar)

# take GBM_refPhi for reference
setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME")
GBM_refPhi <- readRDS("GBM_refPhi.RDS")
length(intersect(rownames(bulk_ar), rownames(GBM_refPhi@phi.cs)))

InstaPrism.res <- InstaPrism(input_type = 'refPhi_cs', filter=FALSE, bulk_Expr = bulk_ar, refPhi_cs = GBM_refPhi, n.core = 16)

setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Instaprism Results")
wb <- createWorkbook()
addWorksheet(wb, "Cell Types")
writeData(wb, "Cell Types", InstaPrism.res@Post.ini.ct@theta, rowNames = TRUE)
addWorksheet(wb, "Cell States")
writeData(wb, "Cell States", InstaPrism.res@Post.ini.cs@theta, rowNames = TRUE)
addWorksheet(wb, "Cell Types Transposed")
writeData(wb, "Cell Types Transposed", t(InstaPrism.res@Post.ini.ct@theta), rowNames = TRUE)
addWorksheet(wb, "Cell States Transposed")
writeData(wb, "Cell States Transposed", t(InstaPrism.res@Post.ini.cs@theta), rowNames = TRUE)

file_name <- "Murat Microarray InstaPrism results.xlsx"
saveWorkbook(wb, file = file_name, overwrite = TRUE)



#Joo ####
setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/GlioVis Datasets")
Joo <- readRDS("Joo.Rds")
df <- Joo[["expr"]] %>%
  dplyr::filter(Histology == "GBM") %>%
  dplyr::select(-c(Sample, Histology, Grade, Recurrence, Subtype, CIMP_status, survival, status)) %>%
  t()

df_linear <- 2^df  # genes in rows, samples in columns
df_tpm <- sweep(df_linear, 2, colSums(df_linear), FUN = "/") * 1e6
df_tpm_matrix <- t(df_tpm)  # transpose: samples as rows
df_tpm_qn <- normalize.quantiles(df_tpm_matrix)
rownames(df_tpm_qn) <- rownames(df_tpm_matrix)
colnames(df_tpm_qn) <- colnames(df_tpm_matrix)
df_tpm_qn_final <- t(df_tpm_qn)

#get into correct data type
bulk <- as.data.frame(df_tpm_qn_final)
bulk_ma <- as.matrix(bulk)
bulk_ar <- array(bulk_ma, dim = c(nrow(bulk_ma), ncol(bulk_ma)))
rownames(bulk_ar) <- rownames(bulk_ma)
colnames(bulk_ar) <- colnames(bulk_ma)
class(bulk_ar)

# take GBM_refPhi for reference
setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME")
GBM_refPhi <- readRDS("GBM_refPhi.RDS")
length(intersect(rownames(bulk_ar), rownames(GBM_refPhi@phi.cs)))

InstaPrism.res <- InstaPrism(input_type = 'refPhi_cs', filter=FALSE, bulk_Expr = bulk_ar, refPhi_cs = GBM_refPhi, n.core = 16)

setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Instaprism Results")
wb <- createWorkbook()
addWorksheet(wb, "Cell Types")
writeData(wb, "Cell Types", InstaPrism.res@Post.ini.ct@theta, rowNames = TRUE)
addWorksheet(wb, "Cell States")
writeData(wb, "Cell States", InstaPrism.res@Post.ini.cs@theta, rowNames = TRUE)
addWorksheet(wb, "Cell Types Transposed")
writeData(wb, "Cell Types Transposed", t(InstaPrism.res@Post.ini.ct@theta), rowNames = TRUE)
addWorksheet(wb, "Cell States Transposed")
writeData(wb, "Cell States Transposed", t(InstaPrism.res@Post.ini.cs@theta), rowNames = TRUE)

file_name <- "Joo Microarray InstaPrism results.xlsx"
saveWorkbook(wb, file = file_name, overwrite = TRUE)




#Ducray ####
setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/GlioVis Datasets")
Ducray <- readRDS("Ducray.Rds")
df <- Ducray[["expr"]] %>%
  dplyr::filter(Histology == "GBM") %>%
  dplyr::select(-c(Sample, Histology, Grade, Recurrence, Subtype, CIMP_status, survival, status)) %>%
  t()

df_linear <- 2^df  # genes in rows, samples in columns
df_tpm <- sweep(df_linear, 2, colSums(df_linear), FUN = "/") * 1e6
df_tpm_matrix <- t(df_tpm)  # transpose: samples as rows
df_tpm_qn <- normalize.quantiles(df_tpm_matrix)
rownames(df_tpm_qn) <- rownames(df_tpm_matrix)
colnames(df_tpm_qn) <- colnames(df_tpm_matrix)
df_tpm_qn_final <- t(df_tpm_qn)

#get into correct data type
bulk <- as.data.frame(df_tpm_qn_final)
bulk_ma <- as.matrix(bulk)
bulk_ar <- array(bulk_ma, dim = c(nrow(bulk_ma), ncol(bulk_ma)))
rownames(bulk_ar) <- rownames(bulk_ma)
colnames(bulk_ar) <- colnames(bulk_ma)
class(bulk_ar)

# take GBM_refPhi for reference
setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME")
GBM_refPhi <- readRDS("GBM_refPhi.RDS")
length(intersect(rownames(bulk_ar), rownames(GBM_refPhi@phi.cs)))

InstaPrism.res <- InstaPrism(input_type = 'refPhi_cs', filter=FALSE, bulk_Expr = bulk_ar, refPhi_cs = GBM_refPhi, n.core = 16)

setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Instaprism Results")
wb <- createWorkbook()
addWorksheet(wb, "Cell Types")
writeData(wb, "Cell Types", InstaPrism.res@Post.ini.ct@theta, rowNames = TRUE)
addWorksheet(wb, "Cell States")
writeData(wb, "Cell States", InstaPrism.res@Post.ini.cs@theta, rowNames = TRUE)
addWorksheet(wb, "Cell Types Transposed")
writeData(wb, "Cell Types Transposed", t(InstaPrism.res@Post.ini.ct@theta), rowNames = TRUE)
addWorksheet(wb, "Cell States Transposed")
writeData(wb, "Cell States Transposed", t(InstaPrism.res@Post.ini.cs@theta), rowNames = TRUE)

file_name <- "Ducray Microarray InstaPrism results.xlsx"
saveWorkbook(wb, file = file_name, overwrite = TRUE)



