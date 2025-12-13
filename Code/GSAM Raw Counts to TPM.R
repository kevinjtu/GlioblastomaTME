# Load necessary libraries####
library(data.table)
library(dplyr)
library(biomaRt)
library(stringr)

#import gsam raw counts data####
setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Expression Data")
gsam_expression <- fread("gsam raw counts.csv")

# Extract Ensembl gene ID
gsam_expression$ensembl_gene_id <- str_extract(gsam_expression$V1, "^ENSG[0-9]+")

# Remove the metadata column
expression_counts <- gsam_expression[, -1, with = FALSE]

# Aggregate counts by gene ID (summing duplicates)
expression_counts_summarized <- expression_counts %>%
  group_by(ensembl_gene_id) %>%
  summarise(across(where(is.numeric), sum), .groups = "drop")

# Extract gene IDs and counts separately
ensembl_ids <- expression_counts_summarized$ensembl_gene_id
raw_counts <- as.data.frame(expression_counts_summarized[, -1])
rownames(raw_counts) <- ensembl_ids

# Use Ensembl with human genes dataset explicitly
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", mirror = "useast")  # "useast" mirror often works

# Try fetching lengths again
gene_lengths <- getBM(
  attributes = c("ensembl_gene_id", "start_position", "end_position"),
  filters = "ensembl_gene_id",
  values = ensembl_ids,
  mart = ensembl
) %>%
  mutate(gene_length_kb = (end_position - start_position + 1) / 1000) %>%
  distinct(ensembl_gene_id, .keep_all = TRUE)

# Merge gene lengths with count matrix
counts_with_length <- expression_counts %>%
  left_join(gene_lengths, by = "ensembl_gene_id") %>%
  filter(!is.na(gene_length_kb))  # remove genes with no length info

# Separate lengths and counts
gene_length_kb <- counts_with_length$gene_length_kb
raw_counts <- counts_with_length %>%
  as.data.frame() %>%
  dplyr::select(-ensembl_gene_id, -start_position, -end_position, -gene_length_kb) %>%
  as.matrix()

# Calculate TPM
# 1. Divide each gene count by its length in kb
rpk <- sweep(raw_counts, 1, gene_length_kb, FUN = "/")

# 2. Scale so that per-sample RPKs sum to 1e6
tpm <- sweep(rpk, 2, colSums(rpk), FUN = "/") * 1e6

# Add back gene names
tpm_df <- as.data.frame(tpm)
tpm_df$ensembl_gene_id <- counts_with_length$ensembl_gene_id

# Optional: reorder columns
tpm_df <- tpm_df %>% dplyr::select(ensembl_gene_id, everything())

#convert ensebly id to hugo id
# Extract Ensembl gene IDs from the TPM data
ensembl_ids <- tpm_df$ensembl_gene_id

# Connect to Ensembl BioMart
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", mirror = "useast")

# Get the mapping: Ensembl ID → HGNC symbol
id_map <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = ensembl_ids,
  mart = ensembl
)

# Ensure only one symbol per Ensembl ID (remove duplicates, keep first)
id_map <- id_map %>% 
  filter(hgnc_symbol != "") %>%
  distinct(ensembl_gene_id, .keep_all = TRUE)

# Merge TPM data with symbols
tpm_annotated <- tpm_df %>%
  left_join(id_map, by = "ensembl_gene_id") %>%
  relocate(hgnc_symbol, .after = ensembl_gene_id)

# Optional: Drop Ensembl column and rename HGNC to gene_id
tpm_final <- tpm_annotated %>%
  dplyr::select(-ensembl_gene_id) %>%
  dplyr::rename(gene_id = hgnc_symbol)

# Preview result
head(tpm_final)

# Save output
write.csv(tpm_final, "gsam_TPM.csv", row.names = FALSE)
