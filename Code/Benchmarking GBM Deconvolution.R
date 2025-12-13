#libraries####
library(data.table)
library(R.utils)
library(ggpointdensity)
library(msigdbr)
library(dplyr)
library(textshape)

library(devtools)
setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/deconvBenchmarking-main")
devtools::install(".")
3
#restart R session
library(deconvBenchmarking)

#Data Input####
options(timeout = 3000)

setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Single Cell Atlas GBM")

pdata_raw <- as.data.frame(fread("IDHwt.GBM.Metadata.SS2.txt")) %>%
  dplyr::slice(-1) %>%
  dplyr::select(NAME, Sample, CellAssignment) %>%
  dplyr::rename(cells = NAME, samples = Sample, cell.types = CellAssignment)

tpm_raw <- as.data.frame(fread("IDHwtGBM.processed.SS2.logTPM.txt"))

tpm_processed = tpm_raw %>% column_to_rownames("GENE")
tpm_unlog = (2^tpm_processed - 1) * 10

# --- Prepare Metadata (scMeta_prepared) ---
scMeta_prepared <- data.table::copy(pdata_raw) # Start with a fresh copy
scMeta_prepared <- as.data.table(scMeta_prepared) # Convert scMeta_prepared to a data.table

# 1. Rename cell type column (This part was successful)
original_cell_type_col <- "cell.types" 
target_cell_type_col_literal <- "colnames_of_cellType" # Literal name expected by function

if (original_cell_type_col %in% colnames(scMeta_prepared)) {
  setnames(scMeta_prepared, old = original_cell_type_col, new = target_cell_type_col_literal)
  print(paste("Renamed column '", original_cell_type_col, "' to '", target_cell_type_col_literal, "' in scMeta_prepared.", sep=""))
} else {
  warning(paste("Expected cell type column '", original_cell_type_col, "' not found. Available columns: ", paste(colnames(scMeta_prepared), collapse=", "), sep=""))
}

# 2. Handle sample column -- ADJUSTED FOR THE NEW ERROR
# The error indicates the function is looking for a column literally named 'colnames_of_sample'.
# Your file's actual sample column (after loading) is 'samples'.
actual_sample_col_in_file <- "samples"
target_sample_col_literal <- "colnames_of_sample" # Literal name the function seems to be looking for

if (actual_sample_col_in_file %in% colnames(scMeta_prepared)) {
  # Rename the 'samples' column to 'colnames_of_sample'
  setnames(scMeta_prepared, old = actual_sample_col_in_file, new = target_sample_col_literal)
  print(paste("Renamed sample column '", actual_sample_col_in_file, "' to '", target_sample_col_literal, "' in scMeta_prepared to match literal expectation from error.", sep=""))
} else {
  warning(paste("Original sample column '", actual_sample_col_in_file, "' not found in scMeta_prepared. Available columns: ", paste(colnames(scMeta_prepared), collapse=", "), sep=""))
}

# 3. Identify and key by the cell identifier column (This part was successful)
actual_cell_id_col_in_file <- "cells" 

if (actual_cell_id_col_in_file %in% colnames(scMeta_prepared)) {
  # Basic check for cell ID consistency (optional but recommended)
  if (!all(colnames(tpm_unlog) %in% scMeta_prepared[[actual_cell_id_col_in_file]])) {
    warning(paste("Caution: Not all cell IDs from expression data (tpm_unlog colnames) are present in the '", actual_cell_id_col_in_file, "' column of scMeta_prepared.", sep=""))
  }
  if (!all(scMeta_prepared[[actual_cell_id_col_in_file]] %in% colnames(tpm_unlog))) {
    warning(paste("Caution: Not all cell IDs from scMeta_prepared$' ", actual_cell_id_col_in_file, "' are present as colnames in tpm_unlog. Consider filtering scMeta_prepared.", sep=""))
  }
  
  setkeyv(scMeta_prepared, actual_cell_id_col_in_file)
  print(paste("Keyed scMeta_prepared by column: '", actual_cell_id_col_in_file, "'.", sep=""))
} else {
  stop(paste("Expected cell identifier column '", actual_cell_id_col_in_file, "' not found in scMeta_prepared. This is crucial. Available columns: ", paste(colnames(scMeta_prepared), collapse=", "), sep=""))
}

#Create benchmarking_init object (ensure arguments match the prepared column names) ####
benchmarking_obj = benchmarking_init(
  scExpr = tpm_unlog,
  scMeta = scMeta_prepared, 
  colnames_of_cellType = target_cell_type_col_literal, # This is "colnames_of_cellType"
  colnames_of_sample = target_sample_col_literal,     # This is "colnames_of_sample"
  
  nbulk = 100,
  fixed_cell_type = 'Malignant', # <--- CORRECTED TO MATCH YOUR DATA
  training_ratio = 0.5,
  split_by = "cell", 
  bulkSimulator_methods = c('homo','semi','heter'),
  heter_cell_type = 'Malignant', # <--- CORRECTED TO MATCH YOUR DATA
  dirichlet_cs_par = 0.05,
  min.subcluster.size = 10,
  refMarkers_methods = c('limma','scran'),
  refMatrix_methods = c('raw','limma','scran'),
  include_tcga = FALSE, 
  #tcga_abbreviation = 'SKCM',
  purity_methods =  c('ESTIMATE', 'ABSOLUTE', 'LUMP', 'IHC', 'CPE')
)

print("benchmarking_init call completed if no errors.")

#Deconvolute pseudobuk data####
#DEBUG: convert into data.frame 
# actual_cell_id_col_in_file should still be "cells"
if (!exists("actual_cell_id_col_in_file") || is.null(actual_cell_id_col_in_file)) {
  actual_cell_id_col_in_file <- "cells" 
  warning(paste("actual_cell_id_col_in_file was not defined, defaulting to:", actual_cell_id_col_in_file))
}

# Get the order of cell IDs from the expression matrix
cell_order_from_expr <- colnames(tpm_unlog)

# Convert scMeta_prepared to data.frame
scMeta_df_temp <- as.data.frame(scMeta_prepared)

# Ensure the cell ID column exists in scMeta_df_temp
if (!(actual_cell_id_col_in_file %in% names(scMeta_df_temp))) {
  stop(paste0("Cell ID column '", actual_cell_id_col_in_file, "' not found in scMeta_prepared. Cannot proceed."))
}

# Match and reorder scMeta_df_temp based on the expression matrix's cell order
# This creates an index for reordering
match_indices <- match(cell_order_from_expr, scMeta_df_temp[[actual_cell_id_col_in_file]])

# Check if all cell IDs from expression matrix were found in metadata
if (any(is.na(match_indices))) {
  missing_cells_in_meta <- cell_order_from_expr[is.na(match_indices)]
  warning(paste0("Some cell IDs from scExpr (tpm_unlog) were not found in scMeta_prepared's '", 
                 actual_cell_id_col_in_file, "' column. Number missing: ", length(missing_cells_in_meta)))
  # If critical, you might stop here or filter tpm_unlog to only include cells present in metadata
  # For now, let's proceed with cells that do match.
  # Filter out NAs from match_indices and adjust cell_order_from_expr
  valid_match_indices <- !is.na(match_indices)
  cell_order_from_expr <- cell_order_from_expr[valid_match_indices]
  match_indices <- match_indices[valid_match_indices]
  
  # Also filter tpm_unlog to ensure consistency if some cells from expr are not in meta
  if (length(cell_order_from_expr) < ncol(tpm_unlog)) {
    cat("INFO: Subsetting tpm_unlog to match cells found in metadata.\n")
    tpm_unlog <- tpm_unlog[, cell_order_from_expr, drop = FALSE] # Ensure it's still a matrix
  }
}

if (length(match_indices) == 0) {
  stop("No matching cell IDs found between scExpr and scMeta. Cannot proceed.")
}

# Reorder the temporary data.frame
scMeta_for_deconv <- scMeta_df_temp[match_indices, , drop = FALSE]

# Now set the rownames using the (now correctly ordered) cell ID column
# The rownames will be in the same order as colnames(tpm_unlog)
rownames(scMeta_for_deconv) <- scMeta_for_deconv[[actual_cell_id_col_in_file]] 
# Alternatively, and perhaps safer if there were any NAs in cell_order_from_expr earlier:
# rownames(scMeta_for_deconv) <- cell_order_from_expr

cat(paste0("DEBUG: Converted and reordered scMeta_prepared to 'scMeta_for_deconv'. Rownames set from '", 
           actual_cell_id_col_in_file, "' column, matching colnames(tpm_unlog).\n"))
cat(paste0("DEBUG: ncol(tpm_unlog): ", ncol(tpm_unlog), ", nrow(scMeta_for_deconv): ", nrow(scMeta_for_deconv), "\n"))

# For verification (optional, can remove after confirming)
print("First 5 colnames(tpm_unlog):")
print(head(colnames(tpm_unlog), 5))
print("First 5 rownames(scMeta_for_deconv):")
print(head(rownames(scMeta_for_deconv), 5))
print(paste("Are first 5 identical after reorder/rownames set?", 
             identical(head(colnames(tpm_unlog),5), head(rownames(scMeta_for_deconv),5)) ))
print(paste("Are all identical after reorder/rownames set?", 
             identical(colnames(tpm_unlog), rownames(scMeta_for_deconv)) ))


#deconvolution with benchmarking object
deconvResults = benchmarking_deconv(
  benchmarking_obj,
  marker_based_methods = c('firstPC','debCAM','TOAST'),
  
  # arguments for regression-based methods
  regression_based_methods = c('nnls', 'MuSiC'),
  #cibersort_path = './', 
  scExpr = tpm_unlog,  
  scMeta = scMeta_for_deconv, # <--- USE THE NEWLY PREPARED DATA.FRAME
  colnames_of_cellType = target_cell_type_col_literal, # This is "colnames_of_cellType"
  colnames_of_sample = target_sample_col_literal,   # This is "colnames_of_sample"
  
  refFree_methods = c('debCAM'),
  
  # arguments for Bayesian-based methods
  Bayesian_methods = 'InstaPrism', 
  key = 'Malignant', 
  
  immunedeconv_methods = NULL,
  n.core = 1 # Example, adjust as needed
)

#Evalutate Deconvolution Performance####
performance = benchmarking_evalu(deconvResults)

summarized_performance = gather_performance(performance)

df = summarized_performance[summarized_performance$class %in% c('homo','semi', 'heter'),]
df = df[df$method %in% c('MarkerBased_firstPC_scran', 'MarkerBased_debCAM_scran',
                         'MarkerBased_TOAST_scran',
                         'RefBased_nnls_raw', 'RefBased_MuSiC_raw', 
                         'Bayesian_InstaPrism_updatedAll'),]
df$reg_method = str_extract(df$method, "(?<=_)[^_]+")

pd <- position_dodge(width = 0.1)
ggplot(df,aes(x = reorder(class,cor,median,decreasing = T), 
              y = cor, group = reg_method)) +
  geom_line(position = pd,aes(color=reg_method)) +
  geom_point(position = pd,aes(color=reg_method,shape = reg_method )) +
  scale_color_manual(values = c('#F46D43','#FDAE61',"#8DA0CB","#E78AC3","#A6D854", "black"))+
  guides(linetype = guide_legend("Transmission"))+
  facet_wrap(~cell_type,scales = 'free',nrow = 2)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size = 7),
        strip.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 7),
        axis.title.x = element_text(size = 7),
        axis.title.y = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7),
        strip.background = element_rect(fill = 'white'))+
  xlab('')+
  ylab('pearson r')