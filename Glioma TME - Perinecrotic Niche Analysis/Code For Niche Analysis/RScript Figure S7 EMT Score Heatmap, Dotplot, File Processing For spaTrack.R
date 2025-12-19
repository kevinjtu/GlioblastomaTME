## Figure S7
## Purpose Of Script: Use the R package UCell to calculate average EMT scores for each tumor cell type in the perinecrotic and non-necrotic niche of each sample. The average EMT scores for each tumor cell type in 
## both niches in each sample were visualized using heatmaps. Dotplots were used to visualize the expression of EMT markers in both niches in all tumor cell types. The python package spaTrack was used to construct detailed 
## tumor cell trajectories for individual Xenium samples and to determine if EMT gene expression changed as specific tumor cells transitioned over time. R version 4.2.2. RStudio version 2022.12.03.353. Python version 3.9.
##
## Author: George D. Dalton
##
## Date: December 15, 2025

## Work in R using UCell to determine EMT scores to make Figure S7A ##

ls() ## list objects currently present in workspace
rm(list=ls()) ## remove objects from current workspace
ls() ## confirm current workspace is empty

## libraries loaded during session
library(remotes) # version 2.4.2
library(Seurat) # version 5.1.0
library(tidyverse) # version 2.0.0
library(sf) # version 1.0.17
library(scCustomize) # version 2.1.2
library(scclusteval) # version 0.0.0.9000
library(ComplexHeatmap) # version 2.14.0
library(BiocParallel) # version 1.32.5
library(ggplot2) # version 3.5.1
library(patchwork) # version 1.3.0
library(SeuratData) # version 0.2.2.9001
library(RColorBrewer) # version 1.1.3
library(UCell) # version 2.2.0

# Set working directory 
setwd("/home/gd2")

# Step 1 Use UCell sample dataset (sample.matrix) to test UCell installation
# Load UCell sample.matrix dataset, this is a sparse matrix of 600 cells and 20729 genes
data(sample.matrix)

# Look at number of rows and columns in dataset
dim(sample.matrix) # should be 20729 600

# Look at first 6 rows of data in the dataset
head(sample.matrix)

# Can now define simple gene sets to get gene signature scores for cells in the sample
gene.sets <- list(Tcell_signature = c("CD2","CD3E","CD3D"),
                  Myeloid_signature = c("SPI1","FCER1G","CSF1R"))

# Run ScoreSignatures_UCell to get signatures scores for all cells
scores <- ScoreSignatures_UCell(sample.matrix, features=gene.sets)

# Look at gene signature scores for first 6 rows of data, this confirms package is installed
head(scores)

# Step 2 Use UCell to determine EMT scores for each tumor cell type in the perinecrotic and non-necrotic niche of each Xenium sample

# Read in a Xenium file (example name: Sample). This is a Seurat file stored as an R object 
Xenium.obj <- readRDS("~/Sample.rds")
# Look at Seurat object contents
Xenium.obj 

# Subset Xenium.obj into sub-objects composed of individual tumor cell types
# Check current cell identities of Xenium.obj
Idents(Xenium.obj) # Levels: SeuratProject

# Set current cell identities of Xenium.obj to “celltype”
Idents(object = Xenium.obj) <- "celltype"

# Confirm current cell identities of Xenium.obj is “celltype”
Idents(Xenium.obj) # should be # of Levels: list of individual cell types

# Subset Xenium.obj into sub-objects composed of individual tumor cell types
# In this example, we do this for Tumor_1 cells
# Subset Xenium.obj into a sub-object with only Tumor_1 cells and store in variable Tumor1
Tumor1 <- subset(Xenium.obj, subset = celltype == "Tumor_1")

# Want EMT signature score for Tumor_1 cells in the perinecrotic niche and non-necrotic niche
# Suppose the perinecrotic niche is niche 4
# To get EMT score for non-necrotic niches, subset Tumor1 by removing perinecrotic niche 4

# Check current identities for Tumor1
Idents(Tumor1) # should be  Levels: Tumor_1

# Remove perinecrotic niche #4 from Tumor1 by subsetting, niche info is in metadata “niches” column
Nonnecrotic_Only <- subset (Tumor1,subset =niches!="4")

# Check cell identities of Nonnecrotic_Only
Idents(Nonnecrotic_Only) # should be  Levels: Tumor_1

# Change cell identities to “niches”
Idents(object = Nonnecrotic_Only) <- "niches"

# Check cell identities
Idents(Nonnecrotic_Only) # should be niches Levels: 1 2 3 5

# Change assay name in Seurat object to RNA
Tumor1cell <- RenameAssays(Nonnecrotic_Only, assay.name = "SCT", new.assay.name = "RNA")

# Get normalized count data from Seurat object
data.matrix <- Tumor1cell@assays$RNA@data

# Look at number of rows and columns in matrix 
dim(data.matrix)
class(data.matrix) # should be a dgCMatrix or Matrix

# Store EMT signature in the variable basic.sign
basic.sign <- list(EMT_signature = c("CAPG", "CDH6", "COL12A1", "CXCL8", "CXCL12", "DCN", 
                                      "FBLN1", "GJA1", "IGFBP3", "IGFBP4", "IL6", "LAMA2", "LOX", 
                                      "PLAUR", "POSTN", "SLIT3", "SPP1", "TGFB1", "TGFBI", "THBS1", 
                                      "TNC", "VCAN", "VEGFA", "VIM"))

# Calculate EMT scores using UCell
scores <- ScoreSignatures_UCell(data.matrix, features=basic.sign)

# Look at first 5 rows of scores
scores[1:5,]

# Convert scores to a data frame
avg.mon <- as.data.frame(scores)

# Save data frame as a CSV file
write.table(avg.mon, sep=",", file = "UCell_Tumor_1.csv", row.names=T)

# The average EMT score is taken for all Tumor 1 cells and shown in a heatmap.
# Average EMT scores are shown for each cell type in perinecrotic & non-necrotic niches in each sample

# Figure 7A Code For Heatmaps
ls() ## list objects currently present in workspace
rm(list=ls()) ## remove objects from current workspace
ls() ## confirm current workspace is empty

## libraries loaded during session
library(ComplexHeatmap) # version 2.14.0

# Load CSV file and store in variable just.raw.counts, load a CSV file for each Xenium sample
# CSV file row 1 is cell type names, col A row 2 is EMT Score (Perinecrotic), col A row 3 is EMT Score (Non-Necrotic), file has average EMT scores for each cell type in the two niches
just.raw.counts<-read.table(file.choose(), sep = ",", header = TRUE, row.names = 1)

# Look at data in just.raw.counts
just.raw.counts

# Store just.raw.counts as variable data_subset
data_subset <- just.raw.counts

# Convert data_subset to a matrix
data_subset <-as.matrix(data_subset)
class(data_subset) # should be “matrix” “array”

# Check dimensions of data_subset
dim(data_subset)

# Look at first 6 rows of data_subset
head(data_subset)

# Make basic heatmap of data_subset (not for publication)
pheatmap(data_subset, cluster_cols=FALSE, cluster_rows=FALSE) # basic heatmap

# Add Colors To The Row And Column Annotations and store in variable my_colour
# Group is for study treatment groups, Pvalue A is p < 0.001, B is p < 0.01, C is p < 0.05, D is not significant
my_colour = list(Group = c(No_Tx="#7570B3", Tx ="deeppink", D2C7="chartreuse"),
                 Pvalue = c(A = "orange", B = "blue", C = "red", D = "black"))

# Read in table with sample group assignments and Pvalue information and store as variable res
res <-read.table(file.choose(), sep = ",", header = TRUE)
class(res) # should be data.frame
head(res) # look at first 6 rows of table
dim(res) # look at table dimensions
res$Group <- factor(res$Group, levels = c("No_Tx", "Tx", "D2C7"))
res$Pvalue <- factor(res$Pvalue, levels = c("A", "B", "C", "D"))
# Store res in the variable annotation_col
annotation_col = res

# Make heatmap for Figure 7A, done for each Xenium sample
pheatmap(data_subset, annotation_colors=my_colour, 
         annotation_col=annotation_col,
         cluster_cols=FALSE, 
         cluster_rows=FALSE,
         gaps_col = c(0),
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
         cellheight=20, cellwidth=30, 
         fontsize = 10, legend = FALSE)

# Make dotplot for Figure 7A, goes below heatmap, make for each Xenium sample
ls() ## list objects currently present in workspace
rm(list=ls()) ## remove objects from current workspace
ls() ## confirm current workspace is empty

## libraries loaded during session
library(SeuratDisk) # version 0.0.0.9021
library(SeuratData) # version 0.2.2.9001
library(Seurat) # version 5.1.0
library(dplyr) # version 1.1.4
library(patchwork) # version 1.3.0
library(ggplot2) # version 3.5.1
library(RColorBrewer) # version 1.1.3

# Set working directory 
setwd(“/home/gd2”)

# Load Xenium/Seurat file for a sample that has been subsetted so it only contains tumor cells
Xenium.obj <- readRDS("~/Sample.rds")
# Look at Seurat object
Xenium.obj 

# Only want cells in perinecrotic niche (e.g. niche 5), subset Xenium.obj so only contains niche 5
# In Xenium.obj, niche info is in “niches” column in metadata
AllCells <- subset (Xenium.obj, subset = niches == "5")
# Look at AllCells
AllCells
# Set DefaultAssay to “SCT”
DefaultAssay(AllCells) <- "SCT"

# define variable "markers.to.plot" with EMT marker genes to look at in AllCells perinecrotic niche
markers.to.plot <- c("VIM", "VEGFA", "TNC", "SPP1", "LOX", "APP")

# set levels to order how cell types appear on dot plot
levels(AllCells) <- c('Tumor_OPC-like_2', 'Tumor_OPC-like_1', 'Tumor_1',
                      'Tumor_OPC-like_3', 'Tumor_AC-like_1', 'Tumor_5',
                      'Tumor_2', 'Tumor_AC-like_2', 'Tumor_3',
                      'Tumor_4')

# Now, make the dotplot
DotPlot(AllCells, features = markers.to.plot) +
  RotatedAxis() + coord_flip() + scale_colour_gradient2(low = "purple", mid = "turquoise", high = "red")

## Figure S7B-G ##
## Initial processing of Seurat object in R to make gene expression & spatial coordinate TSV files
## TSV files used to make anndata object in Python, work done in Jupyter Notebook
## Anndata object used for spaTrack analysis, work done in Jupyter Notebook

## Work done in R to make gene expression, spatial coordinate TSV files

ls() ## list objects currently present in workspace
rm(list=ls()) ## remove objects from current workspace
ls() ## confirm current workspace is empty

## libraries loaded during session
library(Seurat) # version 5.1.0
library(SeuratData) # version 0.2.2.9001
library(pals) # version 1.9
library(gridExtra) # version 2.3
library(SpatialExperiment) # version 1.8.1
library(SummarizedExperiment) # version 1.28.0
library(scuttle) # version 1.8.4
library(scater) # version 1.26.1
library(cowplot) # version 1.1.3

# Set working directory 
setwd(“/home/gd2”)

# Read in a Xenium file (example name: Sample). This is a Seurat file stored as an R object 
Xenium.obj <- readRDS("~/Sample.rds")
# Look at Seurat object contents
Xenium.obj 

# Subset Xenium.obj into a sub-object composed only of tumor cell types

# Remove non-necrotic niche info from Xenium.obj with only tumor cells
# Example, if perinecrotic niche is niche 5
# Set current cell identities of Xenium.obj to “niches”
Idents(object = Xenium.obj) <- "niches"

# Check cell identities
Idents(Xenium.obj) # should be # of Levels: 1 2 3 4 5

# Subset Xenium.obj to remove non-necrotic niche info
Necrotic_Only <- subset (Xenium.obj, subset = niches == "5")

# Confirm Necrotic_Only contains only perinecrotic niche 5
Idents(Necrotic_Only) # should be # of Levels: 5

# Save metadata table of Necrotic_Only in variable h
h <- (Necrotic_Only@meta.data)
head(h) # look at first 6 rows 
avg.mon <- as.data.frame(h) # convert to data frame then save as CSV file
write.table(avg.mon, sep=",", file = " NecroticOnly_metadata.csv", row.names=T)

# Extract gene expression data Necrotic_Only and store as variable gcm
gcm <- Necrotic_Only@assays$SCT@data
dim(gcm) # look at table dimensions
head(gcm) # look at first 6 rows of table
class(gcm) # look at type of object, should be a matrix
# Convert gcm to a data frame stored in variable response1
response1 <- as.data.frame(gcm)
# Transpose data so genes are columns and cell ids are rows
transposed_df <- t(response1) 
# Confirm data is transposed
head(transposed_df) 
# Save this as a TSV file for use in Python analysis
write.table(transposed_df, sep="\t", file="Sample_TransposedGeneExpression.tsv", row.names=TRUE, col.names=NA, quote=FALSE)

# Get spatial centroid info from Necrotic_Only in file and store as csv file, if fov = fov4
y <- (Necrotic_Only@images$fov4)
z <- GetTissueCoordinates(y)
head(z) # look at first 6 rows of table
# Save as a CSV file
write.table(z, sep=",", file = "Necrotic Only_TissueCoords_All.csv", row.names=T)

# Use metadata CSV file and Tissue Coords CSV file for Necrotic_Only to make a CSV file with 
# the following info: col A has cell types and is titled “cluster”, col B has X coords and is titled “x”, col C 
# has Y coords and is titled “y”, col D has cell IDs and is titled “ID”
# Read in this CSV file and store in variable res
res <-read.table(file.choose(), sep = ",", header = TRUE)
class(res) # check class
head(res) # look at first 6 rows
dim(res) # check dimensions for accuracy
# If file looks good, save as a TSV file for use in Python analysis
write.table(res, sep="\t", file="Sample_TissueCoords.tsv", col.names=NA, quote=FALSE)

