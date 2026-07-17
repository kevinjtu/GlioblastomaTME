## Figure S6
## Purpose Of Script: Conduct analysis of each Xenium sample to determine if
## lymphoid-like niche is present. ImageDimPlots were used to show the
## lymphoid-like niche is colocalized with the perinecrotic niche in each 
## sample. Differences in gene expression were compared in immune cells types
## in the lymphoid-like niche and perinecrotic niche of each sample.
##
## Author: George D. Dalton
##
## Date: July 16, 2026

ls() ## list objects currently present in workspace
rm(list=ls()) ## remove objects from current workspace
ls() ## confirm current workspace is empty

library(remotes) # version 2.4.2
library(Seurat) # version 5.1.0
library(tidyverse) # version 2.0.0
library(sf) # version 1.0.17
library(scCustomize) # version 2.1.2
library(scclusteval) # version 0.0.0.9000
library(ComplexHeatmap) # version 2.14.0
library(BiocParallel) # version 1.32.5
library(EnhancedVolcano) # version 1.16.0
library(clusterProfiler) # version 4.6.2
library(plotly) # version 4.10.4
library(magrittr) # version 2.0.3
library(base64enc) # version 0.1.3
library(shiny) # version 1.9.1

# Set working directory
setwd("/home/gd2")

# Read in a Xenium file. This is a Seurat file stored as an R object
Xenium.obj <- readRDS("/home/gd2/Xenium_Sample.rds")
# Look at Seurat object contents
Xenium.obj

# This is an example: first, subset object to include immune cells only 
# and store in variable jj
# can also further subset jj so it contains only the perinecrotic and 
# lymphoid-like niches
jj <- subset(Xenium.obj, subset = celltype == c("Macrophage_1", "Macrophage_2", "Microglia_1_activated", "Microglia_2", "Microglia_3"))

# Determine current cell identities of jj
Idents(object = jj)

# Set current cell identities of jj to "niches"
Idents(object = jj) <- "niches"

# Check that current cell identities of jj are "niches"
Idents(jj)

# Generate a dimplot and use "molecules" to look at abundance of specific 
# genes in plot, "cols" lets you color the two niches of interest if object
# was subsetted to only contain perinecrotic and lymphoid-like niches
ImageDimPlot(jj, fov = "fov7", molecules = c(""), nmols = 20000, cols = c("blue", "yellow"), size = 2.5, axes = TRUE)

# Look at expression of marker genes in immune cell types and compare
# expression of the genes between the two niches
# Set default assay of jj to SCT assay
DefaultAssay(jj) <- "SCT"

# Check current cell identities of jj
Idents(jj)

# Set jj current cell identities to "niches"
Idents(object = jj) <- "niches"

# Confirm jj current cell identities set to "niches"
Idents(jj)

# Generate variables to look at certain gene signatures via dotplot 
features_to_plot_3 <- c("CD68", "CD163", "CD14", "CXCR4") # macrophage markers
features_to_plot_4 <- c("PDCD1", "CD274", "EOMES", "BATF3") # T cell exhaustion markers
features_to_plot_5 <- c("CD19", "CD1C", "CD1B", "KLRB1") # B cell, dendritic cell, NK markers
features_to_plot_6 <- c("GZMB", "GZMA", "GNLY", "CD2", "CD4", "CXCR6", "CD3G", "THEMIS", "TESPA1", "IL7R") # T cell markers
features_to_plot_7 <- c("IFNA1", "IFNA2", "IFNA5", "IFNA6", "IFNA10", "IFNA21", "IFNL1") # interferon markers
features_to_plot_8 <- c("CCL2", "CXCL10", "CXCL12") # inflammatory markers
features_to_plot_9 <- c("VEGFA", "HILPDA") # hypoxia/perinecrotic niche markers
features_to_plot_10 <- c("CD36", "IL6", "CSPG4", "CXCL11", "PCNA", "MKI67", "APOE", "IGFBP4", "STAT3") # senescence markers

# Generate dotplots and compare gene expression between niches
DotPlot(jj, features = features_to_plot_3) + RotatedAxis()
DotPlot(jj, features = features_to_plot_4) + RotatedAxis()
DotPlot(jj, features = features_to_plot_5) + RotatedAxis()
DotPlot(jj, features = features_to_plot_6) + RotatedAxis()
DotPlot(jj, features = features_to_plot_7) + RotatedAxis()
DotPlot(jj, features = features_to_plot_8) + RotatedAxis()
DotPlot(jj, features = features_to_plot_9) + RotatedAxis()
DotPlot(jj, features = features_to_plot_10) + RotatedAxis()

# Run dotplot and store dotplot info for each gene signature in variable dp1
dp1 <- DotPlot(jj, features = features_to_plot_3) + RotatedAxis()

# Extract raw data from variable dp1 and store as variable plot_data
plot_data <- dp1$data

# Look at first six rows of data in variable plot_data
head(plot_data)

# Convert plot_data to data frame and store in variable how1
how1 <- as.data.frame(plot_data) # convert to data frame

# save as CSV, this data will be used to make the heatmap
write.table(how1, sep=",", file = "Sample_T_Cell_Markers.csv", 
            row.names=TRUE, col.names = NA, quote=FALSE)

## To generate the heatmap ##
# Load CSV file and store in variable just.raw.counts, 
# CSV file has data for all genes of interest and all samples
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
pheatmap(data_subset, cluster_cols=FALSE, cluster_rows=FALSE)

# Put Formula To Calculate Z-Score In A Vector - Z-Score Lets You Scale 
# The Rows
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

# Calculate Row Z-Scores For Our Data Matrix
data_subset_norm <- t(apply(data_subset, 1, cal_z_score))

# Convert data_subset_norm to a data frame and store as variable called "t"
t <- as.data.frame(data_subset_norm)

# Generate a CSV file with Z scored gene expression values
write.table(t, sep=",", file = "Heatmap_Scaled_All_8_Samples.csv", 
            row.names=TRUE, col.names = NA, quote=FALSE)

# Run PHEATMAP On Z-Scores (not final published heatmap)
pheatmap(data_subset_norm, cluster_cols=FALSE, cluster_rows=FALSE)

## Prepare Final Heatmap ##

# Add Colors To The Row And Column Annotations and store in variable my_colour
# Group is for study treatment groups, Pvalue denotes perinecrotic or lymphoid-
# like niche, Gene denotes marker genes that were examined
my_colour = list(Pvalue = c(A = "orange", B = "blue"),
                 Group = c(No_Tx="red", Tx ="yellow", D2C7="chartreuse"),
                 Gene = c(cytotoxic = "cadetblue1", marker = "khaki2", exhaustion = "palegreen1", interferon = "snow2", macros = "lightsalmon", Bcell = "mediumorchid1", hypoxia = "plum1", inflammatory = "black", senescence = "green"))

# Read in table with niche ids and treatment information and store as 
# variable res
res <-read.table(file.choose(), sep = ",", header = TRUE)
class(res) # should be data.frame
head(res) # look at first 6 rows of table
dim(res) # look at table dimensions
res$Group <- factor(res$Group, levels = c("No_Tx", "Tx", "D2C7"))
res$Pvalue <- factor(res$Pvalue, levels = c("A", "B"))
# Store res in the variable annotation_col
annotation_col = res


# Read in table with gene information
genes <-read.table(file.choose(), sep = ",", header = TRUE)
class(genes) # should be data.frame
head(genes) # look at first 6 rows of table
dim(genes) # look at table dimensions
genes$Gene <- factor(genes$Gene, levels = c("cytotoxic", "marker", "exhaustion", "interferon", "macros", "Bcell", "hypoxia", "inflammatory", "senescence"))
# Store res in the variable annotation_col
annotation_row = genes

# Set breaks for Z score legend
my_breaks <- c(-2, -0.5, -0.05, 0, .12, 1, 2)
# Set colors for breaks
my_colors <- c("black", "blue", "cyan", "green", "yellow", "red", "white")

# Generate heatmap
pheatmap(data_subset_norm,
         color = my_colors,
         breaks = my_breaks,
         annotation_colors=my_colour,
         annotation_row = annotation_row,
         annotation_col = annotation_col,
         cluster_cols=FALSE, 
         cluster_rows=FALSE, 
         gaps_col = c(8),
         gaps_row = c(4,10,14,21,24,28,30,33,34),
         cellheight = 12, 
         cellwidth = 30,
         fontsize = 13, legend = TRUE)


