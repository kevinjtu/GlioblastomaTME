## Figures 5, S3, S4, S5A,B,D
## Purpose of Script: Data from Xenium samples was stored in a Seurat object and saved as an R object. Xenium samples were loaded into R and analyzed. UMAP plots were made of cell type clusters present in each sample. 
## Niche analysis was conducted using Seurat’s BuildNicheAssay and each sample was organized into 5 niches. A perinecrotic niche was identified in each sample using this analysis. The perinecrotic niche was overlaid on 
## the spatial plot for each sample to show its location. VEGFA and LOX expression were overlaid on the sample spatial plot to show these genes were predominantly expressed within the perinecrotic niche. TMEM173 (STING) 
## expression was also overlaid on the sample spatial plot to show STING expression is found mainly in the non-necrotic niche. The cell type composition of each niche was determined and graphed for each sample. R version 4.2.2. 
## RStudio version 2022.12.03.353.
##
## Author: George D. Dalton
##
## Date: December 15, 2025

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

# Set working directory 
setwd("/home/gd2")

# Read in a Xenium file (example name: Sample). This is a Seurat file stored as an R object 
Xenium.obj <- readRDS("~/Sample.rds")
# Look at Seurat object contents
Xenium.obj

# Generate Figures 5A-C, S3A-B, S4A-B, S5A
DimPlot(Xenium.obj, reduction = "umap", group.by = "celltype", label = F) + DarkTheme()

# Perform Build Niche Analysis on Sample, for this sample fov = fov2, generate 5 niches                                                  
xenium.obj <- BuildNicheAssay(object = Xenium.obj, fov = "fov2", group.by = "celltype",
                              niches.k = 5, neighbors.k = 30)

# Generate a spatial plot showing the 5 niches identified by BuildNicheAssay, store in niche.plot
niche.plot <- ImageDimPlot(xenium.obj, group.by = "niches", size = 2, dark.background = F) + ggtitle("Niches") +
  scale_fill_manual(values = c("#442288", "#6CA2EA", "#B5D33D", "#FED23F", "#EB7D5B"))
# Look at niche.plot figure to see 5 niches overlaid on spatial plot
niche.plot

# Look at DEGs in each niche
# Do this using Seurat’s FindMarkers function
# Set DefaultAssay to SCT
DefaultAssay(xenium.obj) <- "SCT"
# Look at cell identities of xenium.obj
Idents(xenium.obj)
# Change cell identities to “niches”
Idents(object = xenium.obj) <- "niches"
# Confirm new cell identities
Idents(xenium.obj) # should be Levels: 1 2 3 4 5
# Here, we use FindMarkers to look at DEGs in niche 2 versus other niches
# DEGs at Padjust < 0.05 and Log2FC < -0.5 or > 0.5
Niche2.markers <- FindMarkers(xenium.obj, only.pos = FALSE, ident.1 = "2", ident.2 = NULL)
# Look at first 6 rows of DEGs
head(Niche2.markers)
# Convert to data frame
avg.mon <- as.data.frame(Niche2.markers)
# save DEG info as a CSV file, DEGs will have padjust < 0.05 and Log2FC < -0.5 or > 0.5
write.table(Niche2.markers, sep=",", file = "Niche2Markers_All.csv", row.names=T)

# Note: Perinecrotic Niche will have high expression of genes like VEGFA, LOX, HILPDA and low expression of immune-related genes like TMEM173(STING)

# Generate Figures 5D-F, S3C-D, S4C-D, S5B
# Example: In Sample, the perinecrotic niche was niche 4
# Set DefaultAssay to “SCT”
DefaultAssay(xenium.obj) <- "SCT"
# Check cell identities of xenium.obj
Idents(xenium.obj)
# If needed, change cell identities to “niches”
Idents(object = xenium.obj) <- "niches"
# Confirm cell identities are “niches”
Idents(xenium.obj)

# Use ImageDimPlot to make figures

# Perinecrotic niche
ImageDimPlot(xenium.obj, fov = "fov2", nmols = 20000, cols = c("lightgrey", "lightgrey", "lightgrey", "#D6272B", "lightgrey"), size = 1, alpha = 0.3)

# VEGFA
ImageDimPlot(xenium.obj, fov = "fov2", molecules = c("VEGFA"), mols.cols = "blue", mols.size = 0.25, nmols = 20000, cols = c("lightgrey", "lightgrey", "lightgrey", "#D6272B", "lightgrey"), size = 1.0, alpha = 0.3)

# TMEM173
ImageDimPlot(xenium.obj, fov = "fov2", molecules = c("TMEM173"), mols.cols = "#7CAE00", mols.size = 0.25, nmols = 20000, cols = c("lightgrey", "lightgrey", "lightgrey", "#D6272B", "lightgrey"), size = 1.0, alpha = 0.3)

# LOX
ImageDimPlot(xenium.obj, fov = "fov2", molecules = c("LOX"), mols.cols = "#17BECF", mols.size =0.25, nmols = 20000, cols = c("lightgrey", "lightgrey", "lightgrey", "#D6272B", "lightgrey"), size = 1.0, alpha = 0.3)

# For H&E Overlay of perinecrotic niche on H&E image use this
ImageDimPlot(xenium.obj, fov = "fov2", nmols = 20000, cols = c("black", "black", "black", "#D6272B", "black"), size = 1, alpha = 0.3)

# Generate Figures 5J-L, S3G-H, S4G-H, S5D
# Generate a table with numbers of each cell type in each niche
table(xenium.obj$celltype, xenium.obj$niches)
# Store table in variable called response
response <- table(xenium.obj$celltype, xenium.obj$niches)
# Convert to data frame and store in variable response 1
response1 <- as.data.frame(response)
# Save as a CSV file
write.table(response1, sep=",", file="Sample_CellTypes_Per_Niche.csv", row.names=TRUE, col.names=NA, quote=FALSE)

# Organize CSV file: col 1A is “Niche”, col 1B is “Cell_Type”, col 1C is “Percent”. Percent is percentage of each cell type in each niche the total of which is 100%.

# Load new CSV file to make figure of cell composition, store in variable data_summary
data_summary <- read_csv("Sample Sample Cell Composition To Graph CSV.csv")

# Make cell composition figure
ggplot(data_summary, aes(x = factor(Niche), y = Percent, fill = Cell_Type, colour = "black")) + 
  geom_bar(stat = "identity", color = "black") + theme_classic() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
