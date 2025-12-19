## Figure S6
## Purpose Of Script: Conduct Banksy Analysis Of Each Xenium Sample. Determine If Banksy Identified A Perinecrotic Niche (PN) In The Sample. Determine If The Niche Was Identical To The PN Identified By Seurat’s BuildNicheAssay. 
## R version 4.2.2. RStudio version 2022.12.03.353.
##
## Author: George D. Dalton
##
## Date: December 15, 2025

ls() ## list objects currently present in workspace
rm(list=ls()) ## remove objects from current workspace
ls() ## confirm current workspace is empty

## libraries loaded during session
library(Seurat) # version 5.1.0
library(SeuratData) # version 0.2.2.9001
library(SeuratWrappers) # version 0.3.2
library(Banksy) # version 0.1.6
library(ggplot2) # version 3.5.1
library(pals) # version 1.9
library(gridExtra) # version 2.3
library(SpatialExperiment) # version 1.8.1
library(SummarizedExperiment) # version 1.28.0
library(scuttle) # version 1.8.4
library(scater) # version 1.26.1
library(cowplot) # version 1.1.3

# Set working directory 
setwd("/home/gd2")

# Step 1 Read In Xenium file (example name: Sample) that is a Seurat file stored as an R object 
Xenium.obj <- readRDS("~/Sample.rds")
# Look at Seurat object
Xenium.obj # for this example, fov = fov2

# Step 2 Extract expression data from Xenium.obj and store as variable gcm
gcm <- Xenium.obj@assays$SCT@data
# Get dimensions of gcm
dim(gcm)
# Look at first 6 rows of data in gcm
head(gcm)
# gcm should be a matrix
class(gcm)

# Step 3 Get spatial centroid info from spatial object in Xenium.obj and store as csv file
y <- (Xenium.obj@images$fov2)
z <- GetTissueCoordinates(y)
# Save centroid information as a CSV file
write.table(z, sep=",", file = "TissueCoords_All.csv", row.names=T)

# Step 4 Read in data frame with centroid spatial info from file “TissueCoords_All.csv”
# row.names: Col 1A = “cell”, Col 1B = “x”, Col 1C = “y”
res <-read.table(file.choose(), sep = ",", header = TRUE, row.names = 1)
class(res) # res should be a data frame
# Look at first 6 rows of data
head(res)
# Look at dimensions of gcm
dim(res)
is.atomic(res) # this should be false

# Run BANKSY package
bank <- BanksyObject(own.expr = gcm, cell.locs = res)
# Look at Banksy object contents
bank
# will see something like below
# Object of class BanksyObject 
# Assay with 358059 cells 354 features
# Spatial dimensions: sdimx sdimy 
# Metadata names: cell_ID nCount NODG 
# Dimension reductions:

# Continue to run BANKSY
# Normalize
bank <- NormalizeBanksy(bank)
# Scale data
bank <- ScaleBanksy(bank)
# Get BANKSY matrix
bank <- ComputeBanksy(bank, k_geom = 50) 
# Run PCA
bank <- RunBanksyPCA(bank, use_agf = FALSE, lambda = 0.8, npcs = 20, verbose = TRUE)
# Run UMAP
bank <- RunBanksyUMAP(bank, use_agf = FALSE, lambda = 0.8, npcs = 20, verbose = TRUE) 
# Get BANKSY clusters, note: this step can take time
bank <- ClusterBanksy(bank, use_agf = FALSE, lambda = 0.8, npcs = 20, k.neighbors = 50, resolution = 0.4)
# Look at final Banksy object and the info it now contains
bank

# Save BANKSY object as an R object 
saveRDS(bank,"Sample_BANKSY.Robj")

# Take a look at BANKSY clusters on spatial map
plotSpatial(bank, by = 'clust_M0_lam0.8_k50_res0.4', type = 'discrete')

# Take a look at BANKSY clusters on a UMAP plot
plotReduction(bank, by = 'clust_M0_lam0.8_k50_res0.4', type = 'discrete')

# Get BANKSY object metadata and store as variable h
h <- meta.data(bank)
# Convert to data frame stored as variable banksy.meta
banksy.meta <- as.data.frame(h)
# Save data frame as a CSV file
write.table(banksy.meta, sep=",", file = "Sample_BANKSY_Metadata.csv", 
            row.names=TRUE, col.names = NA, quote=FALSE)

# Read in CSV file with BANKSY metadata, column with clusters is clust_M0_lam0.8_k50_res0.4
banksy.meta2 <-read.table(file.choose(), sep = ",", header = TRUE, row.names = 1)
# Look at column names in CSV file
head(banksy.meta2)
# Add BANKSY cluster data To Xenium.obj metadata with column name BanksyCluster
Xenium.obj@meta.data$BanksyCluster <- banksy.meta2$clust_M0_lam0.8_k50_res0.4
# Confirm Xenium.obj has metadata column with BANKSY clusters called BanksyCluster
head(Xenium.obj@meta.data)
# Store Xenium.obj metadata as variable called how
how <- Xenium.obj@meta.data
# Convert how to a data frame stored as how1
how1 <- as.data.frame(how)
# save the new Xenium.obj metadata as a CSV file
write.table(how1, sep=",", file = "Sample_Metadata_With_BANKSY_Clusters.csv", 
            row.names=TRUE, col.names = NA, quote=FALSE)

# Save Updated Xenium.obj as an R object
saveRDS(Xenium.obj,"Sample_With_BANKSY_Info.rds")

# Look at expression of each perinecrotic niche marker gene (e.g., LOX, HILPDA, TMEM173, VEGFA) in each BANSKY identified niche
DefaultAssay(Xenium.obj) <- "SCT"
Idents(Xenium.obj)
Idents(object = Xenium.obj) <- "BanksyCluster"
Idents(Xenium.obj)
DotPlot(Xenium.obj, features = "HILPDA") + RotatedAxis()

# Can also look at perinecrotic niche marker genes (LOX, HILPDA, VEGFA, TMEM173) simultaneously in each BANKSY identified niche 
# Define variable "markers.to.plot" to look at specific markers
markers.to.plot <- c("LOX", "VEGFA", "HILPDA", "TMEM173")
# Plot
DotPlot(Xenium.obj, features = markers.to.plot, cols = c("gray", "tomato"), dot.scale = 8) +
  RotatedAxis()

# If plots suggest BANKSY niche 2 is perinecrotic niche that was also identified with Seurat’s BuildNicheAssay
# Look at DEGs in BANKSY niche 2 compared to other BANKSY niches
# Do this using Seurat’s FindMarkers function
DefaultAssay(Xenium.obj) <- "SCT"
Idents(Xenium.obj)
Idents(object = Xenium.obj) <- "BanksyCluster"
Idents(Xenium.obj)
Niche2.markers <- FindMarkers(Xenium.obj, only.pos = FALSE, ident.1 = "2", ident.2 = NULL)
head(Niche2.markers)
# convert to data frame
avg.mon <- as.data.frame(Niche2.markers)
# save DEG info as a CSV file, DEGs will have padjust < 0.05 and Log2FC < -0.5 or > 0.5
write.table(Niche2.markers, sep=",", file = "Niche2Markers_All.csv", row.names=T)
# Confirm DEGs in Banksy niche 2 have elevated HILPDA, VEGFA, LOX and low TMEM173, etc.

# To make Figure S6 BANKSY plot 
# Determined BANKSY niche 2 is perinecrotic niche identified using Seurat’s BuildNicheAssay
# Perinecrotic niche should look same on spatial plots of BuildNicheAssay and Banksy niches
# Make spatial plot of Banksy perinecrotic niche
# Subset Xenium.obj to make an object with only BANKSY niche 2
Two_Only <- subset (Xenium.obj, subset = BanksyCluster == "2")
Idents(Two_Only)
Idents(object = Two_Only) <- "BanksyCluster"
Idents(Two_Only)
# Plot this
niche.plot <- ImageDimPlot(Two_Only, group.by = "BanksyCluster", size = 1.5, dark.background = F) + ggtitle("Niches")
niche.plot