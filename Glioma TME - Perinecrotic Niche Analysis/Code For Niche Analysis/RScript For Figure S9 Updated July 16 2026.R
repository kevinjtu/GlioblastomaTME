## Figure S9
## Purpose of Script: Data from Xenium samples was stored in a Seurat object 
## and saved as an R object. Xenium samples were loaded into R and analyzed. 
## Cell types were removed from each sample so that each sample only contained
## myeloid cells. Seurat's FindMarkers function was used to conduct differential
## gene expression analysis comparing gene expression in the perinecrotic niche
## to the non-necrotic niches in each sample. Barplots were generated to show
## Log2FoldChange of perinecrotic niche markers (VEGF, LOX) and TMEM173 in the
## perinecrotic niche compared to the non-necrotic niches in each sample.
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
library(tidyverse) # version 2.0.0
library(sf) # version 1.0-17
library(scCustomize) # version 2.1.2
library(scclusteval) # version 0.0.0.9000
library(ggplot2) # version 3.5.1
library(patchwork) # version 1.3.0
library(pals) # version 1.9
library(gridExtra) # version 2.3

# Set working directory
setwd("/home/gd2")

# Initial Processing Of Xenium object to remove all cell types but myeloid cells 
# Example here is Sample 3P
# Read In Xenium file for Sample 3P. This is a Seurat file stored as an R object
Xenium.obj <- readRDS("/home/gd2/Sample.rds")
# Look at Seurat object contents
Xenium.obj

# Eliminate All Cell Types But Myeloid Cells 
# Look at cell identities of Xenium.obj
Idents(Xenium.obj) 
# Change cell identities to "celltype"
Idents(object = Xenium.obj) <- "celltype"
# Confirm new cell identities are "celltype"
Idents(Xenium.obj)
# Remove non-myeloid cell types from Xenium.obj
# Store object in variable called Male
Male <- subset(Xenium.obj, subset = celltype != "AC/Oligo-like_1")
Male <- subset(Male, subset = celltype != "AC/Oligo-like_2")
Male <- subset(Male, subset = celltype != "Neuronal")
Male <- subset(Male, subset = celltype != "Oligodendrocyte")
Male <- subset(Male, subset = celltype != "Tumor_1")
Male <- subset(Male, subset = celltype != "Tumor_AC/OPC_like_1")
Male <- subset(Male, subset = celltype != "Tumor_AC/OPC_like_2")
Male <- subset(Male, subset = celltype != "Tumor_AC/OPC_like_3")
Male <- subset(Male, subset = celltype != "Tumor_AC/OPC_like_4")
Male <- subset(Male, subset = celltype != "Tumor_AC/OPC_like_5")
Male <- subset(Male, subset = celltype != "Tumor_OPC-like_1")
Male <- subset(Male, subset = celltype != "Tumor_OPC-like_2")
Male <- subset(Male, subset = celltype != "Unknown_1")
Male <- subset(Male, subset = celltype != "Unknown_2")
Male <- subset(Male, subset = celltype != "Vascular")
# Look at contents in Male variable
Male

# Confirm Only Myeloid Cells Present In Male variable
# Do this by looking at cell types in metadata
h <- (Male@meta.data) # get metadata from Male and store in variable h
head(h) # look at first 6 rows in variable h
avg.mon <- as.data.frame(h) # convert to data frame
# Save metadata info as a CSV file
write.table(avg.mon, sep=",", file = "Sample3P_Myeloid_Only_metadata.csv", row.names=T)

# Save new Xenium object with myeloid cells only
saveRDS(Male,"Sample_3P_With_MyeloidCells_Only.Robj")

# Look at DEGs in Perinecrotic Niche 1 versus non-necrotic niches
DefaultAssay(Male) <- "SCT" # set default assay to "SCT"
Idents(Male) # check cell identities of xenium object
Idents(object = Male) <- "niches" # change cell identities to niches
Idents(Male) # confirm cell identities are now niches
# Use FindMarkers function to look at DEGs, store output in variable Niche1.markers
Niche1.markers <- FindMarkers(Male, only.pos = FALSE, ident.1 = "1", ident.2 = NULL)
head(Niche1.markers) # look at first 6 rows of data
avg.mon <- as.data.frame(Niche1.markers) # convert to data frame
# Save DEG info as a CSV file
write.table(avg.mon, sep=",", file = "Sample3P_Myeloid_PN_DEGs.csv", row.names=TRUE, col.names = NA, quote=FALSE)

# Load CSV file. 1A is “Sample”, 1B is “Gene”, 1C is “Log2FC”
# Store in variable called data_summary
data_summary <- read_csv("TMEM173 Myeloid (To make graph).csv")

# Define the desired order of genes of interest
desired_order <- c("VEGFA", "LOX", "TMEM173")

# Reorder the factor levels in your data frame
data_summary$Gene <- factor(data_summary$Gene, levels = desired_order)

# Generate Figure S14 as a bar plot
ggplot(data_summary, aes(x = factor(Sample),
                         y = Log2FC,
                         fill = Gene)) +
  geom_bar(stat = "identity", color = "black",
           position = position_dodge(width=0.9)) + theme_classic() + ylim(-2, 4) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.y = element_text(size = 28, face = "bold"), axis.line.x = element_line(size = 2), axis.line.y = element_line(size = 2),
        axis.ticks = element_line(size = 2), axis.ticks.length = unit(0.25, "cm"))
