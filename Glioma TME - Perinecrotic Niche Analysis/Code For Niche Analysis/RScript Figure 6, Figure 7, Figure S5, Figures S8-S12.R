## Figures 6, 7, S5, S8-S12
## Purpose Of Script: Use the R package CellChat to conduct cell-cell communication analysis in perinecrotic and non-necrotic niches in Xenium samples. R version 4.2.2. RStudio version 2022.12.03.353.
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
library(CellChat) # version 2.1.2
library(patchwork) # version 1.3.0
library(NMF) # version 0.28
library(ggalluvial) # version 0.12.5
library(presto) # version 1.0.0
library(ggplot2) # version 3.5.1
library(tidyverse) # version 2.0.0
library(forcats) # version 1.0.0
library(dplyr) # version 1.1.4

# Set working directory 
setwd("/home/gd2")

# Read in a Xenium file (example name: Sample). This is a Seurat file stored as an R object.
Xenium.obj <- readRDS("~/Sample.rds")
# Look at contents of Xenium.obj
Xenium.obj

# For CellChat analysis, create a perinecrotic and non-necrotic cellchat object
# For this sample, suppose niche 4 is perinecrotic 
# Subset Xenium.obj to create a cellchat object that only contains niche 4
Xenium.obj <- subset (Xenium.obj, subset = niches != "1") # remove niche 1
Xenium.obj <- subset (Xenium.obj, subset = niches != "2") # remove niche 2
Xenium.obj <- subset (Xenium.obj, subset = niches != "3") # remove niche 3
Xenium.obj <- subset (Xenium.obj, subset = niches != "5") # remove niche 5
# Look at Xenium.obj, should now contain only perinecrotic  niche 4
Xenium.obj

# To create a cellchat object with only the non-necrotic niches
# Subset Xenium.obj to remove perinecrotic niche 4 
Xenium.obj <- subset(Xenium.obj, subset = niches != "4") # removes niche 4
# Look at Xenium.obj that now contains only the non-necrotic niches
Xenium.obj

# Prepare input data for CellChat analysis

# Get normalized data matrix or gene expression data of spots/cells and store as data.input
# rownames are genes and colnames are cells
data.input = Seurat::GetAssayData(Xenium.obj, slot = "data", assay = "SCT") 
# Look at first 6 rows of table
head(data.input)
# Check dimensions for accuracy
dim(data.input)

# Define the metadata by manually creating a data frame consisting of the cell labels
# Store Xenium.obj metadata in variable called meta
meta <- Xenium.obj@meta.data
# Check meta dimensions
dim(meta)
# Look at first 6 rows of data in meta table
head(meta)
# Check cell identities of Xenium.obj
Idents(Xenium.obj)
# Set Xenium.obj cell identities as celltype
Idents(object = Xenium.obj) <- "celltype"
# Confirm Xenium.obj cell identities are cell types
Idents(Xenium.obj)
# Create a table: first column has cell ids, second columns title “labels” and has cell type, third column 
# titled “samples” and says Perinecrotic since this is the perinecrotic niche
meta = data.frame(labels = Seurat::Idents(Xenium.obj), samples = "Perinecrotic", row.names = colnames(data.input)) 
# Look at meta table
meta
# Run this
meta$samples <- factor(meta$samples)

# Check that all cell types or cell labels are present, this will list them
unique(meta$labels) 

# Check the sample labels
unique(meta$samples) # this will be Perinecrotic

# Recheck first 6 rows of meta table to see col 1 (cell ids), col 2 (labels), col 3 (samples) and data
head(meta)

# Get spatial centroid info from spatial object in file and store as csv file
# This Xenium sample's fov is fov9
y <- (Xenium.obj@images$fov9)
z <- GetTissueCoordinates(y)

# Look at first 6 rows of table
head(z)
# Save table as a CSV file
write.table(z, sep=",", file = "Sample_Perinecrotic.csv", row.names=T)

# Open CSV file on desktop, move cell ids to col 1, move x and y coords to col 2 and col 3, title col 2 
# imagerow, title col 3 imagecol, col 1 has no title
# Read in data frame with centroid spatial info and store in variable res
res <-read.table(file.choose(), sep = ",", header = TRUE, row.names = 1)
class(res) # should be data frame
head(res) # look at first 6 rows of table
dim(res) # check table dimensions
is.atomic(res) # should be false
# Store res in variable called spatial.locs
spatial.locs = res
# Check first 6 rows of data in spatial.locs
head(spatial.locs)
class(spatial.locs) # should be data frame

# Spatial factors of spatial Coordinates, must do this conversion data for Xenium data
conversion.factor = 1
spot.size = 10
spatial.factors = data.frame(ratio=conversion.factor, tol=spot.size/2)
head(spatial.factors)

# Create CellChat object
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels", datatype = "spatial", 
                           coordinates = spatial.locs, spatial.factors = spatial.factors)

# Look at CellChat object contents
cellchat

# Choose CellChat Database, choose human as Xenium samples are human samples
# Store in variable CellChatDB
CellChatDB <- CellChatDB.human

# Show database categories for cellchat analysis, pie chart shows up in R plot window
showDatabaseCategory(CellChatDB)

# Use All CellChatDB except for "Non-protein Signaling" for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB)

# Set the used database in the object
cellchat@DB <- CellChatDB.use

# Pre-Processing Expression Data For Cell-To-Cell Communication Analysis
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
cellchat <- identifyOverExpressedGenes(cellchat, do.fast = FALSE)
cellchat <- identifyOverExpressedInteractions(cellchat, variable.both = F)

# Compute the communication probability and infer cellular communication network
cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1,
                              distance.use = FALSE, interaction.range = 250,
                              contact.dependent = TRUE, contact.range = 100)

# Great information here. Gives a table with many communication categories
# Use table info to subdivide total interactions into communication categories for Figures 6A, S5E, S8A, S9A
df.net <- subsetCommunication(cellchat, slot.name = "net")
# Convert to data frame
avg.meta <- as.data.frame(df.net)
# Save as a CSV file.
write.table(avg.meta, sep=",", file = "Sample_Communication.csv", row.names=TRUE, col.names = NA, quote=FALSE)

# Make table with total interactions for Figures 6A, S5E, S8A, S9A
# Total interactions from this table should equal total interactions in Sample_Communication.csv file
y <- rowSums(cellchat@net$count)
# Look at first 6 rows of the table
head(y)
# Convert to a data frame and store in variable avg.meta1
avg.meta1 <- as.data.frame(y)
# Save as a CSV file. Table can be used to get total interaction for the sample.
write.table(avg.meta1, sep=",", file = "Sample_total_Interactions.csv", row.names=TRUE, col.names = NA, quote=FALSE)

# To make Figures 6A, S5E, S8A, S9A
# Load CSV file. In file, 1A is “Niche”, 1B is “Cell Type”, 1C is “Percent”. Niche is PN or NN, Cell Type is 
# signaling category, Percent is number of total interactions per category in the respective niche.
data_summary <- read_csv("Sample New Composition To Graph.csv")
# Plot data
ggplot(data_summary, aes(x = factor(Niche), y = Percent, fill = Cell_Type, colour = "black")) + 
  geom_bar(stat = "identity", color = "black") + theme_classic() + ylim(0,500) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.y = element_text(size = 28, face = "bold"), axis.line.x = element_line(size = 2), axis.line.y = element_line(size = 2),
        axis.ticks = element_line(size = 2), axis.ticks.length = unit(0.25, "cm"))

# Get table with total interaction strength for Figures 6B, S5F, S8B, S9B
# These bar plots were made in Microsoft Excel
y <- rowSums(cellchat@net$weight)
# Look at first 6 rows of the table
head(y)
# Convert to a data frame and store in variable avg.meta2
avg.meta2 <- as.data.frame(y)
# Save as a CSV file. Table can be used to get total interaction strength for sample. 
write.table(avg.meta2, sep=",", file = "Sample_Interaction_Strength.csv", row.names=TRUE, col.names = NA, quote=FALSE)

# Get table with VEGF and MHC-II signaling output and input to target cells
# Look at cell types and determine degree of outgoing connectivity in specific pathway
cellchat@netP$centr$VEGF$hub
# Look at receiver cell types and amount of input signal they receive
cellchat@netP$centr$VEGF$indeg
# Use output and input to make CSV file of signaling output to receiver cell types
# Load CSV file with signal output to specific targets cells in sample and store in variable data_summary2
data_summary2 <- read_csv("Sample_Output.csv")
# Make barplot of data for Figures 6C, 7A, S5G, S8C, S9C, S10A, S11A, S12A
ggplot(data_summary2, aes(x = factor(Niche), y = y, fill = Cell_Type, colour = "black")) + 
  geom_bar(stat = "identity", color = "black") + theme_classic() + ylim(0,20) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.y = element_text(size = 28, face = "bold"), axis.line.x = element_line(size = 2), axis.line.y = element_line(size = 2),
        axis.ticks = element_line(size = 2), axis.ticks.length = unit(0.25, "cm"))

# Make chord plots for Figures 6D,F,H; 7B,D,F; S5H; S8D,F; S9D,F; S10B,D; S11B,D; S12C
# Chord plots were made to look at specific receivers or targets (e.g. Macrophage)
pathways.show <- c("MHC-II or VEGF") 
par(mfrow=c(1,1), xpd = TRUE) 
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle",
                    vertex.label.cex = 0.05, targets.use = c("Macrophage"))
dev.off()

# Make heatmap of network centrality scores for Figures 6E,G,I; 7C,E,G; S5I; S8E,G; S9E,G; 
# S10C,E; S11C,E; S12D
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") 
par(mfrow=c(1,1))
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 20, height = 4, font.size = 10, measure.name = c("Sender", "Receiver"), measure = c("outdeg", "indeg"))

# Make dotplot of L-R communication probabilities from sender cells to receiver cells in perinecrotic &
# non-necrotic niches. Figures 6J-L; 7H-J; S5J; S8H-I; S9H-I; S10F-G; S11F-G; S12D
# Get ligand-receptor information, communication probabilities, pvalues from this table
df.net <- subsetCommunication(cellchat, slot.name = "net")
# Convert to data frame and store in variable avg.meta
avg.meta <- as.data.frame(df.net)
# Save as a CSV file.
write.table(avg.meta, sep=",", file = "Sample_LR.csv", row.names=TRUE, col.names = NA, quote=FALSE)

# Get relevant information on ligand-receptors, communication probabilities and pvalues and set up table in csv file
# Will plot informatio as a dotplot

# Read in CSV from dotplot folder
d <- read.csv("Sample MHCII.csv")
# Look at table
d
# Code for dotplot, stored in variable p
p <- d %>%
  mutate(x = fct_relevel(x, 
                         "Task1", "Task2", "Task3", 
                         "Task4", "Task5", "Task6",
                         "Task7", "Task8", "Task9",
                         "Task10", "Task11", "Task12",
                         "Task13", "Task14", "Task15",
                         "Task16", "Task17", "Task18")) %>%
  ggplot(aes(x, forcats::fct_rev(y), size = value, fill = p)) +
  geom_point(shape = 21, stroke = 0) +
  scale_x_discrete(position = "top") +
  scale_radius(range = c(1, 9)) +
  scale_fill_gradient(low = "red", high = "blue", breaks = c(0, .5, 1), labels = c("Great", "OK", "Bad"), limits = c(0, 1)) +
  theme_minimal() +
  theme(legend.position = "bottom", 
        panel.grid.major = element_blank(),axis.text = element_text(size = 12, face = "bold", color = "black"),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8)) +
  guides(size = guide_legend(override.aes = list(fill = NA, color = "black", stroke = .5), 
                             label.position = "bottom",
                             title.position = "right", 
                             order = 1),
         fill = guide_colorbar(ticks.colour = NA, title.position = "top", order = 2)) +
  labs(size = "Area = Time Spent", fill = "Score:", x = NULL, y = NULL)
# Show dotplot
p

# Save The CellChat Object
saveRDS(cellchat, file = "CellChat_Sample_Perinecrotic_Only.rds")
