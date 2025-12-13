#libraries####
library(data.table)
library(Biobase)
library(dplyr)
library(textshape)
library(tibble)
library(readxl)
library(GSVA)
library(GSEABase)

#Gene list for subtypes, from Wang 2017####
gbm_subtype_gene_sets <- list(
  Mesenchymal = c("ARPC1B","S100A11","CTSC","GLIPR1","NNMT","VDR","RGS2","CTSB","TGFBI","PLAUR","LY96","BCL3","TNFAIP8","IER3","PRSS23","IL7R","RAB27A","RUNX1","P4HA2","CYP1B1","BACE2","ACPP","FTL","SLPI","RAC2","RARRES1","SYNGR2","THBS1","IL6","CAV1","PI3","CDCP1","ITGB1","LOX","CD72","COL1A2","ANPEP","MMP7","SPAG4","BNC2","NDRG1","CNN2","LUM","PTGS2","COL3A1","COL5A1","SDC1","COL1A1","GPRC5A","COL15A1"),
  Proneural = c("HN1", "RAB33A", "HDAC2", "MYT1", "MTSS1", "HOXD3", "GPR17", "PTTG1", "KLRC3", "HRASLS", "TCP1", "NPPA", "PFDN2", "CA10", "EPHB1", "UGT8", "PAK7", "SLC1A1", "NARF", "DCTN3", "SMPD3", "ZNF804A", "RASL11B", "MYB", "PDGFRA", "ERBB3", "CLGN", "SOX10", "BCL11A", "NMU", "ZNF643", "CDKN1C", "JPH3", "PCDHA9", "IL1RAPL1", "MAST1", "VIPR2", "SIM2", "BAMBI", "PKMYT1", "PLCB4", "SLC17A6", "KLRK1", "CENPJ", "NHLH1", "GABRB3", "KLRC4", "KCNK3", "GRID2", "DACH1"),
  Classical = c("PTPRA", "ELOVL2", "MLC1", "SOX9", "ARNTL", "DENND2A", "BBS1", "ABLIM1", "PAX6", "ZHX3", "USP8", "PLCG1", "CDH4", "RASGRP1", "ACSBG1", "CST3", "BCKDHB", "LHFP", "VAV3", "ACSL3", "EYA2", "SEPT11", "SLC4A4", "SLC20A2", "C14orf159", "CTNND1", "ZFHX4", "SPRY2", "ZNF45", "NCOA1", "PLCE1", "DTNA", "POLRMT", "SALL1", "TYK2", "TJP1", "MEOX2", "FGFR3", "STXBP3", "GRIK1", "GATM", "UPF1", "NPEPL1", "KIAA0494", "RBCK1", "PHKB", "SLC3A2", "PPARGC1A", "PNPLA6", "MYO5C")
)

#TCGA ####
setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Expression Data")
tcga <- as.data.frame(read_excel("TCGA Expression No Dupes.xlsx"))
rownames(tcga) <- tcga[[1]]                 # Use first column as rownames
tcga <- tcga[ , -1]   
tcga_clean <- apply(tcga, 2, as.numeric)
rownames(tcga_clean) <- rownames(tcga)

tcga_clean <- as.matrix(tcga_clean)

ssgsea_param <- ssgseaParam(
  exprData = tcga_clean,     
  geneSets = gbm_subtype_gene_sets,     
  alpha = 0.25,                        
  normalize = TRUE                      
)

ssgsea_scores <- gsva(ssgsea_param)

tcga_subtype <- as.data.frame(apply(ssgsea_scores, 2, function(x) names(which.max(x))))

#TEMPUS####
setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Expression Data")
tempus <- as.data.frame(fread("TEMPUS RNAseq expression.csv"))
colnames(tempus) <- as.character(tempus[1, ])   # Use row 1 as colnames
tempus <- tempus[-1, ]    
rownames(tempus) <- tempus[[1]]                 # Use first column as rownames
tempus <- tempus[ , -1]   
tempus <- log2(tempus + 1)
tempus_clean <- apply(tempus, 2, as.numeric)
rownames(tempus_clean) <- rownames(tempus)

tempus_clean <- as.matrix(tempus_clean)
gbm_subtype_filtered <- lapply(gbm_subtype_gene_sets, function(g) intersect(g, rownames(tempus_clean)))
gbm_subtype_trimmed <- lapply(gbm_subtype_filtered, function(genes) head(genes, 50))

ssgsea_param <- ssgseaParam(
  exprData = tempus_clean,     
  geneSets = gbm_subtype_gene_sets,     
  alpha = 0.25,                        
  normalize = TRUE                      
)

ssgsea_scores <- gsva(ssgsea_param)

tempus_subtype <- as.data.frame(apply(ssgsea_scores, 2, function(x) names(which.max(x))))
