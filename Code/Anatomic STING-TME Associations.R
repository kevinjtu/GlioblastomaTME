#libraries####
library(data.table)
library(dplyr)
library(ggplot2)
library(readxl)
library(ggpubr)
library(mclust)
library(stringr)
library(reshape2)
library(effsize)
library(tidyverse)
library(broom)

#input data####
#input TME data
setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Instaprism Results")
ivygap_tme <- read_excel("IVYGAP RNAseq InstaPrism results.xlsx", sheet = 4)
colnames(ivygap_tme)[1] <- "sample_id"

opc_like_columns <- grep("OPC-like", names(ivygap_tme), value = TRUE)
ac_like_columns <- grep("AC-like", names(ivygap_tme), value = TRUE)
npc_like_columns <- grep("NPC-like", names(ivygap_tme), value = TRUE)
mes_like_columns <- grep("MES-like", names(ivygap_tme), value = TRUE)
tcell_columns <- grep("CD4|CD8", names(ivygap_tme), value = TRUE)
dc_columns <- grep("DC", names(ivygap_tme), value = TRUE)
TAM_columns <- grep("TAM", names(ivygap_tme), value = TRUE)

ivygap_tme$OPC_like <- rowSums(ivygap_tme[opc_like_columns], na.rm = TRUE)
ivygap_tme$AC_like <- rowSums(ivygap_tme[ac_like_columns], na.rm = TRUE)
ivygap_tme$NPC_like <- rowSums(ivygap_tme[npc_like_columns], na.rm = TRUE)
ivygap_tme$MES_like <- rowSums(ivygap_tme[mes_like_columns], na.rm = TRUE)
ivygap_tme$Tcell <- rowSums(ivygap_tme[tcell_columns], na.rm = TRUE)
ivygap_tme$DendriticCell <- rowSums(ivygap_tme[dc_columns], na.rm = TRUE)
ivygap_tme$TAM <- rowSums(ivygap_tme[TAM_columns], na.rm = TRUE)

ivygap_tme <- ivygap_tme[, !(names(ivygap_tme) %in% c(opc_like_columns, ac_like_columns, npc_like_columns, mes_like_columns, tcell_columns, dc_columns, TAM_columns))]
rowSums(ivygap_tme[, 2:27], na.rm = TRUE) #check the TME still adds up to 1

#input anatomic tumor locations
setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/IVYGAP")
ivygap_anatomic <- read_excel("columns-samples.xlsx")
ivygap_anatomic <- ivygap_anatomic %>%
  dplyr::rename("sample_id" = "rna_well_id")

#fnd samples with high expression of STING (sting categoris 1-3, with 3 being the heighest)
setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Expression Data")
ivygap_expression <- fread("IVYGAP.csv")
setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/IVYGAP")
genes <- fread("rows-genes.csv")
genes <- genes[,c(1,4)]
ivygap_expression[,1] <- genes$gene_symbol
colnames(ivygap_expression)[1] <- "sample_id"
ivygap_expression <- as.data.frame(t(ivygap_expression))
ivygap_expression <- rownames_to_column(ivygap_expression, var = "sample_id")
colnames(ivygap_expression) <- ivygap_expression[1, ]
ivygap_expression <- ivygap_expression[-1, ]

mclust <- Mclust(as.data.frame(ivygap_expression$TMEM173)) #geneid for STING1/TMEM173 is 109790
sting_cat <- as.data.frame(mclust$classification)
sting_cat <- dplyr::rename(sting_cat, STING_group = `mclust$classification`)
sting_cat <- cbind(ivygap_expression, sting_cat)
sting_cat <- dplyr::select(sting_cat, sample_id, TMEM173, STING_group)

#marge all info together
ivygap_merged <- merge(ivygap_anatomic, sting_cat, by = "sample_id", all.x = TRUE)
ivygap_merged <- merge(ivygap_merged, ivygap_tme, by = "sample_id", all.x = TRUE)

#Compare cell types between anatomic locations for high sting samples####
#have the order go from Cellular Tumor -> Infiltrating Tumor -> Leading Edge, do not consider the other sites for now
#only have high sting samples (STING category 3)
#Facet by cell_type (columns 15:56)
anatomic_sting <- ivygap_merged %>%
  filter(STING_group == 3) %>%
  #filter(anatomic_location %in% c("Cellular Tumor", "Infiltrating Tumor", "Leading Edge")) %>%
  dplyr::select(sample_id, anatomic_location, 15:41)%>%
  melt(id.vars = c("sample_id", "anatomic_location"), variable.name = "cell_type", value.name = "cell_type_proportion")

anatomic_sting$anatomic_location <- factor(anatomic_sting$anatomic_location, levels = c(
  "Perinecrotic Zone",
  "Pseudopalisading Cells",
  "Cellular Tumor",
  "Hyperplastic blood vessels in cellular tumor",
  "Microvascular Proliferation",
  "Infiltrating Tumor",
  "Leading Edge"
))

anatomic_sting_filter <- anatomic_sting %>%
  filter(cell_type %in% c("Tcell", "TAM"))

mean_sem <- function(x) {
  mean_val <- mean(x, na.rm = TRUE)
  sem <- sd(x, na.rm = TRUE)/sqrt(length(x[!is.na(x)]))
  data.frame(y = mean_val, ymin = mean_val - sem, ymax = mean_val + sem)
}

# Create plot with error bars
ggplot(anatomic_sting_filter, 
       aes(x = anatomic_location, 
           y = cell_type_proportion, 
           group = cell_type, 
           color = cell_type)) +
  stat_summary(fun.data = mean_sem, geom = "errorbar", 
               width = 0.2, size = 0.8) +
  stat_summary(fun = "mean", geom = "point", size = 3, shape = 16) +
  stat_summary(fun = "mean", geom = "line", aes(group = cell_type), size = 1) +
  facet_wrap(~cell_type, scales = "free_y", nrow = 2) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Tumors with High STING Expression", 
       x = "", y = "Proportion of Tumor") +
  theme(legend.position = "none")

#do linear regressions between sting expression (TMEM173 column) and proporiton of cell types, have different regressions for each cell type
anatomic_sting <- ivygap_merged %>%
  dplyr::select(sample_id, anatomic_location, TMEM173, 16:41) %>%
  melt(id.vars = c("sample_id", "anatomic_location", "TMEM173"), variable.name = "cell_type", value.name = "cell_type_proportion") %>%
  dplyr::rename("STING" = "TMEM173") %>%
  mutate(
    STING = as.numeric(STING),
    STING = ifelse(is.na(STING), 0, STING) # Handle potential conversion NAs
  )

result <- anatomic_sting %>%
  group_by(cell_type, anatomic_location) %>%
  nest() %>%
  mutate(
    # Perform Spearman correlation test
    test = map(data, ~ cor.test(.x$cell_type_proportion, .x$STING, 
                                method = "spearman", exact = FALSE)),
    tidied = map(test, tidy)
  ) %>%
  unnest(tidied) %>%
  dplyr::select(
    cell_type,
    anatomic_location,
    rho = estimate,  # Spearman's rho correlation coefficient
    p_value = p.value
  )

final_cell_types <- result %>%
  # Filter to only the anatomic locations of interest
  filter(anatomic_location %in% c("Pseudopalisading Cells", "Infiltrating Tumor")) %>%
  # Create wide format for comparison
  pivot_wider(
    names_from = anatomic_location,
    values_from = c(rho, p_value),  # Changed from adj_p to p_value
    names_glue = "{anatomic_location}_{.value}"
  ) %>%
  # Filter for conditions using raw p-values
  filter(
    # In Cellular Tumor: non-significant OR negative slope
    (`Pseudopalisading Cells_p_value` >= 0.05 | `Pseudopalisading Cells_rho` <= 0) &
      # In Infiltrating Tumor: significant AND positive slope
      (`Infiltrating Tumor_p_value` < 0.05 & `Infiltrating Tumor_rho` > 0)
  ) %>%
  select(cell_type, contains("Pseudopalisading"), contains("Infiltrating"))
