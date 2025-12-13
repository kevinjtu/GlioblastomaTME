#libraries####
library(data.table)
library(dplyr)
library(tidyr)
library(tibble)
library(msigdbr)
library(ggplot2)
library(ggpubr)
library(scales)
library(compositions)
library(mclust)

#import & clean data####
subtype_cols <- c("OPC_like_Cancer", "AC_like_Cancer", "NPC_like_Cancer", "MES_like_Cancer")
setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME")
data <- as.data.frame(fread("GBM Clinical-TME Cell States.csv")) %>%
  mutate(subtype = subtype_cols[
    max.col( as.matrix( across(all_of(subtype_cols)) ),
             ties.method = "first")])

#get TME fractions accounting for cellularity
# define your six tumor‐cell columns
tumor_cells <- c(
  "OPC_like_Cancer",
  "AC_like_Cancer",
  "NPC_like_Cancer",
  "MES_like_Cancer"
)

data2 <- data %>%
  # 1) compute tumor purity as the row‐sum of those six columns
  rowwise() %>%
  mutate(tumor_purity = sum(c_across(all_of(tumor_cells)))) %>%
  ungroup() %>%
  
  # 2) renormalize all non‐tumor cell fractions so they sum to 1
  #    by dividing each by (1 – tumor_purity)
  mutate(
    across(
      `Oligodendrocyte`:`Radial_Glia`,
      ~ . / (1 - tumor_purity)
    )
  ) %>%
  
  # 3) drop the original tumor‐cell columns
  dplyr::select(-all_of(tumor_cells))

# sanity check: the non‐tumor fractions should now sum to 1
data2 %>%
  mutate(check = rowSums(across(`Oligodendrocyte`:`Radial_Glia`))) %>%
  dplyr::select(sample_id, tumor_purity, check) %>%
  head(10)

data2 <- as.data.frame(data2)

df <- data2 %>%
  filter(recurrence == "primary")

cell_cols <- names(df)[ which(names(df) == "Oligodendrocyte"):
                          which(names(df) == "Radial_Glia") ]

#CLR Transformation####
# 1) extract the raw proportions matrix
df_cells <- df %>% dplyr::select(all_of(cell_cols))

# 2) add a tiny pseudocount to avoid zeros
#pcount <- min(df_cells[df_cells > 0], na.rm = TRUE) * 0.5
#df_cells[df_cells == 0] <- pcount

# 3) compute the CLR
clr_mat <- clr(as.matrix(df_cells))

# 4) rename columns to mark them as CLR‐transformed#
colnames(clr_mat) <- paste0(cell_cols, "_clr")

# 5) rebuild your analysis data frame, replacing raw with CLR
df_clr <- df %>%
  dplyr::select(-all_of(cell_cols)) %>%
  bind_cols(as.data.frame(clr_mat))

# 6) update the list of columns you’ll loop over
cell_cols_clr <- paste0(cell_cols, "_clr")


#separate samples by low & high STING####
clusters <- Mclust(df_clr$STING1_tpm, G = 2, modelNames = "V")
df_clr$cluster <- factor(
  clusters$classification,
  levels = c(1, 2),
  labels = c("Low STING", "High STING")
)

df_clr <- df_clr %>%
  filter(!is.na(location)) %>%
  mutate(
    cluster = ntile(STING1_tpm, 4),
    cluster = factor(
      cluster,
      levels = 1:4,
      labels = c("Q1 (lowest)", "Q2", "Q3", "Q4 (highest)")
    )
  )

# Quick check
table(df_clr$cluster)

# 2. Define the ordered region levels
region_levels <- c(
  "Tumor Body",
  "Hyperplastic blood vessels in cellular tumor",
  "Microvascular Proliferation",
  "Tumor Edge"
)

# 3. Map each location into those four regions
df_region <- df_clr %>%
  mutate(
    region = case_when(
      location %in% c("Perinecrotic Zone", 
                      "Pseudopalisading Cells", 
                      "Cellular Tumor") ~ "Tumor Body",
      location == "Hyperplastic blood vessels in cellular tumor" ~
        "Hyperplastic blood vessels in cellular tumor",
      location == "Microvascular Proliferation" ~
        "Microvascular Proliferation",
      location %in% c("Infiltrating Tumor", "Leading Edge") ~ "Tumor Edge",
      TRUE ~ NA_character_
    ),
    region = factor(region, levels = region_levels)
  ) %>%
  filter(!is.na(region))

# 4. Pivot the CLR‐abundance columns to long form
cell_cols <- df_clr %>%
  select(Oligodendrocyte_clr:Radial_Glia_clr) %>%
  names()

df_long <- df_region %>%
  pivot_longer(
    cols      = all_of(cell_cols),
    names_to  = "cell",
    values_to = "abundance"
  )

# 5. Compute mean and SEM by region, quartile‐cluster, and cell
df_summary <- df_long %>%
  group_by(region, cluster, cell) %>%
  summarise(
    mean_abundance = mean(abundance, na.rm = TRUE),
    sem            = sd(abundance,   na.rm = TRUE) / sqrt(n()),
    .groups        = "drop"
  )

# 6. Plot the four regions for each STING quartile
ggplot(df_summary, aes(
  x     = region,
  y     = mean_abundance,
  color = cluster,
  group = cluster
)) +
  geom_line() +
  geom_point(position = position_dodge(width = 0.2), size = 2) +
  geom_errorbar(
    aes(ymin = mean_abundance - sem, ymax = mean_abundance + sem),
    width    = 0.1,
    position = position_dodge(width = 0.2)
  ) +
  facet_wrap(~ cell, scales = "free_y", ncol = 5) +
  labs(
    x     = "Anatomic Region",
    y     = "Mean CLR Abundance",
    color = "STING₁ TPM Quartile",
    title = "Cell Abundance Across Regions by STING₁ TPM Quartile"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text  = element_text(size = 7)
  )