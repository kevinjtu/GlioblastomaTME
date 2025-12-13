#libraries####
library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(RColorBrewer)

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
  )

# sanity check: the non‐tumor fractions should now sum to 1
data2 %>%
  mutate(check = rowSums(across(`Oligodendrocyte`:`Radial_Glia`))) %>%
  dplyr::select(sample_id, tumor_purity, check) %>%
  head(10)

df <- as.data.frame(data2)

# make plots ####
#top plot, tumor purity
# DEFINE CELL TYPE GROUPS
tumor_cols <- c("OPC_like_Cancer", "AC_like_Cancer", "NPC_like_Cancer", "MES_like_Cancer")

# IMPORTANT: Here we manually specify the biologically ordered tme_cols list:
tme_cols <- c(
  # T Cells (from lightest to darkest green)
  "T_Cells_Stress_Signature", 
  "Proliferating_T_Cells", 
  "CD8_Effector_Memory_T_Cells",
  "NK_like_CD8_T_Cells", 
  "Regulatory_T_Cells", 
  "Cytotoxic_CD8_T_Cell",
  "Resting_CD4_T_Cells", 
  "CD4_T_Cell_IFN_Signature",
  
  # B Cells (cyan shades)
  "B_Cell", 
  "Plasma_Cell",
  
  # Myeloid / Innate (purple shades)
  "BDM_Hypoxia_and_MES_like", 
  "BDM_MHC_Expression", 
  "BDM_IFN_Signature",
  "Anti_inflammatory_BDM", 
  "Anti_inflammatory_Monocyte", 
  "Hypoxia_associated_Monocytes", 
  "Naive_Monocyte",
  "Microglia_Aging_Signature", 
  "Pro_inflammatory_Microglia", 
  "Proliferative_Microglia",
  "Dendritic_like_Cell", 
  "Conventional_DC_1", 
  "Conventional_DC_2", 
  "Plasmacytoid_DC",
  
  # Natural Killer
  "Natural_Killer_Cells", 
  
  # Mast Cells
  "Mast_Cell",
  
  # Glia / Neural (blue shades)
  "Oligodendrocyte", 
  "Astrocyte", 
  "OPC", 
  "Radial_Glia", 
  "Neuron", 
  "Pericyte",
  
  # Vasculature (grayscale)
  "Arterial_Vessel", 
  "Tiplike_Vessel", 
  "Capilary_Vessel"
)

# CALCULATE AVERAGE TUMOR CELLULARITY PER STUDY (for ordering)
tumor_avg <- df %>%
  dplyr::select(study, all_of(tumor_cols)) %>%
  group_by(study) %>%
  summarise(across(everything(), mean, .names = "mean_{.col}")) %>%
  mutate(total_tumor = rowSums(across(starts_with("mean_"))))

# Get ordered study levels based on increasing tumor cellularity
study_order <- tumor_avg %>%
  arrange(total_tumor) %>%
  pull(study)

# Pivot tumor data long for plotting
tumor_long <- tumor_avg %>%
  dplyr::select(-total_tumor) %>%
  pivot_longer(
    cols = -study,
    names_to = "cell_type",
    values_to = "proportion",
    names_prefix = "mean_"
  ) %>%
  mutate(study = factor(study, levels = study_order))

# --- Bottom Plot Data: TME Composition ---
tme_avg <- df %>%
  dplyr::select(study, all_of(tme_cols)) %>%
  group_by(study) %>%
  summarise(across(everything(), mean, .names = "mean_{.col}")) %>%
  pivot_longer(
    cols = -study,
    names_to = "cell_type",
    values_to = "proportion",
    names_prefix = "mean_"
  ) %>%
  mutate(
    study = factor(study, levels = study_order), # match study order
    cell_type = factor(cell_type, levels = tme_cols) # enforce biologically meaningful order
  )

# COMPREHENSIVE TME COLOR MAP
tme_color_map <- c(
  "T_Cells_Stress_Signature"     = "#e5f5e0",
  "Proliferating_T_Cells"        = "#c7e9c0",
  "CD8_Effector_Memory_T_Cells"  = "#a1d99b",
  "NK_like_CD8_T_Cells"          = "#74c476",
  "Regulatory_T_Cells"           = "#41ab5d",
  "Cytotoxic_CD8_T_Cell"         = "#238b45",
  "Resting_CD4_T_Cells"          = "#006d2c",
  "CD4_T_Cell_IFN_Signature"     = "#00441b",

  "B_Cell"                       = "#80cdc1",
  "Plasma_Cell"                  = "#35978f",
  
  "BDM_Hypoxia_and_MES_like"     = "#dadaeb",
  "BDM_MHC_Expression"           = "#bcbddc",
  "BDM_IFN_Signature"            = "#9e9ac8",
  "Anti_inflammatory_BDM"        = "#807dba",
  "Anti_inflammatory_Monocyte"   = "#6a51a3",
  "Hypoxia_associated_Monocytes" = "#54278f",
  "Naive_Monocyte"               = "#3f007d",
  
  "Microglia_Aging_Signature"    = "#fcbba1",
  "Pro_inflammatory_Microglia"   = "#fc9272",
  "Proliferative_Microglia"      = "#fb6a4a",
  
  "Dendritic_like_Cell"          = "#ef3b2c",
  "Conventional_DC_1"            = "#cb181d",
  "Conventional_DC_2"            = "#a50f15",
  "Plasmacytoid_DC"              = "#67000d",
  
  "Natural_Killer_Cells"         = "#fa9fb5",
  "Mast_Cell"                    = "#dd3497",
  
  "Oligodendrocyte"              = "#6baed6",
  "Astrocyte"                    = "#4292c6",
  "OPC"                          = "#2171b5",
  "Radial_Glia"                  = "#08519c",
  "Neuron"                       = "#08306b",
  
  "Pericyte"                     = "#d9d9d9",
  "Arterial_Vessel"              = "#bdbdbd",
  "Tiplike_Vessel"               = "#737373",
  "Capilary_Vessel"              = "#525252"
)

# Tumor color map (simple version)
tumor_color_map <- c(
  "OPC_like_Cancer" = "#f6e8c3",
  "AC_like_Cancer"  = "#dfc27d",
  "NPC_like_Cancer" = "#bf812d",
  "MES_like_Cancer" = "#8c510a"
)

# PLOTS
# Subset datasets by platform
df_rnaseq <- df #%>% filter(platform == "rnaseq")
#df_microarray <- df %>% filter(platform == "microarray")

# Function to generate combined tumor + TME plots per platform
generate_platform_plot <- function(df_platform, platform_label) {
  
  # Calculate tumor averages and order studies by tumor cellularity
  tumor_avg <- df_platform %>%
    dplyr::select(study, all_of(tumor_cols)) %>%
    group_by(study) %>%
    summarise(across(everything(), mean, .names = "mean_{.col}")) %>%
    mutate(total_tumor = rowSums(across(starts_with("mean_"))))
  
  study_order <- tumor_avg %>%
    arrange(total_tumor) %>%
    pull(study)
  
  tumor_long <- tumor_avg %>%
    dplyr::select(-total_tumor) %>%
    pivot_longer(
      cols = -study,
      names_to = "cell_type",
      values_to = "proportion",
      names_prefix = "mean_"
    ) %>%
    mutate(study = factor(study, levels = study_order))
  
  # TME data processing
  tme_avg <- df_platform %>%
    dplyr::select(study, all_of(tme_cols)) %>%
    group_by(study) %>%
    summarise(across(everything(), mean, .names = "mean_{.col}")) %>%
    pivot_longer(
      cols = -study,
      names_to = "cell_type",
      values_to = "proportion",
      names_prefix = "mean_"
    ) %>%
    mutate(
      study = factor(study, levels = study_order),
      cell_type = factor(cell_type, levels = tme_cols)
    )
  
  # Tumor plot
  plot_tumor <- ggplot(tumor_long, aes(x = study, y = proportion, fill = cell_type)) +
    geom_bar(stat = "identity", position = "stack", color = "black") +  # add black outline
    scale_fill_manual(values = tumor_color_map, name = "Tumor Subtype") +
    labs(
      title = paste0(platform_label),
      x = NULL, y = "Tumor Cellularity"
    ) +
    theme_classic() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "right"
    ) +
    scale_y_continuous(labels = scales::percent)
  
  # TME plot
  plot_tme <- ggplot(tme_avg, aes(x = study, y = proportion, fill = cell_type)) +
    geom_bar(stat = "identity", position = "stack", color = "black") +  # add black outline
    scale_fill_manual(values = tme_color_map) +
    labs(
      title = NULL,
      x = NULL, y = "TME Composition"
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.x = element_text(angle = 90, hjust = 1),
      legend.position = "right",
    ) +
    guides(fill = guide_legend(ncol = 4)) +
    scale_y_continuous(labels = scales::percent)
  
  # Combine tumor and TME plots vertically
  combined_plot <- plot_tumor / plot_tme + plot_layout(heights = c(0.2, 1))
  return(combined_plot)
}

# Generate plots for each platform
plot_rnaseq <- generate_platform_plot(df_rnaseq, "RNA-seq")
#plot_microarray <- generate_platform_plot(df_microarray, "Microarray")

# Display plots separately
print(plot_rnaseq)
#print(plot_microarray)