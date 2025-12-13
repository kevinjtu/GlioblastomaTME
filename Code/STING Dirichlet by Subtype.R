#libraries####
library(data.table)
library(dplyr)
library(DirichletReg)
library(purrr)
library(tibble)
library(ggplot2)
library(stringr)

#Import & clean data####
subtype_cols <- c("OPC_like_Cancer", "AC_like_Cancer", "NPC_like_Cancer", "MES_like_Cancer")
setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME")
data <- as.data.frame(fread("GBM Clinical-TME Cell States.csv")) %>%
  mutate(subtype = subtype_cols[
    max.col( as.matrix( across(all_of(subtype_cols)) ),
             ties.method = "first")])

tumor_cells <- c("OPC_like_Cancer","AC_like_Cancer",
                 "NPC_like_Cancer","MES_like_Cancer")

data2 <- data %>%
  rowwise() %>%
  mutate(tumor_purity = sum(c_across(all_of(tumor_cells)))) %>%
  ungroup() %>%
  mutate(across(
    # only non-tumor up through Radial_Glia
    Oligodendrocyte:Radial_Glia, 
    ~ . / (1 - tumor_purity)
  )) %>%
  dplyr::select(-all_of(tumor_cells))

df <- data2 %>%
  filter(recurrence == "primary") %>%
  mutate(
    study_f = factor(study),
    STING_s = as.numeric(REACTOME_STING)
  )

non_tumor_cols <- names(df)[
  which(names(df)=="Oligodendrocyte") :
    which(names(df)=="Radial_Glia")
]

#dirichlet regression####
#change refernce group here. #tried uisng arterial vessel as reference group--doesn't work. Only works for pericyte & astrocyte
df_subset <- filter(df, subtype == "AC_like_Cancer")
df_subset$Y <- DR_data(as.matrix(df_subset[ non_tumor_cols ]))
res_AC <- DirichReg(Y ~ REACTOME_STING + study_f, data = df_subset, model="alternative", base = 1)

df_subset <- filter(df, subtype == "MES_like_Cancer")
df_subset$Y <- DR_data(as.matrix(df_subset[ non_tumor_cols ]))
res_MES <- DirichReg(Y ~ REACTOME_STING + study_f, data = df_subset, model="alternative", base = 1)

df_subset <- filter(df, subtype == "OPC_like_Cancer")
df_subset$Y <- DR_data(as.matrix(df_subset[ non_tumor_cols ]))
res_OPC <- DirichReg(Y ~ REACTOME_STING + study_f, data = df_subset, model="alternative", base = 1)

df_subset <- filter(df, subtype == "NPC_like_Cancer")
df_subset$Y <- DR_data(as.matrix(df_subset[ non_tumor_cols ]))
res_NPC <- DirichReg(Y ~ REACTOME_STING + study_f, data = df_subset, model="alternative", base = 1)

res_list <- list(
  AC  = res_AC,
  MES = res_MES,
  OPC = res_OPC,
  NPC = res_NPC
)

combined_df <- imap_dfr(res_list, function(fit, subtype) {
  # pull out the raw coef & SE vectors
  coefs <- fit$coefficients
  ses   <- fit$se
  
  # find just the STING_s rows
  idx   <- grep("STING", names(coefs), fixed = TRUE)
  
  # parse out the cell name from the parameter name
  #   e.g. if names(coefs)[i] == "STING_s.Oligodendrocyte" 
  #         cell will become "Oligodendrocyte"
  cell  <- sub(":REACTOME_STING", "", names(coefs)[idx])
  
  estimate  <- coefs[idx]
  std.error <- ses[idx]
  
  # compute 95% CIs and p-values
  conf.low  <- estimate - 1.96 * std.error
  conf.high <- estimate + 1.96 * std.error
  p.value   <- 2 * pnorm(-abs(estimate / std.error))
  
  tibble(
    subtype    = subtype,
    cell       = cell,
    estimate   = estimate,
    conf.low   = conf.low,
    conf.high  = conf.high,
    p.value    = p.value,
    adj.p.value = p.adjust(p.value, method = "BH")
  )
})

#plot the regressions as a vertical caterpillar plot####
#lineage map
lineages <- c(
  CD8_Effector_Memory_T_Cells  = "Lymphoid",
  Pro_inflammatory_Microglia   = "Myeloid",
  Radial_Glia                  = "Normal Tissue",
  Astrocyte                    = "Normal Tissue",
  OPC                          = "Normal Tissue",
  Neuron                       = "Normal Tissue",
  Oligodendrocyte              = "Normal Tissue",
  Arterial_Vessel              = "Endothelial",
  Capilary_Vessel              = "Endothelial",
  Tiplike_Vessel               = "Endothelial",
  Pericyte                     = "Endothelial",
  Hypoxia_associated_Monocytes = "Myeloid",
  Naive_Monocyte               = "Myeloid",
  BDM_Hypoxia_and_MES_like     = "Myeloid",
  Anti_inflammatory_Monocyte   = "Myeloid",
  Anti_inflammatory_BDM        = "Myeloid",
  Microglia_Aging_Signature    = "Myeloid",
  BDM_IFN_Signature            = "Myeloid",
  Proliferative_Microglia      = "Myeloid",
  BDM_MHC_Expression           = "Myeloid",
  Conventional_DC_1            = "Myeloid",
  Conventional_DC_2            = "Myeloid",
  Dendritic_like_Cell          = "Myeloid",
  Plasmacytoid_DC              = "Myeloid",
  Mast_Cell                    = "Myeloid",
  T_Cells_Stress_Signature     = "Lymphoid",
  Proliferating_T_Cells        = "Lymphoid",
  NK_like_CD8_T_Cells          = "Lymphoid",
  Cytotoxic_CD8_T_Cell         = "Lymphoid",
  Plasma_Cell                  = "Lymphoid",
  B_Cell                       = "Lymphoid",
  Resting_CD4_T_Cells          = "Lymphoid",
  CD4_T_Cell_IFN_Signature     = "Lymphoid",
  Regulatory_T_Cells           = "Lymphoid",
  Natural_Killer_Cells         = "Lymphoid"
)

lin_df <- enframe(lineages, name = "cell", value = "lineage")

df_plot <- combined_df %>%
  # transform log-odds → fold-change ratios
  mutate(
    ratio      = exp(estimate),
    ratio_low  = exp(conf.low),
    ratio_high = exp(conf.high)
  ) %>%
  # join lineage info
  left_join(lin_df, by = "cell") %>%
  # within each subtype, order cells by descending ratio
  group_by(subtype) %>%
  arrange(desc(ratio), .by_group = TRUE) %>%
  mutate(
    cell        = factor(cell, levels = unique(cell)),
    # create text labels
    label       = sprintf("%.2f [%.2f, %.2f]", ratio, ratio_low, ratio_high),
    # flag significance where CI excludes 1
    sig         = (ratio_low > 1) | (ratio_high < 1),
    label_color = ifelse(sig, lineage, "ns"),
    # fixed x-position per facet
    label_x     = max(ratio_high) + 0.05
  ) %>%
  ungroup()

cols <- c(
  Lymphoid        = "#bc3c29",
  Myeloid         = "#e18727",
  Endothelial     = "#0072b5",
  `Normal Tissue` = "#20854e"
)

subtypes <- unique(df_plot$subtype)

# Loop over each subtype and build + print a separate plot
for (s in subtypes) {
  df_sub <- df_plot %>%
    filter(subtype == s) %>%
    arrange(desc(ratio)) %>%
    mutate(cell = factor(cell, levels = cell))  # lock in that order

  p <- ggplot(df_sub, aes(x = ratio, y = cell, color = lineage)) +
    geom_vline(xintercept = 1, linetype = "solid", color = "black") +
    geom_point(size = 2) +
    geom_errorbarh(aes(xmin = ratio_low, xmax = ratio_high), height = 0.2) +
    geom_text(
      aes(x = label_x, label = label, color = label_color),
      hjust = 0,
      size  = 2.5
    ) +
    scale_color_manual(values = cols) +
    scale_x_continuous(
      expand = expansion(mult = c(0.01, 0.15)),
      limits = c(min(df_sub$ratio_low), NA)
    ) +
    labs(
      x     = "μᵢ / μref (fold-change)",
      y     = NULL,
      color = "Lineage",
      title = paste("Subtype:", s)
    ) +
    theme_bw() +
    theme(
      axis.text.y     = element_text(size = 12),
      legend.position = "bottom",
      plot.title      = element_text(hjust = 0.5),
      plot.margin     = margin(5, 5, 5, 5)
    ) +
    coord_cartesian(clip = "off")
  
  # 3. Print the plot
  print(p)
}

