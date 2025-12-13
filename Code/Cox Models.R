#libraries####
library(data.table)
library(dplyr)
library(survival)
library(ggplot2)
library(forcats)
library(compositions)

#Set  survival metric here: "OS", "PFS", or "DFS"
Surv_metric <- "OS"

#import & clean data####
setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME")
data <- as.data.frame(fread("GBM Clinical-TME Cell States.csv"))

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
  filter(recurrence == "primary") %>%
  filter(is.na(location))

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

# 4) rename columns to mark them as CLR‐transformed
colnames(clr_mat) <- paste0(cell_cols, "_clr")

# 5) rebuild your analysis data frame, replacing raw with CLR
df_clr <- df %>%
  dplyr::select(-all_of(cell_cols)) %>%
  bind_cols(as.data.frame(clr_mat))

# 6) update the list of columns you’ll loop over
cell_cols_clr <- paste0(cell_cols, "_clr")

#Cox Model####
time_col  <- paste0(Surv_metric, "_years")
event_col <- paste0(Surv_metric, "_event")

res_list <- list()

df_clr_idh <- df_clr #%>%
  #filter(IDH_status == "mut")

for(cell in cell_cols_clr) {
  df_cell <- df_clr_idh %>%
    mutate(Q = ntile(.data[[cell]], 4)) %>%
    filter(!is.na(Q))
  
  if(sd(df_cell$Q)==0) {
    message("Skipping ", cell, ": no quartile variation.")
    next
  }
  
  surv_fml <- as.formula(paste0(
    "Surv(", time_col, ", ", event_col, ") ~ Q + KarnPerfScore + age +",
    "strata(study)"
  ))
  fit  <- coxph(surv_fml, data = df_cell)
  sfit <- summary(fit)
  
  if(!"Q" %in% rownames(sfit$coefficients)) {
    message("Skipping ", cell, ": no 'Q' term in model.")
    next
  }
  
  hr    <- sfit$coefficients["Q","exp(coef)"]
  lower <- sfit$conf.int     ["Q","lower .95"]
  upper <- sfit$conf.int     ["Q","upper .95"]
  n_obs  <- sfit$n            # total rows used in fitting
  n_evt  <- sfit$nevent       # number of events
  p_val <- sfit$coefficients["Q","Pr(>|z|)"]
  
  res_list[[cell]] <- tibble(
    cell  = cell,
    HR    = hr,
    lower = lower,
    upper = upper,
    n      = n_obs,
    events = n_evt,
    p_value = p_val
  )
}

results_df <- bind_rows(res_list)

# Forest‐style plot ####
n_tot      <- results_df$n[1]
events_tot <- results_df$events[1]

results_df$cell <- gsub("_clr$", "", results_df$cell)

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

# Merge lineage into results_df
results_df_lineage <- results_df %>% 
  merge(lin_df, by = "cell")

# Define colors
cols <- c(
  Lymphoid        = "#20854e",
  Myeloid         = "#bc3c29",
  Endothelial     = "#54278f",
  `Normal Tissue` = "#0072b5"
)

# Create formatted label column
results_df_lineage <- results_df_lineage %>%
  mutate(
    HR_label = sprintf("%.2f (%.2f - %.2f)", HR, lower, upper)
  )

#plot
ggplot(results_df_lineage, aes(x = HR, y = fct_reorder(cell, HR), color = lineage)) +
  geom_vline(xintercept = 1, linetype = "solid") + 
  geom_point(size = 2, shape = 15) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2) +
  
  # Add formatted HR (CI) labels OUTSIDE the plot area
  geom_text(
    aes(x = max(upper) * 1, label = HR_label),
    hjust = 0,
    size = 3
  ) +
  
  annotate(
    "text", 
    x = 1.5, y = -Inf, 
    label = paste0("n=", n_tot, ", events=", events_tot),
    hjust = 1.1, vjust = -0.5, 
    size = 3
  )+
  scale_color_manual(values = cols) +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.2))) + 
  labs(
    x     = "Hazard Ratio",
    y     = NULL,
    title = paste("Cox Model,", Surv_metric),
    color = "Lineage"
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black"),
    #plot.margin = margin(5.5, 80, 5.5, 5.5)  # add extra right margin for text
    legend.position = "bottom"
  ) +
  coord_cartesian(clip = "off")

