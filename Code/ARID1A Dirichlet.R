#libraries####
library(data.table)
library(dplyr)
library(DirichletReg)
library(purrr)
library(tibble)
library(ggplot2)

#Import data & dirichlet regression####
setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME")
data <- as.data.frame(fread("GBM Clinical-TME Cell States.csv"))

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
    ARID1A_s = as.numeric(scale(ARID1A_tpm))
  )

non_tumor_cols <- names(df)[
  which(names(df)=="Oligodendrocyte") :
    which(names(df)=="Radial_Glia")
]

#change refernce group here. #tried uisng arterial vessel as reference group--doesn't work. Only works for pericyte & astrocyte
my_cells <- setdiff(non_tumor_cols, "Arterial_Vessel") 
my_cells <- c(my_cells, "Arterial_Vessel")

df$DR_nt <- DR_data(
  as.matrix(df[ my_cells ]),
  trafo = 0
)

fit_nt <- DirichReg(
  DR_nt ~ ARID1A_s + study_f,
  data    = df,
  control = list(maxit = 1e4, reltol = 1e-8)
)

summary(fit_nt)

#visualize with vertical caterpillar plot####
s <- summary(fit_nt)
coef_mat <- as.data.frame(s$coef.mat)

n_terms <- nrow(coef_mat) / length(s$varnames)  # should be 9

coef_df <- coef_mat %>%
  rownames_to_column("variable") %>%
  mutate(
    cell     = rep(s$varnames, each = n_terms),
    CI_lower = Estimate - 1.96 * `Std. Error`,
    CI_upper = Estimate + 1.96 * `Std. Error`
  ) %>%
  dplyr::select(cell, variable, Estimate, `Std. Error`, `z value`, `Pr(>|z|)`, CI_lower, CI_upper)

df_arid <- coef_df %>%
  filter(grepl("ARID1A_s", variable)) %>%
  transmute(
    cell     = cell,
    estimate = exp(Estimate),
    se       = `Std. Error`,
    CI_lower = exp(Estimate - 1.96 * `Std. Error`),
    CI_upper = exp(Estimate + 1.96 * `Std. Error`)
  )





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

df_plot <- df_arid %>%
  left_join(lin_df, by = "cell") %>%
  arrange(desc(estimate)) %>%
  mutate(
    cell  = factor(cell, levels = cell),
    label = sprintf("%.3f [%.3f, %.3f]", estimate, CI_lower, CI_upper)
  )

cols <- c(
  Lymphoid        = "#bc3c29",
  Myeloid         = "#e18727",
  Endothelial     = "#0072b5",
  `Normal Tissue` = "#20854e"
)

# 1. Flag significance (95% CI that excludes 0)
df_plot <- df_plot %>%
  mutate(
    sig = ifelse(CI_lower > 0 | CI_upper < 0, TRUE, FALSE),
    label_color = ifelse(sig, lineage, "ns")
  )

# 2. Add "ns" to color scale for non-significant labels
cols_full <- c(cols, ns = "grey50")

# 3. Define fixed label x-position
x_pos <- max(df_plot$CI_upper) + 0.05

x_lim_max <- max(df_plot$CI_upper)
df_plot <- df_plot %>%
  mutate(
    label_x = x_lim_max + 0.05  # position outside plot area
  )

ggplot(df_plot, aes(x = estimate, y = cell, color = lineage)) +
  geom_vline(xintercept = 1, linetype = "solid", color = "black") +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper), height = 0.2) +
  geom_text(
    aes(x = label_x, label = label, color = label_color),
    hjust = 0,
    size  = 3
  ) +
  scale_color_manual(values = cols_full) +
  scale_x_continuous(
    expand = expansion(mult = c(0.01, 0.15)),
    limits = c(min(df_plot$CI_lower), x_lim_max + 0.12)
  ) +
  labs(
    x     = "μᵢ/μ_ref",
    y     = NULL,
    color = "Lineage",
    title = "Dirichlet Regression, ARID1A Effect"
  ) +
  theme_bw() +
  theme(
    axis.text.y     = element_text(size = 8),
    plot.margin     = margin(5, 150, 5, 5),  # expand right margin
    legend.position = "bottom",
    plot.title      = element_text(hjust = 0.5)
  ) +
  coord_cartesian(clip = "off")  # allow text to overflow

