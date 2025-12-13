#libraries####
library(data.table)
library(dplyr)
library(survival)
library(ggplot2)
library(forcats)
library(compositions)
library(tidyr)
library(purrr)
library(broom)         # for tidy()
library(metafor)       # for rma()
library(DirichletReg)

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
      `Oligodendrocyte`:`Natural_Killer_Cells`,
      ~ . / (1 - tumor_purity)
    )
  ) %>%
  
  # 3) drop the original tumor‐cell columns
  dplyr::select(-all_of(tumor_cells))

# sanity check: the non‐tumor fractions should now sum to 1
data2 %>%
  mutate(check = rowSums(across(`Oligodendrocyte`:`Natural_Killer_Cells`))) %>%
  dplyr::select(sample_id, tumor_purity, check) %>%
  head(10)

data2 <- as.data.frame(data2)

df <- data2 %>%
  filter(recurrence == "primary", study != "TEMPUS")

cell_cols <- names(df)[ which(names(df) == "Oligodendrocyte"):
                          which(names(df) == "Natural_Killer_Cells") ]

#meta analysis####
clr_mat <- compositions::clr( as.matrix(df[cell_cols]) )
colnames(clr_mat) <- paste0("clr_", cell_cols)

df2 <- bind_cols(df, as.data.frame(clr_mat))

# 2) For each cell type, get per‐study β and SE for STING effect
meta_input <- map_df(cell_cols, function(cell){
  clr_cell <- paste0("clr_", cell)
  
  df2 %>%
    group_by(study) %>%
    # only keep studies with ≥ 3 samples to fit a model
    filter(n() >= 3) %>%
    do({
      fit <- lm(
        formula(sprintf("%s ~ STING1_tpm", clr_cell)),
        data = .
      )
      td  <- broom::tidy(fit)
      # extract the STING1_tpm term
      st <- td %>% filter(term=="STING1_tpm")
      tibble(
        cell      = cell,
        study     = unique(.$study),
        beta      = st$estimate,
        se        = st$std.error
      )
    }) %>%
    ungroup()
})


# 3) Meta‐analysis per cell
meta_results <- meta_input %>%
  group_by(cell) %>%
  nest() %>%
  mutate(
    meta_obj = map(data, ~ rma(yi = beta, sei = se, data = .x, method="REML"))
  )

# 4) Summarize and forest‐plot
meta_summaries <- meta_results %>%
  mutate(
    forest_plot = pmap(
      list(obj    = meta_obj,
           df_cell = data,
           cl     = cell),
      function(obj, df_cell, cl) {
        n_studies <- nrow(df_cell)
        forest(obj,
               slab = df_cell$study,
               xlab = paste("STING →", cl),
               main = paste0(cl, " (N=", n_studies, ")"))
      }
    )
  )

# And to generate forest plots for each:
walk(meta_summaries$forest_plot, invisible)
