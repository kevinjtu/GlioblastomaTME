# libraries ####
library(compositions)   
library(dplyr)
library(tidyr)
library(pheatmap)
library(ggplot2)
library(stats)
library(purrr)
library(DirichletReg)

# 1. Import data & DEFINE YOUR COLUMN SETS ---------------------------------------------
setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME")
data <- as.data.frame(fread("GBM Clinical-TME Cell States.csv"))

tumor.cols <- c("OPC_like_Cancer","AC_like_Cancer","NPC_like_Cancer","MES_like_Cancer")

# select all other columns that are TME (here by negative selection; adjust as needed)
tme.cols   <- setdiff(colnames(data), c(tumor.cols, 
                                        # non-cell metadata
                                        "sample_id","patient_id","age","sex","race","grade",
                                        "OS_event","OS_years","DFS_event","DFS_years",
                                        "PFS_event","PFS_years","recurrence","location",
                                        "IDH_status","study","ARID1A_tpm","STING1_tpm"
))

# 2) prepare a tiny pseudocount + renormalize so each row sums to 1 ####
pseudo    <- 1e-6
tme.mat   <- data[tme.cols] + pseudo
tme.norm  <- tme.mat / rowSums(tme.mat)

# 3) turn it into a DirichletRegData object ####
Y <- DirichletReg::DR_data(tme.norm)

# 4) assemble a modeling data.frame ####
mod.df      <- data
mod.df[tme.cols] <- tme.norm
mod.df$Y    <- Y

# 5) fit Dirichlet regression: ####
#    response Y ~ the four tumor fractions + any covariates you like
fit.dirich <- DirichletReg::DirichReg(
  Y ~ OPC_like_Cancer + AC_like_Cancer + NPC_like_Cancer + MES_like_Cancer 
  + age + sex + IDH_status,     # ← add/drop covariates as desired
  data = mod.df
)


summary(fit.dirich)

# 6) organize results for visualization####
beta    <- fit.dirich$coefficients       # named vector of all coefficients
se      <- fit.dirich$se                 # matching named vector of standard errors
coef_nm <- names(beta)                   # names like "muMES.like.Cancer_AC.like.Cancer" or "phi_depth"

coefs_df <- data.frame(
  term     = coef_nm,
  estimate = as.numeric(beta),
  se       = as.numeric(se),
  stringsAsFactors = FALSE
)

terms  <- coefs_df$term

# 1) matrix: mu vs phi
coefs_df$matrix <- ifelse(
  startsWith(terms, "mu"), 
  "mu",
  ifelse(startsWith(terms, "phi"), "phi", NA_character_)
)

# 2) body = everything after "mu" or "phi"
body <- sub("^(mu|phi)", "", terms)

# 3) component = text before first "_"
coefs_df$component <- sub("^([^_]+)_.*$", "\\1", body)

# 4) predictor = text after first "_", or "(Intercept)" if none
has_und <- grepl("_", body)
predictor <- rep("(Intercept)", length(body))
predictor[has_und] <- sub("^.*_(.+)$", "\\1", body[has_und])
coefs_df$predictor <- predictor

# 5) cleanup: drop `body`, reorder if you like
coefs_df$term  <- NULL

coefs_df <- coefs_df %>%
  mutate(
    z       = estimate / se,
    p.value = 2 * pnorm(-abs(z))
  )

# quick check
head(coefs_df, 10)

# 7) Compute CIs and significance ####
coefs_plot <- coefs_df %>%
  # focus on the "mean" portion of the Dirichlet model
  #filter(matrix == "mu") %>%
  mutate(
    lower = estimate - 1.96 * se,
    upper = estimate + 1.96 * se,
    sig   = .data$`p.value` < 0.05
  )