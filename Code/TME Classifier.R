#libraries####
library(caret)
library(randomForest)
library(pROC)
library(glmnet)
library(data.table)
library(dplyr)
library(compositions)
library(readxl)
library(ggplot2)
library(survival)
library(survminer)
library(ggeffects)
library(reshape2)

set.seed(42)

#import & clean data####
subtype_cols <- c("OPC_like_Cancer", "AC_like_Cancer", "NPC_like_Cancer", "MES_like_Cancer")
setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME")
data <- as.data.frame(fread("GBM Clinical-TME Cell States.csv")) %>%
  dplyr::mutate(subtype = subtype_cols[
    max.col( as.matrix( across(all_of(subtype_cols)) ),
             ties.method = "first")])

#get TME fractions accounting for cellularity
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
  filter(recurrence == "primary")

cell_cols <- names(df)[ which(names(df) == "Oligodendrocyte"):
                          which(names(df) == "Natural_Killer_Cells") ]

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


#annotate df for classifier
setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Heatmap - TME Signatures")
clusters <- read_excel("GBM TME Cluster Assignments.xlsx") %>%
  dplyr::select(patient_id, TME_signature)

df_class <- df_clr %>%
  filter(study == "TCGA", IDH_status == "wt", recurrence == "primary") %>%
  dplyr::select(patient_id, all_of(cell_cols_clr)) %>%
  merge(clusters, by = "patient_id") %>%
  dplyr::select(-patient_id)

# make sure your outcome is a factor####
df_class$TME_signature <- factor(df_class$TME_signature, levels=c("Immunosuppressed","Immunoactive"))

# split into train/test
train_idx <- createDataPartition(df_class$TME_signature, p = 0.7, list = FALSE)
train_df  <- df_class[train_idx, ]
test_df   <- df_class[-train_idx, ]

#train the classifier
ctrl <- trainControl(
  method           = "cv",
  number           = 5,
  classProbs       = TRUE,
  summaryFunction  = twoClassSummary,    # enables ROC
  savePredictions  = TRUE
)

rf_mod <- train(
  TME_signature ~ .,
  data    = train_df,
  method  = "rf",
  metric  = "ROC",
  trControl = ctrl
)

print(rf_mod)

# get predicted probabilities and classes
probs <- predict(rf_mod, test_df, type = "prob")[, "Immunoactive"]
preds <- predict(rf_mod, test_df)

# confusion matrix
confusionMatrix(preds, test_df$TME_signature)

# ROC & AUC
roc_obj <- roc(test_df$TME_signature, probs, levels=c("Immunosuppressed","Immunoactive"))
auc(roc_obj)
plot(roc_obj)
varImp(rf_mod) %>% plot(top = 20)

#apply the classifier to the metadata-set####
classified_df <- df_clr %>%
  filter(
    !study %in% c("IVYGAP"),
    recurrence  == "primary") %>%  
  dplyr::select(sample_id, all_of(cell_cols_clr)) %>%
  column_to_rownames("sample_id")

new_probs <- predict(rf_mod, classified_df, type="prob")[, "Immunoactive"]
new_call  <- predict(rf_mod, classified_df)

classified_df$TME_signature <- new_call
classified_df <- classified_df %>%
  rownames_to_column("sample_id") %>%
  dplyr::select(sample_id, TME_signature)

#write csv of the data####
setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Heatmap - TME Signatures")
write.csv(classified_df, "GBM TME Classifier Results.csv", row.names = TRUE)

#subtype x ecotype table chi square####
survival <- df_clr %>%
  merge(classified_df, by = "sample_id") %>%
  filter(platform == "rnaseq")

observed <- table(survival$subtype, survival$TME_signature)
observed_t <- t(observed)
chisq_result <- chisq.test(observed_t)
print(chisq_result)
residuals <- chisq_result$stdres

# Reshape for ggplot
df <- melt(residuals)
colnames(df) <- c("TME_signature", "Subtype", "Residual")  # flipped order due to transpose

# Plot heatmap
ggplot(df, aes(x = Subtype, y = TME_signature, fill = Residual)) +
  geom_tile(color = "black") +  # black outline for each square
  scale_fill_gradient2(low = "#BC3C29", mid = "white", high = '#20854E', midpoint = 0,
                       limits = c(-max(abs(df$Residual)), max(abs(df$Residual))),
                       name = "Std. Residual") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),      # remove gridlines
    panel.border = element_blank(),    # remove plot border
    axis.ticks = element_blank(),      # remove axis ticks
    axis.text.x = element_text(angle = 45, hjust = 1)  # optional: angle x labels
  ) +
  labs(title = "Chi-square Standardized Residuals",
       x = "",
       y = "")

#survival cox model####
survival <- df_clr %>%
  merge(classified_df, by = "sample_id") %>%
  filter(platform == "rnaseq")
  #filter(study %in% c("TCGA", "CGGA Batch 2", "GLASS"))

survival_clean <- survival %>%
  filter(!is.na(OS_event), !is.na(OS_years)) %>%
  mutate(
    subtype = factor(subtype),
    TME_signature = factor(TME_signature),
    study = factor(study)
  )

cox_data <- subset(survival, !is.na(OS_years) & !is.na(OS_event))

cox_data$subtype <- factor(cox_data$subtype, ordered = FALSE)
cox_data$subtype <- relevel(cox_data$subtype, ref = "AC_like_Cancer")

cox_model <- coxph(Surv(OS_years, OS_event) ~ subtype * TME_signature + age + strata(study), data = cox_data)
summary(cox_model)

# Predict hazard ratios for the interaction term: subtype * TME_signature
interaction_plot_data <- ggpredict(cox_model, terms = c("TME_signature", "subtype"))

# Line plot with confidence intervals
ggplot(interaction_plot_data, aes(x = x, y = predicted, group = group, color = group)) +
  geom_line(size = 1.2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.2, color = NA) +
  labs(
    x = "TME Signature",
    y = "Predicted Hazard Ratio",
    color = "Subtype",
    fill = "Subtype",
    title = "Interaction Effects of Subtype and TME Signature on Survival"
  ) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )