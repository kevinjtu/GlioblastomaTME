#libraries####
library(InstaPrism)
library(data.table)
library(Biobase)
library(dplyr)
library(textshape)
library(tibble)
library(readxl)
library(openxlsx)
library(AnnotationDbi)
library(limma)
library(msigdbr)
library(clusterProfiler)
library(GSVA)
library(stringr)
library(survival)
library(survminer)
library(ordinal)
library(brms)
library(tidyr)

#STING niche Enrichment####
# Define niche gene sets
niche_gsea <- list(
  niche_all = c(  # For overall enrichment
    "HILPDA", "LOX", "VEGFA", "CAV1", "CXCL8", "VIM", "TNC",
    "BCAN", "OLIG1", "KDR", "FLT1", "TMEM173", "PECAM1"
  ),
  niche_up = c(   # For directionality
    "HILPDA", "LOX", "VEGFA", "CAV1", "CXCL8", "VIM", "TNC"
  ),
  niche_down = c( # For directionality
    "BCAN", "OLIG1", "KDR", "FLT1", "TMEM173", "PECAM1"
  )
)

# load the dataset, neo1
setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Expression Data")
neoadjuvant1 <- read_excel("Nassiri Nanostring.xlsx") %>%
  column_to_rownames(var = "Gene")

#get into correct data type
bulk_ma <- apply(neoadjuvant1, 2, as.numeric)
rownames(bulk_ma) <- rownames(neoadjuvant1)
bulk_ma <- as.matrix(bulk_ma)
if ("STING1" %in% rownames(bulk_ma)) {
  rownames(bulk_ma)[rownames(bulk_ma) == "STING1"] <- "TMEM173"}

params <- ssgseaParam(
  exprData = bulk_ma,
  geneSets = niche_gsea,
  normalize = FALSE
)

ssgsea_eset_neo1 <- as.data.frame(t(gsva(
  params
)))
ssgsea_eset_neo1$patient_id <- rownames(ssgsea_eset_neo1)

ssgsea_eset_neo1 <- ssgsea_eset_neo1 %>%
  dplyr::rename(sample_id = patient_id) %>%                            # Rename column
  dplyr::mutate(patient_id = str_extract(sample_id, "^[^.]+"))  

#clinical####
setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME/Clinical Data")

clinical <- read_excel("Nassari Clinical.xlsx") %>%
  dplyr::rename(patient_id = "Subject ID") %>%
  merge(ssgsea_eset_neo1, by = "patient_id")

#box plots####
# Step 1: Filter, rename, and reorder response categories
clinical_filtered <- clinical %>%
  filter(!grepl("\\.2$", sample_id)) %>%
  mutate(`Best Response` = recode(`Best Response`,
                                  "CR" = "Complete Response",
                                  "PR" = "Partial Response",
                                  "SD" = "Stable Disease",
                                  "PD" = "Progressive Disease"),
         `Best Response` = factor(`Best Response`, 
                                  levels = c("Complete Response", 
                                             "Partial Response", 
                                             "Stable Disease", 
                                             "Progressive Disease")))

# Step 2: Define pairwise comparisons
comparisons <- combn(levels(clinical_filtered$`Best Response`), 2, simplify = FALSE)

# Step 3: Plot
ggplot(clinical_filtered, aes(x = `Best Response`, y = niche_all)) +
  geom_boxplot(fill = "skyblue", alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.6, color = "black") +
  stat_compare_means(comparisons = comparisons, method = "wilcox.test", label = "p.signif") +
  labs(
    title = "Niche Score by Best Response",
    x = "Best Response",
    y = "niche_all"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid = element_blank()
  )

#ordinal regressions####
# Step 1: Clean and prepare the data
clinical_model <- clinical %>%
  filter(!grepl("\\.2$", sample_id)) %>%
  rename(dose = "DNX- 2401 dose (vp)") %>%
  rename(PDL1 = "PD-L11") %>%
  mutate(
    # Recode response as ordered factor
    Best_Response = recode(`Best Response`,
                           "CR" = "Complete Response",
                           "PR" = "Partial Response",
                           "SD" = "Stable Disease",
                           "PD" = "Progressive Disease"),
    
    Best_Response = factor(Best_Response,
                           levels = c("Complete Response", "Partial Response",
                                      "Stable Disease", "Progressive Disease"),
                           ordered = TRUE),
    
    # Binary Responder status: everything but PD is a responder
    Responder_Status = ifelse(Best_Response == "Progressive Disease",
                              "Non-Responder", "Responder"),
    Responder_Status = factor(Responder_Status, levels = c("Responder", "Non-Responder")),
    
    # Quartile niche_all
    niche_all_quartile = ntile(niche_all, 3)
  )

# Step 2: Fit an binary logistic regression
model_glm <- glm(Responder_Status ~ niche_all_quartile,
                 data = clinical_model,
                 family = binomial(link = "logit"))
summary(model_glm)

ggplot(clinical_model, aes(x = factor(niche_all_quartile), fill = Responder_Status)) +
  geom_bar(position = "fill") +
  ylab("Proportion") +
  xlab("niche_all Quartile") +
  ggtitle("Distribution of Best Response Across niche_all Quartiles") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

model_brm <- brm(
  Responder_Status ~ niche_all_quartile,
  data = clinical_model,
  family = bernoulli(link = "logit"),
  prior = set_prior("normal(0, 1)", class = "b")  # shrink large estimates
)

summary(model_brm)
posterior <- as_draws_df(model_brm)
mean(posterior$b_niche_all_quartile > 0)
#before/after treatment####
# Filter samples with .1 or .2 suffix
paired_samples <- clinical[grepl("\\.1$|\\.2$", clinical$sample_id), ]

# Extract base IDs (before .1 or .2)
paired_samples$base_id <- sub("\\.\\d$", "", paired_samples$sample_id)

# Keep only base IDs that have both .1 and .2
valid_ids <- names(table(paired_samples$base_id))[table(paired_samples$base_id) == 2]
paired_samples <- paired_samples[paired_samples$base_id %in% valid_ids, ]

paired_wide <- pivot_wider(paired_samples[, c("base_id", "sample_id", "niche_all")],
                           names_from = sample_id, values_from = niche_all)

# Identify .1 and .2 columns
sample_cols <- names(paired_wide)[-1]  # exclude 'base_id'
sample_ids_1 <- grep("\\.1$", sample_cols, value = TRUE)
sample_ids_2 <- sub("\\.1$", ".2", sample_ids_1)

# Keep only pairs where both .1 and .2 exist
valid_pairs <- sample_ids_1[sample_ids_2 %in% names(paired_wide)]
sample_ids_2 <- sub("\\.1$", ".2", valid_pairs)

# Extract matched values from each pair
vals_1 <- as.numeric(unlist(paired_wide[valid_pairs]))
vals_2 <- as.numeric(unlist(paired_wide[sample_ids_2]))

# Remove NA pairs
valid_idx <- complete.cases(vals_1, vals_2)
vals_1 <- vals_1[valid_idx]
vals_2 <- vals_2[valid_idx]

# Perform paired t-test
t_result <- t.test(vals_1, vals_2, paired = TRUE)
print(t_result)

long_df <- data.frame(
  value = c(vals_1, vals_2, (vals_1 + vals_2)/2),
  group = rep(c(".1", ".2", "Mean"), each = length(vals_1))
)

# Plot
ggplot(long_df, aes(x = group, y = value, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, alpha = 0.4, size = 1) +
  labs(title = "Paired Sample Values and Their Means",
       y = "niche_all") +
  scale_fill_manual(values = c(".1" = "#66c2a5", ".2" = "#fc8d62", "Mean" = "#8da0cb")) +
  theme_minimal(base_size = 14)