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

#Combile data####
setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME")
combined_standard <- read_excel("GBM TME-Clinical Data.xlsx")
print(names(combined_standard))
new_names <- names(combined_standard)
new_names <- str_replace_all(new_names, "\\-", "")       # Remove hyphens
new_names <- str_replace_all(new_names, "\\/", "_")       # Remove slashes
new_names <- str_replace_all(new_names, " ", "_")         # Replace spaces with underscores
new_names <- str_replace_all(new_names, "\\+", "")         # Remove plus signs
combined_standard <- setNames(combined_standard, new_names)
print(names(combined_standard))

opc_like_columns <- grep("OPClike", names(combined_standard), value = TRUE)
ac_like_columns <- grep("AClike", names(combined_standard), value = TRUE)
npc_like_columns <- grep("NPClike", names(combined_standard), value = TRUE)
mes_like_columns <- grep("MESlike", names(combined_standard), value = TRUE)

combined_standard$OPC_like <- rowSums(combined_standard[opc_like_columns], na.rm = TRUE)
combined_standard$AC_like <- rowSums(combined_standard[ac_like_columns], na.rm = TRUE)
combined_standard$NPC_like <- rowSums(combined_standard[npc_like_columns], na.rm = TRUE)
combined_standard$MES_like <- rowSums(combined_standard[mes_like_columns], na.rm = TRUE)
combined_standard <- combined_standard[, !(names(combined_standard) %in% c(opc_like_columns, ac_like_columns, npc_like_columns, mes_like_columns))]

setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME")
sting <- fread("STING1__mRNA_expression_(RNA_Seq_V2_RSEM).txt")
sting <- na.omit(sting)

mclust <- Mclust(sting$`STING1: mRNA expression (RNA Seq V2 RSEM)`)
sting_cat <- as.data.frame(mclust$classification)
sting_cat <- rename(sting_cat, STING_group = "mclust$classification")
sting_cat <- cbind(sting, sting_cat)
sting_cat <- rename(sting_cat, sample_id = "Sample ID")

combined_standard <- merge(combined_standard, sting_cat, by.x = "sample_id", all.x = TRUE)
#Cohens's D and T test to find signifiacant differences####
sting_comp <- as.data.frame(combined_standard[,c(1, 14:55, 59)])

cohens_d_results <- list()
for (col_name in colnames(sting_comp)[2:43]) {
  group1 <- sting_comp[sting_comp$STING_group == 1, col_name]
  group2 <- sting_comp[sting_comp$STING_group == 2, col_name]
    cohen_d <- cohen.d(group1, group2, na.rm = TRUE)
    cohens_d_results[[col_name]] <- cohen_d$estimate
}

cohens_d_df <- data.frame(
  Cell_Type = names(cohens_d_results),
  Cohens_d = unlist(cohens_d_results)
)

# Initialize an empty data frame to store t test results
t_test_results <- data.frame(
  Cell_Type = character(),
  t_statistic = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

for (col_name in colnames(sting_comp)[2:43]) {
  group1 <- sting_comp[sting_comp$STING_group == 1, col_name]
  group2 <- sting_comp[sting_comp$STING_group == 2, col_name]
    t_test <- t.test(group1, group2, na.rm = TRUE)
    t_test_results <- rbind(
    t_test_results,
    data.frame(
      Cell_Type = col_name,
      t_statistic = t_test$statistic,
      p_value = t_test$p.value
    )
  )
}
t_test_results$adjusted_p_value <- p.adjust(t_test_results$p_value, method = "BH")

effect_df <- merge(cohens_d_df, t_test_results, by = "Cell_Type")

#store the Cell_Type that have cohen's D > |0.8| and p-value < 0.05
significant_cells <- effect_df %>%
  filter(abs(Cohens_d) > 0.8 & p_value < 0.05)
sig_cells <- significant_cells$Cell_Type

#plotting significance cells on box plots####
long_data <- melt(sting_comp, id.vars = c("sample_id", "STING_group"), 
                  measure.vars = colnames(sting_comp)[2:43],
                  variable.name = "Cell_Type", value.name = "Proportion")

sig_cells <- c("Endo_arterial", "Pericyte", "TAMMG_aging_sig", 
               "TAMBDM_antiinfl", "Mono_antiinfl", "AC_like")
long_data_plot <- long_data %>%
  filter(Cell_Type %in% sig_cells) %>%
  mutate(STING_group = factor(ifelse(STING_group == 1, "Silenced", "Activated"), levels = c("Silenced", "Activated")))

annotation_data <- significant_cells %>%
  mutate(
    label = paste0("Cohen's d: ", round(Cohens_d, 2), "\nP: ", formatC(p_value, format = "e", digits = 2))
  ) %>%
  select(Cell_Type, label)

# Merge annotation data with the main plotting data
plot_data <- long_data_plot %>%
  inner_join(annotation_data, by = c("Cell_Type"))

# Create the plot
ggplot(plot_data, aes(x = factor(STING_group), y = Proportion, fill = factor(STING_group))) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.7) +
  geom_jitter(position = position_jitter(width = 0.2, height = 0), size = 1, alpha = 0.6) +
  facet_wrap(~ Cell_Type, scales = "free_y", ncol = 3) +
  geom_text(
    data = annotation_data,
    aes(
      x = 1.5, y = Inf, label = label  # x=1.5 centers text; y=Inf places it at the top
    ),
    vjust = 1.5, hjust = 0.5, inherit.aes = FALSE
  ) +
  theme_bw() +
  labs(
    x = "STING",
    y = "Fraction of Tumor",
    fill = "STING Group"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 10, face = "bold", hjust = 0),  # Move facet label to the left
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none"  # Move the legend to the bottom
  )