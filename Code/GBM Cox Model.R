#libraries####
library(data.table)
library(Biobase)
library(dplyr)
library(textshape)
library(tibble)
library(readxl)
library(survival)
library(gridExtra)
library(ggplot2)
library(rms)
library(forestmodel)
library(stringr)

surv_metric <- "OS" # options are OS, DFS

# Load the cell type data####
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

combined_standard <- combined_standard %>% mutate_at(vars(7:10), as.numeric) #set survival data to numeric
cox_results <- data.frame(Cell_Type = character(),
                          Hazard_Ratio = numeric(),
                          CI_Lower = numeric(),
                          CI_Upper = numeric(),
                          P_Value = numeric(),
                          stringsAsFactors = FALSE)

# cox regression model test at cell state with 4-quantile ####
cell_types <- unique(colnames(combined_standard[,14:22]))

for (cell in cell_types) {
    print(cell)
    
    # Create quartiles for the 'cell' variable
    combined_standard$cell_quartiles <- cut(combined_standard[[cell]], 
                                               breaks = quantile(combined_standard[[cell]], probs = seq(0, 1, 0.25), na.rm = TRUE), 
                                               labels = FALSE)
    
    # Convert 'cell_quartiles' to a numeric
    combined_standard$cell_quartiles <- as.integer(combined_standard$cell_quartiles)
    
    # Define the formula with 'cell_quartiles'
    formula <- as.formula(paste("Surv(", paste0(surv_metric, "_years"), ", ", paste0(surv_metric, "_event)"), 
                                " ~ cell_quartiles + age + sex"))
    
    # Fit the Cox proportional hazards model
    cox_model <- coxph(formula, data = combined_standard)
    summary_cox <- summary(cox_model)
    hr <- coef(cox_model)["cell_quartiles"] # Hazard Ratio
    ci <- confint(cox_model, level = 0.95)["cell_quartiles", ] # Confidence Interval
    p_value <- summary_cox$coefficients["cell_quartiles", "Pr(>|z|)"] # P-Value
    
    # Store results in the data frame
    cox_results_cell <- data.frame(Cell_Type = cell,
                                   Hazard_Ratio = hr,
                                   CI_Lower = ci[1],
                                   CI_Upper = ci[2],
                                   P_Value = p_value)
    
    cox_results <- rbind(cox_results, cox_results_cell)
  }

cox_results$P_Value_BH <- p.adjust(cox_results$P_Value, method = "BH")

#plot cox regression model results####
#plot er statuses together
cox_results <- data.table(cox_results)
cox_results <- dplyr::rename(cox_results, c("Cell_State" = "Cell_Type"))
cox_results$Cell_Type <- as.character(cox_results$Cell_Type)
cox_results$Cell_Type[grep("TAM", cox_results$Cell_State)] <- "Tumor Associated Macrophage"
cox_results$Cell_Type[grep("Endo", cox_results$Cell_State)] <- "Endothelial"
cox_results$Cell_Type[grep("Tip", cox_results$Cell_State)] <- "Endothelial"
cox_results$Cell_Type[grep("DC", cox_results$Cell_State)] <- "Dendritic Cell"
cox_results$Cell_Type[grep("Mono", cox_results$Cell_State)] <- "Monocyte"
cox_results$Cell_Type[grep("Prolif", cox_results$Cell_State)] <- "T Cell"
cox_results$Cell_Type[grep("Reg", cox_results$Cell_State)] <- "T Cell"
cox_results$Cell_Type[grep("CD", cox_results$Cell_State)] <- "T Cell"
cox_results$Cell_Type[grep("Stress", cox_results$Cell_State)] <- "T Cell"
cox_results$Cell_Type[grep("B_cell", cox_results$Cell_State)] <- "B Cell"
cox_results$Cell_Type[grep("Plasma", cox_results$Cell_State)] <- "B Cell"
cox_results$Cell_Type[grep("NK", cox_results$Cell_State)] <- "Granulocyte"
cox_results$Cell_Type[grep("Mast", cox_results$Cell_State)] <- "Granulocyte"
cox_results$Cell_Type[grep("OPC", cox_results$Cell_State)] <- "Oligodendrocyte Precursor Cell"
cox_results$Cell_Type[grep("RG", cox_results$Cell_State)] <- "Radial Glia"
cox_results$Cell_Type[grep("AC", cox_results$Cell_State)] <- "Malignant"
cox_results$Cell_Type[grep("MES", cox_results$Cell_State)] <- "Malignant"
cox_results$Cell_Type[grep("OPClike", cox_results$Cell_State)] <- "Malignant"
cox_results$Cell_Type[grep("NPC", cox_results$Cell_State)] <- "Malignant"

#Plot HR
cox_results$Cell_State <- factor(cox_results$Cell_State, levels = cox_results$Cell_State[order(cox_results$Hazard_Ratio, decreasing = TRUE)])
ggplot(cox_results, aes(x = Hazard_Ratio, y = Cell_State, fill = Cell_Type)) +
  geom_point(size = 2, shape = 22, position = position_dodge(width = 0.5)) +
  geom_errorbarh(aes(xmin = CI_Lower, xmax = CI_Upper), height = 0, position = position_dodge(width = 0.5)) +
  labs(
    title = paste("Forest Plot for", surv_metric),
    x = "Hazard Ratio",
    y = "Cell State",
    fill = "Cell Type"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 9),
    plot.title = element_text(hjust = 0.5)
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  guides(fill = guide_legend(override.aes = list(size = 5)))


