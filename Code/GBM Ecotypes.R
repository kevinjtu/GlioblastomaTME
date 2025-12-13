#libraries####
library(data.table)
library(stringr)
library(NbClust)
library(readxl)
library(dplyr)
library(ComplexHeatmap)
library(reshape2)
library(circlize)
library(compositions)
library(survminer)
library(survival)
library(ggsurvfit)
library(survival)
library(Metrics)
library(compositions)

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

opc_like_columns <- grep("OPClike", names(combined_standard), value = TRUE)
ac_like_columns <- grep("AClike", names(combined_standard), value = TRUE)
npc_like_columns <- grep("NPClike", names(combined_standard), value = TRUE)
mes_like_columns <- grep("MESlike", names(combined_standard), value = TRUE)

combined_standard$OPC_like <- rowSums(combined_standard[opc_like_columns], na.rm = TRUE)
combined_standard$AC_like <- rowSums(combined_standard[ac_like_columns], na.rm = TRUE)
combined_standard$NPC_like <- rowSums(combined_standard[npc_like_columns], na.rm = TRUE)
combined_standard$MES_like <- rowSums(combined_standard[mes_like_columns], na.rm = TRUE)
combined_standard <- combined_standard[, !(names(combined_standard) %in% c(opc_like_columns, ac_like_columns, npc_like_columns, mes_like_columns))]

#combined_standard <- combined_standard %>%
#  mutate(malignant = OPC_like + AC_like + NPC_like + MES_like) %>%
#  dplyr::select(-OPC_like, -AC_like, -NPC_like, -MES_like)%>%
#  mutate(normalizing_factor = 1 - malignant)
#combined_standard <- combined_standard %>%
#  mutate(across(Oligodendrocyte:RG, ~ .x / normalizing_factor))
#combined_standard <- combined_standard %>%
#  dplyr::select(c(-normalizing_factor, -malignant))

#compositions package, Aitchison simplex transformation (centered log ration)####
df_clr <- as.data.frame(clr(acomp(combined_standard[,14:55])))
zscore <- as.data.frame(scale(df_clr))
data <- cbind(combined_standard[,c(1:13)], zscore)

#heatmap colors####
#colors for annotations
subtype_colors <- c(
  "PN" = "#EB411E",
  "CL" = "#66A61E",
  "MS" = "#F4E83B"
  )

cell_colors <- c(
  "Normal" = "#8B4513",           # SaddleBrown
  "Dendritic_Cell" = "#1E90FF",  # DodgerBlue
  "Macrophages_Microglia" = "#FF69B4", # HotPink
  "Monocyte" = "#FF4500",        # OrangeRed
  "Innate_Immune" = "#708090",   # SlateGray
  "T_Cell" = "#DC143C",          # Crimson
  "B_Cell" = "#FFD700",          # Gold
  "Lymphoid" = "#32CD32",        # LimeGreen
  "Vascular" = "#00CED1",         # DarkTurquoise
  "Malignant" = "purple"         # prple
)

#heatmap data####
data_t <- t(data)
heatmapdata <- data_t[14:55,]
rownames <- rownames(heatmapdata)
heatmapdata <- apply(heatmapdata, 2, as.numeric)
rownames(heatmapdata) <- rownames
colnames(heatmapdata) <- data_t[1,]

#heatmap annotations####
##column annotations
res <- NbClust(data = heatmapdata, distance = "canberra", min.nc = 2, max.nc = 8, 
               method = "ward.D2", index = "ratkowsky", alphaBeale = 0.1)
res$Best.nc

hc <- hclust(dist(data[,14:55], method = "euclidean"), method = "ward.D2")
clusters <- cutree(hc, k = 3)
clusters <- as.data.frame(clusters)
clusters$sample_id <- data$sample_id
clusters <- merge(clusters, data[,1:13], by = "sample_id")
table(clusters$clusters, clusters$Subtype)

cluster_order <- order.dendrogram(as.dendrogram(hc))
clusters <- clusters[cluster_order, ]

# subtype
subtype <- HeatmapAnnotation(
  subtype = clusters$Subtype,
  col = list(subtype = subtype_colors) # Define colors for each subtype
)

##row annotations
# Creating a matching dataframe from cell states to cell types
cell_states <- rownames(heatmapdata)
cell_types <- c(
  "Normal",  # Matches "Oligodendrocyte"
  "Vascular",         # Matches "Pericyte"
  "Normal",           # Matches "Neuron"
  "Normal",        # Matches "Astrocyte"
  "Vascular",      # Matches "Endo arterial"
  "Vascular",      # Matches "Tip-like"
  "Vascular",      # Matches "Endo capilar"
  "Macrophages_Microglia",          # Matches "TAM-BDM hypoxia/MES"
  "Dendritic_Cell",               # Matches "cDC2"
  "Macrophages_Microglia",          # Matches "TAM-MG aging sig"
  "Macrophages_Microglia",          # Matches "TAM-BDM anti-infl"
  "Monocyte",          # Matches "Mono anti-infl"
  "Monocyte",          # Matches "Mono hypoxia"
  "Monocyte",          # Matches "Mono naive"
  "Macrophages_Microglia",          # Matches "TAM-MG pro-infl I"
  "Macrophages_Microglia",          # Matches "TAM-BDM MHC"
  "Macrophages_Microglia",          # Matches "TAM-MG pro-infl II"
  "Macrophages_Microglia",          # Matches "TAM-MG prolif"
  "Macrophages_Microglia",          # Matches "TAM-BDM INF"
  "Innate_Immune",          # Matches "Mast"
  "Dendritic_Cell",               # Matches "DC3"
  "Dendritic_Cell",               # Matches "DC1"
  "Dendritic_Cell",               # Matches "DC2"
  "Dendritic_Cell",               # Matches "cDC1"
  "Dendritic_Cell",               # Matches "pDC"
  "T_Cell",         # Matches "Stress sig"
  "T_Cell",         # Matches "Prolif T"
  "T_Cell",         # Matches "CD8 EM"
  "T_Cell",         # Matches "CD8 NK sig"
  "T_Cell",         # Matches "Reg T"
  "T_Cell",         # Matches "CD8 cytotoxic"
  "B_Cell",         # Matches "B cell"
  "T_Cell",         # Matches "CD4 rest"
  "T_Cell",         # Matches "CD4 INF"
  "B_Cell",         # Matches "Plasma B"
  "Innate_Immune",         # Matches "NK"
  "Normal",              # Matches "OPC"
  "Normal",      # Matches "RG"
  "Malignant",   
  "Malignant",   
  "Malignant",   
  "Malignant"
)

cell_matching <- data.frame(cell_states, cell_types)

hc_row <- hclust(dist(heatmapdata, method = "euclidean"), method = "ward.D2")
order_labels <- order.dendrogram(as.dendrogram(hc_row))
original_labels <- rownames(heatmapdata)
cell_annotation <- data.frame(order = order_labels, cell_states = original_labels[order_labels])
merged_df <- cell_annotation %>%
  merge(cell_matching, by = "cell_states")
ordered_merged_df <- merged_df[match(cell_annotation$cell_states, merged_df$cell_states),]

cell_type <- rowAnnotation(
  cell_types = ordered_merged_df$cell_types,
  col = list(cell_types = cell_colors))

#plot data into heatmap####
#fix the ordering of the metadata to the heatmapdata (make new heatmapdata2)
heatmapdata2<- heatmapdata[ordered_merged_df$cell_states,]
heatmapdata2<- heatmapdata2[,clusters$sample_id]

#heatmap plot
heatmap <- Heatmap(heatmapdata2,
                   name = "z-score",
                   clustering_distance_columns = "euclidean",
                   clustering_method_columns = "ward.D2",
                   clustering_distance_row = "euclidean",
                   clustering_method_row = "ward.D2",
                   show_row_names = TRUE,
                   show_column_names = FALSE,
                   row_split = 4,
                   column_split = 3,
                   top_annotation = subtype,
                   right_annotation = cell_type,
                   show_row_dend = FALSE)
draw(heatmap)

#ecotype survival KM plot####
clusters <- clusters[,1:2]
clusters_surv <- merge(combined_standard, clusters, by = "sample_id")
clusters_surv <- clusters_surv %>%
  mutate(across(7:10, as.numeric))%>%
  dplyr::select(-c(14:55))

#modify er status if you want er-subtype specific plots. Comment out for full survial curves
clusters_surv_subtype <- filter(clusters_surv, Subtype == "PN")
sample_sizes <- clusters_surv_subtype %>%
  group_by(clusters) %>%
  dplyr::summarize(n = n())
sample_sizes
legend_labels <- sample_sizes %>%
  mutate(label = paste0(clusters, " (n = ", n, ")")) %>%
  pull(label)
km_fit <- survfit(Surv(DFS_years, DFS_event) ~ clusters, data = clusters_surv_subtype)
logrank <- survdiff(Surv(DFS_years, DFS_event) ~ clusters, data = clusters_surv_subtype)
plot <- survfit2(Surv(DFS_years, DFS_event) ~ clusters, 
                 data = clusters_surv_subtype) %>% 
  ggsurvfit() +
  labs(x = "Years", y = "Survival probability", color = "Ecotype") +
  geom_text(x = Inf, y = Inf, 
            label = paste("Logrank P =",  format(logrank$pvalue, digits = 2, scientific = TRUE)),
            hjust = 1, vjust = 1, 
            size = 4) +
  ggtitle("DFS, Proneural")+
  ylim(0, 1)+
  theme_minimal()+
  theme(
    axis.title = element_text(size = 14), 
    axis.text = element_text(size = 12), 
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12),
    plot.title = element_text(size = 16)  
  ) +
  scale_color_manual(values = c("1" = "#ffd97d", 
                                "2" = "#FFADAD", 
                                "3" = "#A0C4FF"),
                     labels = legend_labels) +
  scale_x_continuous(limits = c(0, 7),
                     breaks = seq(0, 10, by = 1))  



plot$layers[[1]]$aes_params$size <- 1

print(plot) 

