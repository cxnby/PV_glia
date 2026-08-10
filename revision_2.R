##### Libraries ---------
# Add at the top of your script, before any plotting:
library(showtext)
showtext_auto()
font_add("Arial", "/System/Library/Fonts/Supplemental/Arial.ttf")
# Load the library
library(extrafont)

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(c("BiocNeighbors", "Biobase", "BiocGenerics"))
# devtools::install_github("jinworks/CellChat")

# Import system fonts (this may take a few minutes)
# You only need to run this command once per machine
# font_import()

# Load the imported fonts for the current R session
loadfonts()



library(Seurat)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(dplyr)
library(clusterProfiler)
library(patchwork)
library(ggplot2)
library(CellChat)
library(circlize)
library(NMF)
library(ggalluvial)
library(patchwork)  # for `+` operator combining plots
library(reticulate)
library(future)
library(pheatmap)
library(grid)
library(gtools)
library(pals)
library(tidyr)
library(reshape2)
library(stringr)
library(purrr)
library(readr)
library(tibble)
library(scales)
library(RColorBrewer)
library(gridtext)
library(DOSE)                 # Disease Ontology (DO)



library("ReactomePA")
library(speckle)
library(biomaRt)
library(msigdbr)
library(GSEABase)


library(pheatmap)
library(ComplexHeatmap)

library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(tibble)
library(purrr)
library(patchwork)
library(ggnewscale)  # multiple fill scales



#py_install("umap-learn")
options(stringsAsFactors = FALSE)
mem.maxVSize(vsize = Inf)
select <- dplyr::select

########################### PNN DATA ###################### ----------------
##### Set everything up and load the PNN data --------
setwd("/Users/cnbr/Downloads/Seurat_CellChat")

pv_cortical_filtered_pnn <- readRDS("~/Downloads/Seurat_CellChat/pv_cortical_subclustered_final_pnn.rds")

if (!dir.exists("revision/cellchat_new")) {
  dir.create("revision/cellchat_new", recursive = TRUE)
}

setwd("~/Downloads/Seurat_CellChat/revision/cellchat_new")

##### Quick reorder 
pv_cortical_filtered_pnn <- RenameIdents(
  pv_cortical_filtered_pnn,
  "Layer 6 corticothalamic neurons" = "Layer 6a corticothalamic neurons",
  "Atypical excitatory neurons"     = "Layer 6b neurons",
  "Newly mature oligodendrocytes" = "Newly formed oligodendrocytes"
)

new_order <- c(
  "Layer 2/3 IT neurons",
  "Layer 4 sensory neurons",
  "Layer 5a IT neurons",
  "Layer 5b PT neurons",
  "Layer 5/6 IT neurons",
  "Layer 6a corticothalamic neurons",   # renamed
  "Layer 6b neurons",                   # renamed + moved here
  "Deep-layer extratelencephalic neurons",
  "Corticospinal neurons (Type I)",
  "Corticospinal neurons (Type II)",
  "PV+ interneurons",
  "SST+ interneurons",
  "VIP+ interneurons",
  "Astrocytes",
  "Oligodendrocyte precursor cells",
  "Newly formed oligodendrocytes",
  "Perineuronal oligodendrocytes",
  "Myelinating oligodendrocytes",
  "Microglia",
  "Endothelial cells",
  "Leptomeningeal cells",
  "Meningeal fibroblasts"
)

pv_cortical_filtered_pnn@active.ident <- factor(pv_cortical_filtered_pnn@active.ident, levels = new_order)

##### Plot UMAP with defined colors --------
cluster_ids <- levels(Idents(pv_cortical_filtered_pnn))
colors <- setNames(
  c(
    "#6B5B95",  # 1 deep lavender
    "#45B8AC",  # 2 mellow aqua
    "#955251",  # 3 dusty rose
    "#4E84C4",  # 4 soft steel blue
    "#B565A7",  # 5 soft mauve
    "#88B04B",  # 6 muted chartreuse
    "#7B6888",  # 10 muted indigo
    "#C3447A",  # 7 muted fuchsia
    "#009B77",  # 8 subdued teal
    "#EFC050",  # 9 warm mustard
    "#7FCDCD",  # 11 misty cyan
    "#DD4124",  # 12 brick red
    "#5B5EA6",  # 13 slate purple
    "#E07A5F",  # 14 muted coral
    "#4BACC6",  # 15 soft cerulean
    "#E8A0BF",  # 16 soft blush
    "#ea8a33",  # some orange for perineuronal oligos
    "#9B2335",  # 17 faded burgundy
    "#C17BAE",  # 18 dusty orchid
    "#DECF3F",  # 19 olive gold
    "#789262",  # 20 sage green
    "#BC243C"   # 21 dark raspberry
  ),
  cluster_ids
)



##### Differential expression analysis -------
# Load required libraries
library(ggplot2)
library(ggrepel)

# Parameters
wt_label  <- "PV-Cre/tdTom"
ko_label  <- "PV-Cre/tdTom/Dnmt1 loxP2"

# Check exact genotype names
table(pv_cortical_filtered_pnn$Genotype)

# Create subsets
wt_cortical_pnn <- subset(pv_cortical_filtered_pnn, subset = Genotype == "PV-Cre/tdTom")
ko_cortical_pnn <- subset(pv_cortical_filtered_pnn, subset = Genotype == "PV-Cre/tdTom/Dnmt1 loxP2")

# Validate subsets and check cell counts
dim(wt_cortical_pnn)  
dim(ko_cortical_pnn)

# For wt object
table(Idents(wt_cortical_pnn))

# For ko object
table(Idents(ko_cortical_pnn))

# 1. Get cluster counts for WT and KO
wt_cortical_cluster_counts_pnn <- as.data.frame(table(Idents(wt_cortical_pnn)))
colnames(wt_cortical_cluster_counts_pnn) <- c("Cluster", "CellNumber_WT")

ko_cortical_cluster_counts_pnn <- as.data.frame(table(Idents(ko_cortical_pnn)))
colnames(ko_cortical_cluster_counts_pnn) <- c("Cluster", "CellNumber_KO")

# 2. Get total number of cells in each dataset
total_wt_cortical_pnn <- ncol(wt_cortical_pnn)  # or dim(wt)[1]
total_ko_cortical_pnn <- ncol(ko_cortical_pnn)  # or dim(ko)[1]


##### Gfap expression in astroyctes ------
gene_name <- "Gfap"

# Expression matrices
expr_mat_wt <- GetAssayData(wt_cortical_pnn, layer = "data")
expr_mat_ko <- GetAssayData(ko_cortical_pnn, layer = "data")

# Counts of Gfap+ cells within Astrocytes cluster
n_wt_Astro <- sum(expr_mat_wt[gene_name, WhichCells(wt_cortical_pnn, idents = "Astrocytes")] > 0)
n_ko_Astro <- sum(expr_mat_ko[gene_name, WhichCells(ko_cortical_pnn, idents = "Astrocytes")] > 0)

# Proportions normalized to total cells per genotype
prop_wt_Astro <- n_wt_Astro / total_wt_cortical_pnn
prop_ko_Astro <- n_ko_Astro / total_ko_cortical_pnn

# Plotting data
plot_data <- data.frame(
  Cell_Type = rep("Astrocyte", 2),
  Condition = c("WT", "KO"),
  Proportion = c(prop_wt_Astro, prop_ko_Astro)
)

plot_data$Condition <- factor(plot_data$Condition, levels = c("WT", "KO"))

# Statistical test
astro_test <- prop.test(c(n_wt_Astro, n_ko_Astro), c(total_wt_cortical_pnn, total_ko_cortical_pnn))
p_value_adj <- astro_test$p.value  # single comparison, no FDR correction needed

annotation_data <- plot_data |>
  dplyr::summarise(y_pos = max(Proportion) + 0.002) |>
  dplyr::mutate(
    p_value = p_value_adj,
    significance = dplyr::case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE            ~ "ns"
    )
  )

# Plot
p1 <- ggplot(plot_data, aes(x = Cell_Type, y = Proportion, fill = Condition)) +
  geom_bar(
    stat = "identity",
    position = position_dodge2(width = 0.8, padding = 0.15, preserve = "single"),
    alpha = 1,
    width = 0.65
  ) +
  scale_fill_manual(values = c("WT" = "#c6c6c6", "KO" = "#7FCDCD")) +
  geom_text(
    data = annotation_data,
    aes(x = "Astrocyte", y = y_pos + 0.001, label = paste("p =", format.pval(p_value, digits = 2))),
    inherit.aes = FALSE, hjust = 0.5, vjust = 1, size = 4
  ) +
  labs(
    x = "Cell Type",
    y = "Proportion of Gfap+ cells (of total)",
    fill = "Genotype"
  ) +
  scale_y_continuous(labels = label_percent(accuracy = 0.1),
                     expand = expansion(mult = c(0, 0.05))) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 14, colour = "black"),
    axis.text.y = element_text(size = 14, colour = "black"),
    axis.title = element_text(size = 16, face = "bold"),
    plot.title = element_text(size = 18, hjust = 0.5, face = "bold"),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    legend.position = "top",
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 14)
  )

print(p1)

ggsave(
  "Gfap-population_Astrocyte.pdf",
  plot = p1,
  width = 5,
  height = 6,
  dpi = 600,
  bg = "transparent"
)

##### Gfap expression normalized to all cells ------------

gene_name <- "Gfap"

# Expression matrices
expr_mat_wt <- GetAssayData(wt_cortical_pnn, layer = "data")
expr_mat_ko <- GetAssayData(ko_cortical_pnn, layer = "data")

# Counts of Gfap+ cells across ALL cells (not restricted to Astrocytes cluster)
n_wt_Gfap_all <- sum(expr_mat_wt[gene_name, ] > 0)
n_ko_Gfap_all <- sum(expr_mat_ko[gene_name, ] > 0)

# Proportions normalized to total cells per genotype
prop_wt_Gfap_all <- n_wt_Gfap_all / total_wt_cortical_pnn
prop_ko_Gfap_all <- n_ko_Gfap_all / total_ko_cortical_pnn

# Plotting data
plot_data <- data.frame(
  Cell_Type = rep("All cells", 2),
  Condition = c("WT", "KO"),
  Proportion = c(prop_wt_Gfap_all, prop_ko_Gfap_all)
)

plot_data$Condition <- factor(plot_data$Condition, levels = c("WT", "KO"))

# Statistical test
gfap_all_test <- prop.test(c(n_wt_Gfap_all, n_ko_Gfap_all), c(total_wt_cortical_pnn, total_ko_cortical_pnn))
p_value_adj <- gfap_all_test$p.value  # single comparison, no FDR correction needed

annotation_data <- plot_data |>
  dplyr::summarise(y_pos = max(Proportion) + 0.002) |>
  dplyr::mutate(
    p_value = p_value_adj,
    significance = dplyr::case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE            ~ "ns"
    )
  )

# Plot
p2 <- ggplot(plot_data, aes(x = Cell_Type, y = Proportion, fill = Condition)) +
  geom_bar(
    stat = "identity",
    position = position_dodge2(width = 0.8, padding = 0.15, preserve = "single"),
    alpha = 1,
    width = 0.65
  ) +
  scale_fill_manual(values = c("WT" = "#c6c6c6", "KO" = "#7FCDCD")) +
  geom_text(
    data = annotation_data,
    aes(x = "All cells", y = y_pos + 0.001, label = paste("p =", format.pval(p_value, digits = 2))),
    inherit.aes = FALSE, hjust = 0.5, vjust = 1, size = 4
  ) +
  labs(
    x = "Cell Type",
    y = "Proportion of Gfap+ cells (of total)",
    fill = "Genotype"
  ) +
  scale_y_continuous(labels = label_percent(accuracy = 0.1),
                     expand = expansion(mult = c(0, 0.05))) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 14, colour = "black"),
    axis.text.y = element_text(size = 14, colour = "black"),
    axis.title = element_text(size = 16, face = "bold"),
    plot.title = element_text(size = 18, hjust = 0.5, face = "bold"),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    legend.position = "top",
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 14)
  )

print(p2)

ggsave(
  "Gfap-population_AllCells.pdf",
  plot = p2,
  width = 5,
  height = 6,
  dpi = 600,
  bg = "transparent"
)

##### Gfap statistics ----
##### Summary table: Total cells, Astrocyte cells, and Gfap+ proportions (WT vs KO) 
library(dplyr)
library(gridExtra)
library(grid)

# ---- Astrocyte-specific counts (reuses n_wt_Astro / n_ko_Astro from earlier script) 
n_wt_Astro <- sum(expr_mat_wt[gene_name, WhichCells(wt_cortical_pnn, idents = "Astrocytes")] > 0)
n_ko_Astro <- sum(expr_mat_ko[gene_name, WhichCells(ko_cortical_pnn, idents = "Astrocytes")] > 0)

total_astro_wt <- length(WhichCells(wt_cortical_pnn, idents = "Astrocytes"))
total_astro_ko <- length(WhichCells(ko_cortical_pnn, idents = "Astrocytes"))

prop_wt_Astro <- n_wt_Astro / total_wt_cortical_pnn
prop_ko_Astro <- n_ko_Astro / total_ko_cortical_pnn

# ---- All-cell Gfap+ counts (reuses n_wt_Gfap_all / n_ko_Gfap_all from earlier script) 
n_wt_Gfap_all <- sum(expr_mat_wt[gene_name, ] > 0)
n_ko_Gfap_all <- sum(expr_mat_ko[gene_name, ] > 0)

prop_wt_Gfap_all <- n_wt_Gfap_all / total_wt_cortical_pnn
prop_ko_Gfap_all <- n_ko_Gfap_all / total_ko_cortical_pnn

# ---- Statistical tests for both comparisons 
astro_test <- prop.test(c(n_wt_Astro, n_ko_Astro), c(total_wt_cortical_pnn, total_ko_cortical_pnn))
all_test   <- prop.test(c(n_wt_Gfap_all, n_ko_Gfap_all), c(total_wt_cortical_pnn, total_ko_cortical_pnn))

# ---- Build the summary data frame 
summary_table <- data.frame(
  Metric = c(
    "Total cells (n)",
    "Astrocyte cells (n)",
    "Gfap+ cells within Astrocytes (n)",
    "Gfap+ proportion within Astrocytes (%)",
    "Gfap+ cells across all cells (n)",
    "Gfap+ proportion across all cells (%)"
  ),
  WT = c(
    format(total_wt_cortical_pnn, big.mark = ","),
    format(total_astro_wt, big.mark = ","),
    format(n_wt_Astro, big.mark = ","),
    sprintf("%.2f", prop_wt_Astro * 100),
    format(n_wt_Gfap_all, big.mark = ","),
    sprintf("%.2f", prop_wt_Gfap_all * 100)
  ),
  KO = c(
    format(total_ko_cortical_pnn, big.mark = ","),
    format(total_astro_ko, big.mark = ","),
    format(n_ko_Astro, big.mark = ","),
    sprintf("%.2f", prop_ko_Astro * 100),
    format(n_ko_Gfap_all, big.mark = ","),
    sprintf("%.2f", prop_ko_Gfap_all * 100)
  ),
  p_value = c(
    "-",
    "-",
    "-",
    format.pval(astro_test$p.value, digits = 3),
    "-",
    format.pval(all_test$p.value, digits = 3)
  ),
  stringsAsFactors = FALSE
)

print(summary_table)

# ---- Save as CSV 
write.csv(summary_table, file = "Gfap_summary_table_WT_vs_KO.csv", row.names = FALSE)

# ---- Render as a formatted table and save as PDF 
table_theme <- ttheme_default(
  core = list(
    fg_params = list(hjust = 0, x = 0.05, fontsize = 11),
    bg_params = list(fill = c("#f7f7f7", "#ffffff"))
  ),
  colhead = list(
    fg_params = list(fontface = "bold", fontsize = 12),
    bg_params = list(fill = "#d9d9d9")
  )
)

table_grob <- tableGrob(summary_table, rows = NULL, theme = table_theme)

pdf("Gfap_summary_table_WT_vs_KO.pdf", width = 9, height = 3.5)
grid.draw(table_grob)
dev.off()
##### FULL DATASET WT/KO UMAP — single overlaid panel ----
pv_cortical_filtered_pnn$Genotype_short <- recode(pv_cortical_filtered_pnn$Genotype,
                                                  "PV-Cre/tdTom"              = "Ctrl",
                                                  "PV-Cre/tdTom/Dnmt1 loxP2" = "KO"
)

p_wt_ko_overlay <- DimPlot(pv_cortical_filtered_pnn,
                           reduction = "umap",
                           group.by  = "Genotype_short",
                           cols      = c("Ctrl" = "#c6c6c6", "KO" = "#7FCDCD"),
                           pt.size   = 0.1,
                           alpha     = 1,
                           shuffle   = TRUE) +  # avoids one group hiding under the other
  ggtitle(NULL) +
  theme(
    legend.position = "right",
    legend.text     = element_text(size = 13, face = "bold")
  )

pdf("FullDataset_UMAP_WT_KO_overlay.pdf",
    width = 8, height = 6)
p_wt_ko_overlay
dev.off()



######## ATACseq
suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(ggplot2)
  library(ragg)
})

# -------------------- Paths 
data_dir <- "/Users/cnbr/Downloads/Seurat_scATAC-seq"
objmale <- readRDS(
  file.path(data_dir, "scATAC_07_male_integrated.rds")
)

# -------------------- Load curated annotated object 

stopifnot("cell_type" %in% colnames(objmale@meta.data))

# -------------------- Original ATAC label order 
# 18 categories present after curated integration; Oligodendrocytes (collapsed
# from NFOL/Perineuronal/Myelinating) takes NFOL's original color (#E8A0BF),
# nothing else shifted or renamed.
atacneworder <- c(
  "Layer 2/3 IT neurons",
  "Layer 4 sensory neurons",
  "Layer 5a IT neurons",
  "Layer 5b PT neurons",
  "Layer 5/6 IT neurons",
  "Layer 6a corticothalamic neurons",
  "Layer 6b neurons",
  "Deep-layer extratelencephalic neurons",
  "Corticospinal neurons (Type II)",
  "PV+ interneurons",
  "SST+ interneurons",
  "VIP+ interneurons",
  "Astrocytes",
  "Oligodendrocyte precursor cells",
  "Oligodendrocytes",
  "Microglia",
  "Endothelial cells",
  "Meningeal fibroblasts"
)

hexcolors <- c(
  "#6B5B95",  # Layer 2/3 IT neurons
  "#45B8AC",  # Layer 4 sensory neurons
  "#955251",  # Layer 5a IT neurons
  "#4E84C4",  # Layer 5b PT neurons
  "#B565A7",  # Layer 5/6 IT neurons
  "#88B04B",  # Layer 6a corticothalamic neurons
  "#7B6888",  # Layer 6b neurons
  "#C3447A",  # Deep-layer extratelencephalic neurons
  "#EFC050",  # Corticospinal neurons (Type II) -- keeps its own original color
  "#7FCDCD",  # PV+ interneurons
  "#DD4124",  # SST+ interneurons
  "#5B5EA6",  # VIP+ interneurons
  "#E07A5F",  # Astrocytes
  "#4BACC6",  # Oligodendrocyte precursor cells
  "#E8A0BF",  # Oligodendrocytes -- takes NFOL's color, per collapse rule
  "#C17BAE",  # Microglia
  "#DECF3F",  # Endothelial cells
  "#BC243C"   # Meningeal fibroblasts
)

stopifnot(length(hexcolors) == length(atacneworder))  # guard against future mismatch

rnacolorsatac <- setNames(hexcolors, atacneworder)

# Keep only cell types that actually exist after your curated integration.
atac_types <- atacneworder[
  atacneworder %in% unique(as.character(objmale$cell_type))
]

objmale$cell_type <- factor(
  objmale$cell_type,
  levels = atac_types
)

ataccolors <- rnacolorsatac[atac_types]

# These should be zero because cluster 18 was removed before saveRDS().
message("Leptomeningeal cells: ",
        sum(objmale$cell_type == "Leptomeningeal cells", na.rm = TRUE))

message("Corticospinal Type I cells: ",
        sum(objmale$cell_type == "Corticospinal neurons Type I", na.rm = TRUE))

print(sort(table(objmale$cell_type), decreasing = TRUE))

# -------------------- Plot: same approach as your script 
patacumap <- DimPlot(
  objmale,
  group.by = "cell_type",
  cols = ataccolors,
  pt.size = 0.1,
  raster = FALSE,
  label = FALSE
) +
  ggtitle("snATAC-seq (projected), male") +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title = element_text(face = "bold"),
    legend.position = "right",
    legend.text = element_text(size = 9),
    legend.key.size = grid::unit(0.4, "cm")
  ) +
  guides(
    colour = guide_legend(
      ncol = 1,
      override.aes = list(size = 3.2, alpha = 1)
    )
  )

print(patacumap)

# -------------------- Save 
ggsave(
  file.path(data_dir, "snATAC_male_projected_umap.pdf"),
  patacumap,
  width = 12,
  height = 9,
  device = grDevices::pdf,
  useDingbats = FALSE
)

ggsave(
  file.path(data_dir, "snATAC_male_projected_umap.png"),
  patacumap,
  width = 12,
  height = 9,
  dpi = 600,
  device = ragg::agg_png,
  background = "white"
)



#### Oligo dotplot v2 ----

oligo_clusters <- c(
  "Oligodendrocyte precursor cells",
  "Newly formed oligodendrocytes",
  "Perineuronal oligodendrocytes",
  "Myelinating oligodendrocytes"
)

pv_oligo_only <- subset(pv_cortical_filtered_pnn, idents = oligo_clusters)
pv_oligo_only@active.ident <- factor(pv_oligo_only@active.ident, levels = oligo_clusters)

cat("Cell counts per oligodendrocyte cluster:\n")
print(table(Idents(pv_oligo_only)))

# ---- Lean curated panel: core lineage markers only
curated_markers <- c(
  "Pdgfra", "Cspg4",
  "Enpp6", "Bcas1",
  "Pdgfc",
  "Mbp", "Mog", "Plp1"
)

# ---- PNOL-distinguishing genes (statistically significant vs MOL and/or NFOL)
pnol_distinguishing <- c("Anln", "Ctnna3", "Trf", "Rftn1", "Kif21b")

all_features <- rownames(pv_oligo_only[["RNA"]])
match_to_features <- function(gvec, features) {
  lut <- setNames(features, toupper(features))
  matched <- lut[toupper(gvec)]
  unique(unname(matched[!is.na(matched)]))
}
missing_report <- function(requested, matched) {
  req_up <- toupper(requested)
  feat_up <- toupper(matched)
  missing <- requested[!(req_up %in% feat_up)]
  if (length(missing) > 0) {
    message("Not found in object: ", paste(missing, collapse = ", "))
  }
}

curated_matched <- match_to_features(curated_markers, all_features)
missing_report(curated_markers, curated_matched)
curated_final <- curated_markers[toupper(curated_markers) %in% toupper(curated_matched)]

pnol_matched <- match_to_features(pnol_distinguishing, all_features)
missing_report(pnol_distinguishing, pnol_matched)
pnol_final <- pnol_distinguishing[toupper(pnol_distinguishing) %in% toupper(pnol_matched)]

oligo_marker_panel_combined <- c(curated_final, pnol_final)

# ---- Order genes by peak cluster expression
avg_exp <- AverageExpression(
  pv_oligo_only,
  features = oligo_marker_panel_combined,
  assays = "RNA",
  slot = "data"
)$RNA

avg_exp <- avg_exp[, oligo_clusters, drop = FALSE]

peak_cluster <- apply(avg_exp, 1, function(x) oligo_clusters[which.max(x)])
peak_value <- apply(avg_exp, 1, max)

gene_order_df <- data.frame(
  gene = rownames(avg_exp),
  peak_cluster = factor(peak_cluster, levels = oligo_clusters),
  peak_value = peak_value
) %>%
  arrange(peak_cluster, desc(peak_value))

oligo_marker_panel_final <- gene_order_df$gene

cat("Final lean marker panel (", length(oligo_marker_panel_final), " genes):\n", sep = "")
print(gene_order_df)

# ---- Blue gradient, consistent with your other DotPlots
pal_blues <- colorRampPalette(c(
  "#f7fbff", "#deebf7", "#c6dbef", "#9ecae1",
  "#6baed6", "#4292c6", "#2171b5", "#2166ac", "#084594"
))(100)

cmin <- -1.5
cmax <- 2.5

p_oligo_curated <- DotPlot(
  pv_oligo_only,
  features = oligo_marker_panel_final,
  dot.scale = 5,
  col.min = cmin,
  col.max = cmax
) +
  scale_colour_gradientn(
    colours = pal_blues,
    limits = c(cmin, cmax),
    oob = scales::squish,
    name = "Avg exp scaled"
  ) +
  labs(title = "Oligodendrocyte lineage marker panel", x = NULL, y = NULL) +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, face = "italic"),
    axis.text.y = element_text(color = "black"),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

print(p_oligo_curated)

ggsave(
  filename = "oligodendrocyte_lean_markers_dotplot.pdf",
  plot = p_oligo_curated,
  width = 9, height = 3.5
)




#### Oligo panel ----
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)

#### Setup: process pv_oligo_only from scratch

pv_oligo_only <- NormalizeData(pv_oligo_only)
pv_oligo_only <- FindVariableFeatures(pv_oligo_only, nfeatures = 2000)
pv_oligo_only <- ScaleData(pv_oligo_only)
pv_oligo_only <- RunPCA(pv_oligo_only, npcs = 15)
pv_oligo_only <- FindNeighbors(pv_oligo_only, dims = 1:10)
pv_oligo_only <- RunUMAP(pv_oligo_only, dims = 1:10)

pv_oligo_only$Genotype_short <- recode(pv_oligo_only$Genotype,
                                       "PV-Cre/tdTom"              = "Ctrl",
                                       "PV-Cre/tdTom/Dnmt1 loxP2" = "KO"
)

DefaultAssay(pv_oligo_only) <- "RNA"

oligo_colors <- c(
  "Oligodendrocyte precursor cells" = "#4BACC6",
  "Perineuronal oligodendrocytes"   = "#ea8a33",
  "Newly formed oligodendrocytes"   = "#E8A0BF",
  "Myelinating oligodendrocytes"    = "#9B2335"
)

short_names <- c(
  "Oligodendrocyte precursor cells" = "OPC",
  "Perineuronal oligodendrocytes"   = "PNOL",
  "Newly formed oligodendrocytes"   = "NFOL",
  "Myelinating oligodendrocytes"    = "MOL"
)

pal_blues <- colorRampPalette(c(
  "#f7fbff", "#deebf7", "#c6dbef", "#9ecae1",
  "#6baed6", "#4292c6", "#2171b5", "#2166ac", "#084594"
))(100)

#### Panel A: color-coded overview UMAP (all clusters visible)

p_overview <- DimPlot(pv_oligo_only, reduction = "umap", pt.size = 0.3) +
  scale_colour_manual(values = oligo_colors) +
  ggtitle("Oligodendrocyte lineage sub-types") +
  theme_classic() +
  theme(
    plot.title      = element_text(hjust = 0.5, face = "bold", size = 11),
    legend.position = "right",
    axis.line       = element_line(color = "black", linewidth = 0.4)
  )

#### Panel B: per-cluster Ctrl/KO/Other split panels (4 plots)

cluster_plot_list <- list()

for (cl in oligo_clusters) {
  
  pv_oligo_only$plot_group <- ifelse(
    as.character(Idents(pv_oligo_only)) == cl,
    as.character(pv_oligo_only$Genotype_short),
    "Other"
  )
  
  pv_oligo_only$plot_group <- factor(pv_oligo_only$plot_group,
                                     levels = c("Other", "Ctrl", "KO"))
  
  p <- DimPlot(pv_oligo_only,
               reduction = "umap",
               group.by  = "plot_group",
               order     = c("KO", "Ctrl", "Other"),
               pt.size   = 0.2,
               alpha     = 1) +
    scale_colour_manual(
      values = c("Other" = "grey96",
                 "Ctrl"  = "#c6c6c6",
                 "KO"    = "#7FCDCD"),
      labels = c("Other" = "Other clusters",
                 "Ctrl"  = "Ctrl",
                 "KO"    = "KO")
    ) +
    ggtitle(cl) +
    theme_classic() +
    theme(
      plot.title      = element_text(hjust = 0.5, face = "bold", size = 11),
      legend.position = "right",
      axis.line       = element_line(color = "black", linewidth = 0.4)
    ) +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))
  
  cluster_plot_list[[cl]] <- p
}

#### Panel C: Enpp6 and Pdgfc expression UMAPs

feat_genes <- c("Enpp6", "Pdgfc")

feat_plot_list <- lapply(feat_genes, function(g) {
  FeaturePlot(pv_oligo_only,
              features  = g,
              reduction = "umap",
              pt.size   = 0.3,
              order     = TRUE) +
    scale_colour_gradientn(colours = pal_blues) +
    ggtitle(g) +
    theme_classic() +
    theme(
      plot.title      = element_text(hjust = 0.5, face = "bold", size = 11, face2 = "italic"),
      legend.position = "right",
      axis.line       = element_line(color = "black", linewidth = 0.4)
    )
})
names(feat_plot_list) <- feat_genes

#### Combine into one Ext Data Fig 8 patchwork

p_extdata8 <- (
  p_overview +
    cluster_plot_list[[1]] + cluster_plot_list[[2]] +
    cluster_plot_list[[3]] + cluster_plot_list[[4]] +
    feat_plot_list[["Enpp6"]] + feat_plot_list[["Pdgfc"]]
) +
  plot_layout(ncol = 2) +
  plot_annotation(
    title = "Oligodendrocyte lineage: clusters, genotype, and marker expression",
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 15))
  )

pdf("ExtDataFig8_OligoUMAP_Full.pdf", width = 11, height = 13)
print(p_extdata8)
dev.off()

cat("Saved: ExtDataFig8_OligoUMAP_Full.pdf\n")

### UNUSED BELOW vvvvvvv --------
#### Oligo dotplot v1----
# ---- Define the oligodendrocyte-lineage clusters only 
oligo_clusters <- c(
  "Oligodendrocyte precursor cells",
  "Newly formed oligodendrocytes",
  "Perineuronal oligodendrocytes",
  "Myelinating oligodendrocytes"
)

# ---- Subset the Seurat object to just these clusters 
pv_oligo_only <- subset(pv_cortical_filtered_pnn, idents = oligo_clusters)
pv_oligo_only@active.ident <- factor(pv_oligo_only@active.ident, levels = oligo_clusters)

cat("Cell counts per oligodendrocyte cluster:\n")
print(table(Idents(pv_oligo_only)))

# ---- Curated canonical marker gene panel, ordered by biological stage 
curated_markers <- c(
  # Pan-lineage / progenitor (OPC) markers
  "Pdgfra", "Cspg4", "Olig1", "Olig2", "Sox10",
  # Newly formed / premyelinating oligodendrocyte markers
  "Enpp6", "Bcas1", "Tcf7l2", "Gpr17",
  # Perineuronal oligodendrocyte markers (retain progenitor-like receptor signaling)
  "Pdgfc", "Pdgfrb", "Cldn10", "Cxcr4",
  # Mature myelinating oligodendrocyte markers
  "Mbp", "Mog", "Plp1", "Mag", "Cnp", "Mobp"
)

# ---- Verify which curated genes are present in the object (case-insensitive matching) 
all_features <- rownames(pv_oligo_only[["RNA"]])
match_to_features <- function(gvec, features) {
  lut <- setNames(features, toupper(features))
  matched <- lut[toupper(gvec)]
  unique(unname(matched[!is.na(matched)]))
}
missing_report <- function(requested, matched) {
  req_up <- toupper(requested)
  feat_up <- toupper(matched)
  missing <- requested[!(req_up %in% feat_up)]
  if (length(missing) > 0) {
    message("Not found in object: ", paste(missing, collapse = ", "))
  }
}

curated_matched <- match_to_features(curated_markers, all_features)
missing_report(curated_markers, curated_matched)
curated_final <- curated_markers[toupper(curated_markers) %in% toupper(curated_matched)]

# ---- Data-driven marker discovery restricted to oligodendrocyte clusters 
oligo_markers <- FindAllMarkers(
  pv_oligo_only,
  min.pct = 0.25,        # Gene must be detected in >=25% of cells
  logfc.threshold = 0.25,# Minimum fold-change threshold
  only.pos = FALSE,      # Keep both up- and down-regulated markers
  test.use = "wilcox"    # Default Wilcoxon Rank Sum test
)

write.csv(oligo_markers, file = "oligodendrocyte_cluster_marker_genes.csv", row.names = FALSE)

top5_data_driven <- oligo_markers %>%
  group_by(cluster) %>%
  filter(p_val_adj < 0.05) %>%
  filter(!grepl("^ENSMUSG", gene)) %>%
  filter(!grepl("^Gm", gene)) %>%
  filter(!grepl("Rik$", gene)) %>%
  slice_max(avg_log2FC, n = 5) %>%
  pull(gene) %>%
  unique()

# ---- Combine curated + data-driven, removing redundancy (case-insensitive) 
curated_upper <- toupper(curated_final)
data_driven_unique <- top5_data_driven[!(toupper(top5_data_driven) %in% curated_upper)]

oligo_marker_panel_combined <- c(curated_final, data_driven_unique)

# ---- Determine each gene's "peak" cluster and expression strength for sorting 
avg_exp <- AverageExpression(
  pv_oligo_only,
  features = oligo_marker_panel_combined,
  assays = "RNA",
  slot = "data"
)$RNA

# Ensure column order matches factor levels
avg_exp <- avg_exp[, oligo_clusters, drop = FALSE]

# For each gene, find which cluster has max average expression, and that max value
peak_cluster <- apply(avg_exp, 1, function(x) oligo_clusters[which.max(x)])
peak_value <- apply(avg_exp, 1, max)

gene_order_df <- data.frame(
  gene = rownames(avg_exp),
  peak_cluster = factor(peak_cluster, levels = oligo_clusters),
  peak_value = peak_value
) %>%
  arrange(peak_cluster, desc(peak_value))

oligo_marker_panel_final <- gene_order_df$gene

cat("Final combined marker panel, sorted by peak cluster (", length(oligo_marker_panel_final), " genes):\n", sep = "")
print(gene_order_df)

# ---- Blue color gradient matching your existing DotPlots 
pal_blues <- colorRampPalette(c(
  "#f7fbff", "#deebf7", "#c6dbef", "#9ecae1",
  "#6baed6", "#4292c6", "#2171b5", "#2166ac", "#084594"
))(100)

cmin <- -1.5
cmax <- 2.5

# ---- Build the DotPlot with the sorted combined panel 
p_oligo_curated <- DotPlot(
  pv_oligo_only,
  features = oligo_marker_panel_final,
  dot.scale = 5,
  col.min = cmin,
  col.max = cmax
) +
  scale_colour_gradientn(
    colours = pal_blues,
    limits = c(cmin, cmax),
    oob = scales::squish,
    name = "Avg exp scaled"
  ) +
  labs(title = "Oligodendrocyte lineage marker panel", x = NULL, y = NULL) +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, face = "italic"),
    axis.text.y = element_text(color = "black"),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

print(p_oligo_curated)

# ---- Save as high-res PDF 
ggsave(
  filename = "oligodendrocyte_curated_markers_dotplot.pdf",
  plot = p_oligo_curated,
  width = 12.5, height = 3.5
)


#### OLIGO UMAP v2 ------
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)

# ---- Blue gradient palette for marker expression 
pal_blues <- colorRampPalette(c(
  "#f7fbff", "#deebf7", "#c6dbef", "#9ecae1",
  "#6baed6", "#4292c6", "#2171b5", "#2166ac", "#084594"
))(100)

oligo_clusters <- c("Oligodendrocyte precursor cells",
                    "Newly formed oligodendrocytes",
                    "Perineuronal oligodendrocytes",
                    "Myelinating oligodendrocytes")

oligo_only <- subset(pv_cortical_filtered_pnn, idents = oligo_clusters)

oligo_only <- NormalizeData(oligo_only)
oligo_only <- FindVariableFeatures(oligo_only, nfeatures = 2000)
oligo_only <- ScaleData(oligo_only)
oligo_only <- RunPCA(oligo_only, npcs = 15)
oligo_only <- FindNeighbors(oligo_only, dims = 1:10)
oligo_only <- RunUMAP(oligo_only, dims = 1:10)

oligo_only$Genotype_short <- recode(oligo_only$Genotype,
                                    "PV-Cre/tdTom"              = "Ctrl",
                                    "PV-Cre/tdTom/Dnmt1 loxP2" = "KO"
)

oligo_colors <- c(
  "Oligodendrocyte precursor cells" = "#4BACC6",
  "Perineuronal oligodendrocytes"   = "#ea8a33",
  "Newly formed oligodendrocytes"   = "#E8A0BF",
  "Myelinating oligodendrocytes"    = "#9B2335"
)

short_names <- c(
  "Oligodendrocyte precursor cells" = "OPC",
  "Perineuronal oligodendrocytes"   = "PNOL",
  "Newly formed oligodendrocytes"   = "NFOL",
  "Myelinating oligodendrocytes"    = "MOL"
)

# 1) Overview panels: cluster-colored UMAP + genotype UMAP
p_oligo_clusters <- DimPlot(oligo_only, reduction = "umap",
                            label = FALSE, repel = TRUE, pt.size  = 0.3) +
  scale_colour_manual(values = oligo_colors) +
  ggtitle("Oligodendrocyte lineage sub-types") +
  theme_classic() +
  theme(
    plot.title      = element_text(hjust = 0.5, face = "bold", size = 11),
    legend.position = "right",
    axis.line       = element_line(color = "black", linewidth = 0.4)
  )

p_oligo_genotype <- DimPlot(oligo_only, reduction = "umap",
                            group.by = "Genotype_short",
                            cols     = c("Ctrl" = "#c6c6c6", "KO" = "#7FCDCD"),
                            pt.size  = 0.3) +
  ggtitle("Oligodendrocyte UMAP: Ctrl vs KO") +
  theme_classic() +
  theme(
    plot.title      = element_text(hjust = 0.5, face = "bold", size = 11),
    legend.position = "right",
    axis.line       = element_line(color = "black", linewidth = 0.4)
  )


# 2) Per-cluster WT/KO panels


cluster_plot_list <- list()

for (cl in oligo_clusters) {
  
  oligo_only$plot_group <- ifelse(
    as.character(Idents(oligo_only)) == cl,
    as.character(oligo_only$Genotype_short),
    "Other"
  )
  
  oligo_only$plot_group <- factor(oligo_only$plot_group,
                                  levels = c("Other", "Ctrl", "KO"))
  
  p <- DimPlot(oligo_only,
               reduction = "umap",
               group.by  = "plot_group",
               order     = c("KO", "Ctrl", "Other"),
               pt.size   = 0.2,
               alpha     = 1) +
    scale_colour_manual(
      values = c("Other" = "grey96",
                 "Ctrl"  = "#c6c6c6",
                 "KO"    = "#7FCDCD"),
      labels = c("Other" = "Other clusters",
                 "Ctrl"  = "Ctrl",
                 "KO"    = "KO")
    ) +
    ggtitle(cl) +
    theme_classic() +
    theme(
      plot.title      = element_text(hjust = 0.5, face = "bold", size = 11),
      legend.position = "right",
      axis.line       = element_line(color = "black", linewidth = 0.4)
    ) +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))
  
  cluster_plot_list[[cl]] <- p
}

p_combined_percluster <- (
  cluster_plot_list[[1]] + cluster_plot_list[[2]] +
    cluster_plot_list[[3]] + cluster_plot_list[[4]]
) +
  plot_layout(ncol = 2) +
  plot_annotation(
    title = "Oligodendrocyte lineage: Ctrl vs KO per sub-type",
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
  )

pdf("OligoOnly_UMAP_percluster_WT_KO.pdf", width = 11, height = 8)
print(p_combined_percluster)
dev.off()


# 3) Expanded marker gene sets (testing additional candidates)

marker_sets <- list(
  OPC  = c("Pdgfra"),
  NFOL = c("Enpp6", "Bcas1", "Tcf7l2", "Gpr17", "Itpr2"),
  PNOL = c("Pdgfc"),
  MOL  = c("Mbp", "Plp1", "Mog", "Mag", "Mobp", "Cnp")
)

genes_in_obj <- rownames(oligo_only)
marker_sets  <- lapply(marker_sets, function(g) intersect(g, genes_in_obj))

DefaultAssay(oligo_only) <- "RNA"

marker_plot_list <- list()

for (stage in names(marker_sets)) {
  genes <- marker_sets[[stage]]
  if (length(genes) == 0) next
  
  for (gene in genes) {
    p <- FeaturePlot(oligo_only,
                     features  = gene,
                     reduction = "umap",
                     pt.size   = 0.3,
                     order     = TRUE) +
      scale_colour_gradientn(colours = pal_blues) +
      ggtitle(gene) +
      theme_classic() +
      theme(
        plot.title      = element_text(hjust = 0.5, face = "bold", size = 11),
        legend.position = "right",
        axis.line       = element_line(color = "black", linewidth = 0.4)
      )
    
    marker_plot_list[[gene]] <- p
  }
}

# 4) Combined patchwork: overview + per-cluster + ALL marker candidates

n_marker_panels <- length(marker_plot_list)
n_total_panels  <- 6 + n_marker_panels   # 2 overview + 4 per-cluster + markers

p_combined_all <- wrap_plots(
  c(list(p_oligo_clusters, p_oligo_genotype),
    cluster_plot_list,
    marker_plot_list),
  ncol = 2
) +
  plot_annotation(
    title = "Oligodendrocyte lineage: clusters, genotype, and candidate marker expression",
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))
  )

pdf("OligoOnly_UMAP_FULL_combined.pdf",
    width = 11, height = 4.2 * ceiling(n_total_panels / 2))
print(p_combined_all)
dev.off()

# 5) DotPlot summary across ALL candidate markers (blue gradient)

dot_genes <- unlist(marker_sets, use.names = FALSE)

p_dotplot <- DotPlot(oligo_only,
                     features = dot_genes,
                     cols     = c("#f7fbff", "#084594")) +
  RotatedAxis() +
  ggtitle("Candidate marker expression across oligodendrocyte sub-types") +
  theme(
    plot.title  = element_text(hjust = 0.5, face = "bold", size = 12),
    axis.text.x = element_text(face = "italic")
  )

pdf("OligoOnly_DotPlot_Markers.pdf", width = 10, height = 5)
print(p_dotplot)
dev.off()
#### Some oligo plot ----
cat("All combined oligo UMAP + marker patchworks saved (no individual panel PDFs).\n")

# ---- Blue gradient palette for marker expression 
pal_blues <- colorRampPalette(c(
  "#f7fbff", "#deebf7", "#c6dbef", "#9ecae1",
  "#6baed6", "#4292c6", "#2171b5", "#2166ac", "#084594"
))(100)

oligo_clusters <- c("Oligodendrocyte precursor cells",
                    "Newly formed oligodendrocytes",
                    "Perineuronal oligodendrocytes",
                    "Myelinating oligodendrocytes")

oligo_only <- subset(pv_cortical_filtered_pnn, idents = oligo_clusters)

oligo_only <- NormalizeData(oligo_only)
oligo_only <- FindVariableFeatures(oligo_only, nfeatures = 2000)
oligo_only <- ScaleData(oligo_only)
oligo_only <- RunPCA(oligo_only, npcs = 15)
oligo_only <- FindNeighbors(oligo_only, dims = 1:10)
oligo_only <- RunUMAP(oligo_only, dims = 1:10)

oligo_only$Genotype_short <- recode(oligo_only$Genotype,
                                    "PV-Cre/tdTom"              = "Ctrl",
                                    "PV-Cre/tdTom/Dnmt1 loxP2" = "KO"
)

oligo_colors <- c(
  "Oligodendrocyte precursor cells" = "#4BACC6",
  "Perineuronal oligodendrocytes"   = "#ea8a33",
  "Newly formed oligodendrocytes"   = "#E8A0BF",
  "Myelinating oligodendrocytes"    = "#9B2335"
)

short_names <- c(
  "Oligodendrocyte precursor cells" = "OPC",
  "Perineuronal oligodendrocytes"   = "PNOL",
  "Newly formed oligodendrocytes"   = "NFOL",
  "Myelinating oligodendrocytes"    = "MOL"
)


# 1) Overview panels: cluster-colored UMAP + genotype UMAP

p_oligo_clusters <- DimPlot(oligo_only, reduction = "umap",
                            label = FALSE, repel = TRUE, pt.size  = 0.3) +
  scale_colour_manual(values = oligo_colors) +
  ggtitle("Oligodendrocyte lineage sub-types") +
  theme_classic() +
  theme(
    plot.title      = element_text(hjust = 0.5, face = "bold", size = 11),
    legend.position = "right",
    axis.line       = element_line(color = "black", linewidth = 0.4)
  )

p_oligo_genotype <- DimPlot(oligo_only, reduction = "umap",
                            group.by = "Genotype_short",
                            cols     = c("Ctrl" = "#c6c6c6", "KO" = "#7FCDCD"),
                            pt.size  = 0.3) +
  ggtitle("Oligodendrocyte UMAP: Ctrl vs KO") +
  theme_classic() +
  theme(
    plot.title      = element_text(hjust = 0.5, face = "bold", size = 11),
    legend.position = "right",
    axis.line       = element_line(color = "black", linewidth = 0.4)
  )


# 2) Per-cluster WT/KO panels


cluster_plot_list <- list()

for (cl in oligo_clusters) {
  
  oligo_only$plot_group <- ifelse(
    as.character(Idents(oligo_only)) == cl,
    as.character(oligo_only$Genotype_short),
    "Other"
  )
  
  oligo_only$plot_group <- factor(oligo_only$plot_group,
                                  levels = c("Other", "Ctrl", "KO"))
  
  p <- DimPlot(oligo_only,
               reduction = "umap",
               group.by  = "plot_group",
               order     = c("KO", "Ctrl", "Other"),
               pt.size   = 0.2,
               alpha     = 1) +
    scale_colour_manual(
      values = c("Other" = "grey96",
                 "Ctrl"  = "#c6c6c6",
                 "KO"    = "#7FCDCD"),
      labels = c("Other" = "Other clusters",
                 "Ctrl"  = "Ctrl",
                 "KO"    = "KO")
    ) +
    ggtitle(cl) +
    theme_classic() +
    theme(
      plot.title      = element_text(hjust = 0.5, face = "bold", size = 11),
      legend.position = "right",
      axis.line       = element_line(color = "black", linewidth = 0.4)
    ) +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))
  
  cluster_plot_list[[cl]] <- p
}

for (cl in oligo_clusters) {
  pdf(paste0("OligoOnly_UMAP_", short_names[cl], "_WT_KO.pdf"),
      width = 6, height = 5)
  print(cluster_plot_list[[cl]])
  dev.off()
}


# 3) Marker gene expression UMAPs (blue gradient)


oligo_markers <- c(
  "Pdgfra" = "Oligodendrocyte precursor cells",
  "Enpp6"  = "Newly formed oligodendrocytes",
  "Pdgfc"  = "Perineuronal oligodendrocytes",
  "Mbp"    = "Myelinating oligodendrocytes"
)

DefaultAssay(oligo_only) <- "RNA"

marker_plot_list <- list()

for (gene in names(oligo_markers)) {
  
  p <- FeaturePlot(oligo_only,
                   features  = gene,
                   reduction = "umap",
                   pt.size   = 0.3,
                   order     = TRUE) +
    scale_colour_gradientn(colours = pal_blues) +
    ggtitle(paste0(gene, " (", oligo_markers[gene], ")")) +
    theme_classic() +
    theme(
      plot.title      = element_text(hjust = 0.5, face = "bold", size = 11),
      legend.position = "right",
      axis.line       = element_line(color = "black", linewidth = 0.4)
    )
  
  marker_plot_list[[gene]] <- p
}

for (gene in names(oligo_markers)) {
  pdf(paste0("OligoOnly_UMAP_Marker_", gene, ".pdf"),
      width = 6, height = 5)
  print(marker_plot_list[[gene]])
  dev.off()
}


# 4) ONE combined patchwork: overview + per-cluster + markers (10 panels)


p_combined_all <- (
  p_oligo_clusters + p_oligo_genotype +
    cluster_plot_list[[1]] + cluster_plot_list[[2]] +
    cluster_plot_list[[3]] + cluster_plot_list[[4]] +
    marker_plot_list[["Pdgfra"]] + marker_plot_list[["Enpp6"]] +
    marker_plot_list[["Pdgfc"]] + marker_plot_list[["Mbp"]]
) +
  plot_layout(ncol = 2) +
  plot_annotation(
    title = "Oligodendrocyte lineage: clusters, genotype, and marker expression",
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))
  )

pdf("OligoOnly_UMAP_FULL_combined.pdf", width = 14, height = 20)
print(p_combined_all)
dev.off()


# 5) DotPlot summary (also switched to blue gradient)


p_dotplot <- DotPlot(oligo_only,
                     features = names(oligo_markers),
                     cols     = c("#f7fbff", "#084594")) +
  RotatedAxis() +
  ggtitle("Marker expression across oligodendrocyte sub-types") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12))

pdf("OligoOnly_DotPlot_Markers.pdf", width = 7, height = 5)
print(p_dotplot)
dev.off()

cat("All oligo UMAP plots (including full combined patchwork) saved.\n")

#### OLIGO UMAP v3 ------
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)

# ---- Blue gradient palette for marker expression 
pal_blues <- colorRampPalette(c(
  "#f7fbff", "#deebf7", "#c6dbef", "#9ecae1",
  "#6baed6", "#4292c6", "#2171b5", "#2166ac", "#084594"
))(100)

oligo_clusters <- c("Oligodendrocyte precursor cells",
                    "Newly formed oligodendrocytes",
                    "Perineuronal oligodendrocytes",
                    "Myelinating oligodendrocytes")

oligo_only <- subset(pv_cortical_filtered_pnn, idents = oligo_clusters)

oligo_only <- NormalizeData(oligo_only)
oligo_only <- FindVariableFeatures(oligo_only, nfeatures = 2000)
oligo_only <- ScaleData(oligo_only)
oligo_only <- RunPCA(oligo_only, npcs = 15)
oligo_only <- FindNeighbors(oligo_only, dims = 1:10)
oligo_only <- RunUMAP(oligo_only, dims = 1:10)

oligo_only$Genotype_short <- recode(oligo_only$Genotype,
                                    "PV-Cre/tdTom"              = "Ctrl",
                                    "PV-Cre/tdTom/Dnmt1 loxP2" = "KO"
)

oligo_colors <- c(
  "Oligodendrocyte precursor cells" = "#4BACC6",
  "Perineuronal oligodendrocytes"   = "#ea8a33",
  "Newly formed oligodendrocytes"   = "#E8A0BF",
  "Myelinating oligodendrocytes"    = "#9B2335"
)

short_names <- c(
  "Oligodendrocyte precursor cells" = "OPC",
  "Perineuronal oligodendrocytes"   = "PNOL",
  "Newly formed oligodendrocytes"   = "NFOL",
  "Myelinating oligodendrocytes"    = "MOL"
)


# 1) Overview panels: cluster-colored UMAP + genotype UMAP


p_oligo_clusters <- DimPlot(oligo_only, reduction = "umap",
                            label = FALSE, repel = TRUE, pt.size  = 0.3) +
  scale_colour_manual(values = oligo_colors) +
  ggtitle("Oligodendrocyte lineage sub-types") +
  theme_classic() +
  theme(
    plot.title      = element_text(hjust = 0.5, face = "bold", size = 11),
    legend.position = "right",
    axis.line       = element_line(color = "black", linewidth = 0.4)
  )

p_oligo_genotype <- DimPlot(oligo_only, reduction = "umap",
                            group.by = "Genotype_short",
                            cols     = c("Ctrl" = "#c6c6c6", "KO" = "#7FCDCD"),
                            pt.size  = 0.3) +
  ggtitle("Oligodendrocyte UMAP: Ctrl vs KO") +
  theme_classic() +
  theme(
    plot.title      = element_text(hjust = 0.5, face = "bold", size = 11),
    legend.position = "right",
    axis.line       = element_line(color = "black", linewidth = 0.4)
  )


# 2) Per-cluster WT/KO panels


cluster_plot_list <- list()

for (cl in oligo_clusters) {
  
  oligo_only$plot_group <- ifelse(
    as.character(Idents(oligo_only)) == cl,
    as.character(oligo_only$Genotype_short),
    "Other"
  )
  
  oligo_only$plot_group <- factor(oligo_only$plot_group,
                                  levels = c("Other", "Ctrl", "KO"))
  
  p <- DimPlot(oligo_only,
               reduction = "umap",
               group.by  = "plot_group",
               order     = c("KO", "Ctrl", "Other"),
               pt.size   = 0.2,
               alpha     = 1) +
    scale_colour_manual(
      values = c("Other" = "grey96",
                 "Ctrl"  = "#c6c6c6",
                 "KO"    = "#7FCDCD"),
      labels = c("Other" = "Other clusters",
                 "Ctrl"  = "Ctrl",
                 "KO"    = "KO")
    ) +
    ggtitle(cl) +
    theme_classic() +
    theme(
      plot.title      = element_text(hjust = 0.5, face = "bold", size = 11),
      legend.position = "right",
      axis.line       = element_line(color = "black", linewidth = 0.4)
    ) +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))
  
  cluster_plot_list[[cl]] <- p
}

p_combined_percluster <- (
  cluster_plot_list[[1]] + cluster_plot_list[[2]] +
    cluster_plot_list[[3]] + cluster_plot_list[[4]]
) +
  plot_layout(ncol = 2) +
  plot_annotation(
    title = "Oligodendrocyte lineage: Ctrl vs KO per sub-type",
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
  )

pdf("OligoOnly_UMAP_percluster_WT_KO.pdf", width = 11, height = 8)
print(p_combined_percluster)
dev.off()


# 3) Expanded marker gene sets (testing additional candidates)


marker_sets <- list(
  OPC  = c("Pdgfra"),
  NFOL = c("Enpp6", "Bcas1", "Tcf7l2", "Gpr17", "Itpr2"),
  PNOL = c("Pdgfc"),
  MOL  = c("Mbp", "Plp1", "Mog", "Mag", "Mobp", "Cnp")
)

genes_in_obj <- rownames(oligo_only)
marker_sets  <- lapply(marker_sets, function(g) intersect(g, genes_in_obj))

DefaultAssay(oligo_only) <- "RNA"

marker_plot_list <- list()

for (stage in names(marker_sets)) {
  genes <- marker_sets[[stage]]
  if (length(genes) == 0) next
  
  for (gene in genes) {
    p <- FeaturePlot(oligo_only,
                     features  = gene,
                     reduction = "umap",
                     pt.size   = 0.3,
                     order     = TRUE) +
      scale_colour_gradientn(colours = pal_blues) +
      ggtitle(gene) +
      theme_classic() +
      theme(
        plot.title      = element_text(hjust = 0.5, face = "bold", size = 11),
        legend.position = "right",
        axis.line       = element_line(color = "black", linewidth = 0.4)
      )
    
    marker_plot_list[[gene]] <- p
  }
}


# 4) Combined patchwork: overview + per-cluster + ALL marker candidates


n_marker_panels <- length(marker_plot_list)
n_total_panels  <- 6 + n_marker_panels   # 2 overview + 4 per-cluster + markers

p_combined_all <- wrap_plots(
  c(list(p_oligo_clusters, p_oligo_genotype),
    cluster_plot_list,
    marker_plot_list),
  ncol = 2
) +
  plot_annotation(
    title = "Oligodendrocyte lineage: clusters, genotype, and candidate marker expression",
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))
  )

pdf("OligoOnly_UMAP_FULL_combined.pdf",
    width = 11, height = 4.2 * ceiling(n_total_panels / 2))
print(p_combined_all)
dev.off()


# 5) DotPlot summary across ALL candidate markers (blue gradient)


dot_genes <- unlist(marker_sets, use.names = FALSE)

p_dotplot <- DotPlot(oligo_only,
                     features = dot_genes,
                     cols     = c("#f7fbff", "#084594")) +
  RotatedAxis() +
  ggtitle("Candidate marker expression across oligodendrocyte sub-types") +
  theme(
    plot.title  = element_text(hjust = 0.5, face = "bold", size = 12),
    axis.text.x = element_text(face = "italic")
  )

pdf("OligoOnly_DotPlot_Markers.pdf", width = 10, height = 5)
print(p_dotplot)
dev.off()

cat("All combined oligo UMAP + marker patchworks saved (no individual panel PDFs).\n")
#### OLIGO MARKERS ----
library(Seurat)
library(dplyr)
library(tidyr)

DefaultAssay(oligo_only) <- "RNA"
Idents(oligo_only) <- factor(Idents(oligo_only), levels = oligo_clusters)

all_markers <- FindAllMarkers(
  oligo_only,
  only.pos = FALSE,
  min.pct = 0.1,
  logfc.threshold = 0.1,
  test.use = "wilcox"
)

all_markers <- all_markers %>%
  mutate(direction = ifelse(avg_log2FC > 0, "up", "down"))

top50_per_cluster <- all_markers %>%
  group_by(cluster, direction) %>%
  arrange(p_val_adj, desc(abs(avg_log2FC))) %>%
  slice_head(n = 50) %>%
  ungroup() %>%
  select(cluster, direction, gene, avg_log2FC, pct.1, pct.2, p_val, p_val_adj)

write.csv(top50_per_cluster, "OligoOnly_Top50_Markers_PerCluster.csv", row.names = FALSE)

pairwise_results <- list()

cluster_pairs <- combn(oligo_clusters, 2, simplify = FALSE)

for (pair in cluster_pairs) {
  cl1 <- pair[1]
  cl2 <- pair[2]
  
  de <- FindMarkers(
    oligo_only,
    ident.1 = cl1,
    ident.2 = cl2,
    min.pct = 0.1,
    logfc.threshold = 0.1,
    test.use = "wilcox"
  )
  
  de$gene <- rownames(de)
  de$comparison <- paste0(short_names[cl1], "_vs_", short_names[cl2])
  de$direction <- ifelse(de$avg_log2FC > 0,
                         paste0("up_in_", short_names[cl1]),
                         paste0("up_in_", short_names[cl2]))
  
  top_up <- de %>% arrange(p_val_adj, desc(avg_log2FC)) %>% slice_head(n = 50)
  top_down <- de %>% arrange(p_val_adj, avg_log2FC) %>% slice_head(n = 50)
  
  pairwise_results[[paste0(short_names[cl1], "_vs_", short_names[cl2])]] <-
    bind_rows(top_up, top_down) %>%
    select(comparison, direction, gene, avg_log2FC, pct.1, pct.2, p_val, p_val_adj)
}

pairwise_combined <- bind_rows(pairwise_results)

write.csv(pairwise_combined, "OligoOnly_Top50_Markers_Pairwise.csv", row.names = FALSE)

cat("Saved: OligoOnly_Top50_Markers_PerCluster.csv and OligoOnly_Top50_Markers_Pairwise.csv\n")
#### PCT check ----
DefaultAssay(oligo_only) <- "RNA"
Idents(oligo_only) <- factor(Idents(oligo_only), levels = oligo_clusters)

genes_to_check <- c("Pdgfc", "Enpp6")

genes_to_check <- intersect(genes_to_check, rownames(oligo_only))
missing_genes <- setdiff(c("Pdgfc", "Enpp6"), genes_to_check)
if (length(missing_genes) > 0) {
  warning("Not found in object: ", paste(missing_genes, collapse = ", "))
}

pct_results <- list()

for (gene in genes_to_check) {
  for (cl in oligo_clusters) {
    res <- FindMarkers(
      oligo_only,
      ident.1 = cl,
      features = gene,
      min.pct = 0,
      logfc.threshold = 0,
      test.use = "wilcox"
    )
    res$gene <- gene
    res$cluster <- cl
    pct_results[[paste0(gene, "_", cl)]] <- res
  }
}

pct_summary <- bind_rows(pct_results) %>%
  select(gene, cluster, pct.1, pct.2, avg_log2FC, p_val_adj) %>%
  arrange(gene, cluster)

write.csv(pct_summary, "OligoOnly_PctCheck_Pdgfc_Enpp6.csv", row.names = FALSE)

print(pct_summary)
#### OLIGO MARKERS v2 one-vs-all -----
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)

pal_blues <- colorRampPalette(c(
  "#f7fbff", "#deebf7", "#c6dbef", "#9ecae1",
  "#6baed6", "#4292c6", "#2171b5", "#2166ac", "#084594"
))(100)

DefaultAssay(oligo_only) <- "RNA"
Idents(oligo_only) <- factor(Idents(oligo_only), levels = oligo_clusters)

oligo_colors <- c(
  "Oligodendrocyte precursor cells" = "#4BACC6",
  "Perineuronal oligodendrocytes"   = "#ea8a33",
  "Newly formed oligodendrocytes"   = "#E8A0BF",
  "Myelinating oligodendrocytes"    = "#9B2335"
)

# Panel A: four-cluster overview UMAP
p_A <- DimPlot(oligo_only, reduction = "umap", pt.size = 0.3) +
  scale_colour_manual(values = oligo_colors) +
  ggtitle("Oligodendrocyte lineage sub-types") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
    legend.position = "right",
    axis.line = element_line(color = "black", linewidth = 0.4)
  )

# Panel B/C: Enpp6 and Pdgfc FeaturePlots (matched color scale)
feat_genes <- c("Enpp6", "Pdgfc")
p_feat_list <- lapply(feat_genes, function(g) {
  FeaturePlot(oligo_only, features = g, reduction = "umap",
              pt.size = 0.3, order = TRUE) +
    scale_colour_gradientn(colours = pal_blues, limits = c(0, NA)) +
    ggtitle(g) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 11, face2 = "italic"),
      legend.position = "right",
      axis.line = element_line(color = "black", linewidth = 0.4)
    )
})
p_BC <- wrap_plots(p_feat_list, ncol = 2)

# Panel D: PNOL-specific DEGs, EXCLUDING Pdgfc (the selection marker)
pnol_markers <- FindMarkers(
  oligo_only,
  ident.1 = "Perineuronal oligodendrocytes",
  min.pct = 0.1,
  logfc.threshold = 0.25,
  test.use = "wilcox"
)

pnol_markers$gene <- rownames(pnol_markers)

pnol_markers_filtered <- pnol_markers %>%
  filter(gene != "Pdgfc") %>%
  filter(avg_log2FC > 0) %>%
  arrange(p_val_adj, desc(avg_log2FC)) %>%
  slice_head(n = 15)

write.csv(pnol_markers_filtered, "PNOL_DEGs_excluding_Pdgfc.csv", row.names = FALSE)

p_D <- DotPlot(
  oligo_only,
  features = pnol_markers_filtered$gene,
  cols = c("#f7fbff", "#084594")
) +
  RotatedAxis() +
  ggtitle("Top PNOL-enriched genes (excluding selection marker Pdgfc)") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
    axis.text.x = element_text(face = "italic", size = 9)
  )

# Combined Ext Data Fig 8
p_extdata8 <- (p_A | p_BC) /
  p_D +
  plot_layout(heights = c(1, 1.2)) +
  plot_annotation(
    title = "Extended Data Fig. 8: Oligodendrocyte lineage sub-type validation",
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
  )

pdf("ExtDataFig8_OligoLineage_Validation.pdf", width = 12, height = 12)
print(p_extdata8)
dev.off()

# Also save the pct.1/pct.2 table for Pdgfc/Enpp6 for the legend/rebuttal text
pct_check <- data.frame()
for (gene in feat_genes) {
  for (cl in oligo_clusters) {
    res <- FindMarkers(oligo_only, ident.1 = cl, features = gene,
                       min.pct = 0, logfc.threshold = 0, test.use = "wilcox")
    res$gene <- gene
    res$cluster <- cl
    pct_check <- bind_rows(pct_check, res)
  }
}
write.csv(pct_check, "ExtDataFig8_PctCheck_table.csv", row.names = FALSE)

cat("Saved: ExtDataFig8_OligoLineage_Validation.pdf, PNOL_DEGs_excluding_Pdgfc.csv, ExtDataFig8_PctCheck_table.csv\n")
#### OLIGO MARKERS v3 -----
library(Seurat)
library(dplyr)

DefaultAssay(oligo_only) <- "RNA"
Idents(oligo_only) <- factor(Idents(oligo_only), levels = oligo_clusters)

pnol_vs_nfol <- FindMarkers(
  oligo_only,
  ident.1 = "Perineuronal oligodendrocytes",
  ident.2 = "Newly formed oligodendrocytes",
  min.pct = 0.1,
  logfc.threshold = 0.1,
  test.use = "wilcox"
)
pnol_vs_nfol$gene <- rownames(pnol_vs_nfol)
pnol_vs_nfol$comparison <- "PNOL_vs_NFOL"
pnol_vs_nfol$direction <- ifelse(pnol_vs_nfol$avg_log2FC > 0, "up_in_PNOL", "up_in_NFOL")

pnol_vs_mol <- FindMarkers(
  oligo_only,
  ident.1 = "Perineuronal oligodendrocytes",
  ident.2 = "Myelinating oligodendrocytes",
  min.pct = 0.1,
  logfc.threshold = 0.1,
  test.use = "wilcox"
)
pnol_vs_mol$gene <- rownames(pnol_vs_mol)
pnol_vs_mol$comparison <- "PNOL_vs_MOL"
pnol_vs_mol$direction <- ifelse(pnol_vs_mol$avg_log2FC > 0, "up_in_PNOL", "up_in_MOL")

pnol_vs_nfol_top100 <- pnol_vs_nfol %>%
  filter(gene != "Pdgfc") %>%
  arrange(p_val_adj, desc(abs(avg_log2FC))) %>%
  slice_head(n = 100) %>%
  select(comparison, direction, gene, avg_log2FC, pct.1, pct.2, p_val, p_val_adj)

pnol_vs_mol_top100 <- pnol_vs_mol %>%
  filter(gene != "Pdgfc") %>%
  arrange(p_val_adj, desc(abs(avg_log2FC))) %>%
  slice_head(n = 100) %>%
  select(comparison, direction, gene, avg_log2FC, pct.1, pct.2, p_val, p_val_adj)

write.csv(pnol_vs_nfol_top100, "PNOL_vs_NFOL_Top100.csv", row.names = FALSE)
write.csv(pnol_vs_mol_top100, "PNOL_vs_MOL_Top100.csv", row.names = FALSE)

pnol_specific_both <- intersect(
  pnol_vs_nfol_top100 %>% filter(direction == "up_in_PNOL") %>% pull(gene),
  pnol_vs_mol_top100 %>% filter(direction == "up_in_PNOL") %>% pull(gene)
)

cat("Genes upregulated in PNOL vs BOTH neighbors (top 100 each):\n")
print(pnol_specific_both)

pnol_specific_df <- data.frame(gene = pnol_specific_both)
write.csv(pnol_specific_df, "PNOL_genes_specific_vs_both_neighbors_top100.csv", row.names = FALSE)

cat("Saved: PNOL_vs_NFOL_Top100.csv, PNOL_vs_MOL_Top100.csv, PNOL_genes_specific_vs_both_neighbors_top100.csv\n")