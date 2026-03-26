##### Libraries ---------
# Add at the top of your script, before any plotting:
library(showtext)
showtext_auto()
font_add("Arial", "/System/Library/Fonts/Supplemental/Arial.ttf")
# Load the library
library(extrafont)


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

##### Quick reorder ------
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

#### Plot UMAP with defined colors --------
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

DimPlot(pv_cortical_filtered_pnn, reduction = "umap",repel = TRUE) + 
  scale_colour_manual(values = colors) +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 3)))

umap_2 <- DimPlot(pv_cortical_filtered_pnn, reduction = "umap",repel = TRUE) + 
  scale_colour_manual(values = colors) +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 3)))

png(filename = "umap_pnn.png", width = 10, height = 7, units = "in", res = 600)
umap_2
dev.off()

##### Find marker genes ----------------
# Find marker genes
all_markers_cortical_pnn <- FindAllMarkers(pv_cortical_filtered_pnn,
                                           min.pct = 0.25,          # Gene must be detected in ≥25% cells
                                           logfc.threshold = 0.25,  # Minimum fold change threshold
                                           only.pos = FALSE,         # Keep both upregulated and downregulated markers
                                           test.use = "wilcox"      # Default Wilcoxon Rank Sum test
)

#### Visualize marker genes -------
# Dot blot with top 5 marker genes
top5_markers_pnn <- all_markers_cortical_pnn %>%
  group_by(cluster) %>%
  filter(p_val_adj < 0.05) %>%  # Adjusted p-value cutoff
  filter(!grepl("^ENSMUSG", gene)) %>% 
  slice_max(avg_log2FC, n = 5) %>%
  pull(gene) %>%                        # Extract gene names
  unique() 


p <- DotPlot(pv_cortical_filtered_pnn, 
             features = top5_markers_pnn,
             dot.scale = 3,
             col.min = -1.5,                       # Set color scale limits
             col.max = 2.5
)

# With lncRNAs
png("top5_markers_with_lncRNAs_pnn.png", width = 30, height = 10, units = "in", res = 300, type = "quartz")
print(p + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, face = "italic")))
dev.off()

# Dot blot with top 5 marker genes without Gms
top5_markers_nogm_norik_pnn <- all_markers_cortical_pnn %>%
  group_by(cluster) %>%
  filter(p_val_adj < 0.05) %>%  # Adjusted p-value cutoff
  filter(!grepl("^ENSMUSG", gene)) %>% 
  filter(!grepl("^Gm", gene)) %>% 
  filter(!grepl("Rik$", gene)) %>%
  slice_max(avg_log2FC, n = 5) %>%
  pull(gene) %>%                        # Extract gene names
  unique() 


p <- DotPlot(pv_cortical_filtered_pnn, 
             features = top5_markers_nogm_norik_pnn,
             dot.scale = 3,
             col.min = -1.5,                       # Set color scale limits
             col.max = 2.5
)


# Without Gm/Rik genes
png("top5_markers_without_gms_riks_pnn_45deg.png", width = 30, height = 10, units = "in", res = 300, type = "quartz")
print(p + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, face = "italic")))
dev.off()


##### Differential expression analysis (only for NFOL and PNOL) -------
# Load required libraries
library(ggplot2)
library(ggrepel)

# Parameters
wt_label  <- "PV-Cre/tdTom"
ko_label  <- "PV-Cre/tdTom/Dnmt1 loxP2"
logfc_thr <- 0.26
padj_thr  <- 0.05

# Check metadata
table(pv_cortical_filtered_pnn$Genotype)

# Only these two clusters
clusters <- c("Perineuronal oligodendrocytes", "Newly formed oligodendrocytes")

# Helper: sanitize cluster id for filenames
sanitize_id <- function(x) gsub("[^A-Za-z0-9_]", "_", x)

# Helper: build volcano plot with different labeling modes
build_volcano <- function(de, clust, label_mode = c("top10", "signif", "none")) {
  label_mode <- match.arg(label_mode)
  
  # Base label column
  de$label <- ""
  
  # Determine which genes to label
  if (label_mode == "top10") {
    # Masks for significant UP and DOWN based on thresholds
    is_up   <- (de$avg_log2FC >= logfc_thr) & (de$p_val_adj < padj_thr)
    is_down <- (de$avg_log2FC <= -logfc_thr) & (de$p_val_adj < padj_thr)
    
    # Safety: enforce disjoint sets
    if (any(is_up & is_down, na.rm = TRUE)) {
      warning("Rows classified as both UP and DOWN; check thresholds/data.")
      is_down[is_up] <- FALSE
    }
    
    # Top-10 by adjusted p-value within each side (de is pre-sorted by p_val_adj)
    up_df   <- de[is_up, , drop = FALSE]
    down_df <- de[is_down, , drop = FALSE]
    
    if (nrow(up_df) > 0) {
      up_df <- head(up_df, 10)
      de$label[match(up_df$gene, de$gene)] <- up_df$gene
    }
    if (nrow(down_df) > 0) {
      down_df <- head(down_df, 10)
      de$label[match(down_df$gene, de$gene)] <- down_df$gene
    }
  } else if (label_mode == "signif") {
    de$label[de$diffexpressed != "NO"] <- de$gene
  } # "none" leaves labels empty
  
  # Base plot
  p <- ggplot(
    de,
    aes(x = avg_log2FC, y = neg_log10_padj, color = diffexpressed)
  ) +
    geom_point(alpha = 0.6, size = 1.5, show.legend = TRUE) +
    scale_color_manual(
      values = c("UP" = "#b2182b", "DOWN" = "#2166ac", "NO" = "grey70"),
      labels = c("DOWN" = "Downregulated", "NO" = "Not significant", "UP" = "Upregulated")
    ) +
    geom_vline(xintercept = c(-logfc_thr, logfc_thr), linetype = "dashed", color = "grey60", linewidth = 0.2) +
    geom_hline(yintercept = -log10(padj_thr), linetype = "dashed", color = "grey60", linewidth = 0.2) +
    labs(
      title = paste0("Cluster ", clust, ": KO vs WT"),
      x = "Log2 Fold Change",
      y = "-Log10(Adjusted P-value)",
      color = "Differential Expression"
    ) +
    theme_classic() +
    theme(
      plot.title  = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.title  = element_text(size = 12),
      axis.text   = element_text(size = 10),
      legend.position = "right",
      axis.line   = element_line(color = "black", linewidth = 0.5)
    ) +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))
  
  # Add labels depending on mode
  if (label_mode == "top10" && any(nzchar(de$label))) {
    set.seed(123)
    p <- p +
      ggrepel::geom_text_repel(
        aes(label = label),
        max.overlaps = 20,
        size = 3,
        box.padding = 0.5,
        point.padding = 0.3,
        show.legend = FALSE,
        segment.color = "grey30",
        segment.size  = 0.3,
        min.segment.length = 0
      )
  } else if (label_mode == "signif" && any(nzchar(de$label))) {
    set.seed(123)
    p <- p +
      ggrepel::geom_text_repel(
        aes(label = label),
        max.overlaps = Inf,
        size = 2.5,
        box.padding = 0.1,
        point.padding = 0.1,
        show.legend = FALSE,
        segment.color = NA,
        segment.size  = 0
      )
  }
  p
}

# Iterate specified clusters only, save in working directory
for (clust in clusters) {
  cells  <- WhichCells(pv_cortical_filtered_pnn, idents = clust)
  if (length(cells) == 0) {
    warning("No cells found for cluster: ", clust)
    next
  }
  subobj <- subset(pv_cortical_filtered_pnn, cells = cells)
  Idents(subobj) <- "Genotype"
  
  de <- FindMarkers(
    object = subobj,
    ident.1 = ko_label,
    ident.2 = wt_label,
    min.pct = 0.1,
    logfc.threshold = logfc_thr
  )
  
  # Ensure gene names
  de$gene <- rownames(de)
  
  # Handle possible older Seurat column name
  if (!"avg_log2FC" %in% colnames(de) && "avg_logFC" %in% colnames(de)) {
    de$avg_log2FC <- de$avg_logFC
  }
  
  # Remove Ensembl-like genes starting with ENSMUSG
  de <- de[!grepl("^ENSMUSG", de$gene, ignore.case = FALSE), , drop = FALSE]
  
  # Sort by adjusted p-value
  if (nrow(de) > 0) {
    de <- de[order(de$p_val_adj), , drop = FALSE]
  }
  
  # Classification for colors (significance)
  de$diffexpressed <- "NO"
  de$diffexpressed[de$avg_log2FC >= logfc_thr & de$p_val_adj < padj_thr] <- "UP"
  de$diffexpressed[de$avg_log2FC <= -logfc_thr & de$p_val_adj < padj_thr] <- "DOWN"
  de$diffexpressed <- factor(de$diffexpressed, levels = c("UP", "DOWN", "NO"))
  
  # y-axis: -log10 padj, with safe handling for zeros
  de$neg_log10_padj <- -log10(de$p_val_adj)
  finite_vals <- de$neg_log10_padj[is.finite(de$neg_log10_padj)]
  if (length(finite_vals) == 0) finite_vals <- 0
  max_finite <- max(finite_vals)
  de$neg_log10_padj[!is.finite(de$neg_log10_padj)] <- max_finite + 10
  
  # Filenames (save into working directory)
  safe_clust <- sanitize_id(clust)
  filename_csv <- paste0("cluster_", safe_clust, "_KO_vs_WT_all_genes_filtered.csv")
  filename_top10 <- paste0("cluster_", safe_clust, "_KO_vs_WT_volcano_top10_filtered.png")
  filename_sig   <- paste0("cluster_", safe_clust, "_KO_vs_WT_volcano_sig_labeled_filtered.pdf")
  filename_clean <- paste0("cluster_", safe_clust, "_KO_vs_WT_volcano_clean_filtered.png")
  
  # Write CSV
  write.csv(de, file = filename_csv, row.names = FALSE)
  message("Saved filtered DE results for cluster ", clust, " to ", filename_csv)
  
  # Build and save three plot variants (saved in working directory)
  p_top10 <- build_volcano(de, clust, label_mode = "top10")
  ggsave(filename = filename_top10, plot = p_top10, width = 10, height = 8, dpi = 300, units = "in", device = "png")
  
  p_sig <- build_volcano(de, clust, label_mode = "signif")
  ggsave(filename = filename_sig, plot = p_sig, width = 10, height = 8, units = "in", device = "pdf")
  
  p_clean <- build_volcano(de, clust, label_mode = "none")
  ggsave(filename = filename_clean, plot = p_clean, width = 10, height = 8, dpi = 300, units = "in", device = "png")
}


##### Subset by Genotype -----------------------------------------------------------
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

# 3. Calculate proportion of each cluster in WT and KO
wt_cortical_cluster_counts_pnn$Proportion_WT <- wt_cortical_cluster_counts_pnn$CellNumber_WT / total_wt_cortical_pnn
ko_cortical_cluster_counts_pnn$Proportion_KO <- ko_cortical_cluster_counts_pnn$CellNumber_KO / total_ko_cortical_pnn

# 4. Merge the tables
merged_cortical_counts_pnn <- merge(wt_cortical_cluster_counts_pnn, ko_cortical_cluster_counts_pnn, by = "Cluster", all = TRUE)
merged_cortical_counts_pnn[is.na(merged_cortical_counts_pnn)] <- 0  # Replace NAs with 0

# 5. Calculate relative difference in proportions (as percent)
merged_cortical_counts_pnn$RelativeDifference <- (merged_cortical_counts_pnn$Proportion_KO - merged_cortical_counts_pnn$Proportion_WT) / merged_cortical_counts_pnn$Proportion_WT * 100

# 6. Handle Inf/NaN (if Proportion_WT is zero)
merged_cortical_counts_pnn$RelativeDifference[!is.finite(merged_cortical_counts_pnn$RelativeDifference)] <- NA

merged_cortical_counts_pnn$Proportion_WT <- round(merged_cortical_counts_pnn$Proportion_WT * 100, 2)
merged_cortical_counts_pnn$Proportion_KO <- round(merged_cortical_counts_pnn$Proportion_KO * 100, 2)
merged_cortical_counts_pnn$RelativeDifference <- round(merged_cortical_counts_pnn$RelativeDifference, 2)

# Assuming merged_counts from your previous code
merged_cortical_counts_pnn$p_value <- NA

for (i in 1:nrow(merged_cortical_counts_pnn)) {
  k_WT_pnn <- merged_cortical_counts_pnn$CellNumber_WT[i]
  k_KO_pnn <- merged_cortical_counts_pnn$CellNumber_KO[i]
  # Only test if both groups have at least one cell
  if ((k_WT_pnn + k_KO_pnn) > 0) {
    test <- prop.test(
      x = c(k_WT_pnn, k_KO_pnn),
      n = c(total_wt_cortical_pnn, total_ko_cortical_pnn),
      alternative = "two.sided"
    )
    merged_cortical_counts_pnn$p_value[i] <- test$p.value
  }
}
merged_cortical_counts_pnn$adj_p_value <- p.adjust(merged_cortical_counts_pnn$p_value, method = "fdr")
write.csv(merged_cortical_counts_pnn, file = "cluster_relative_difference_with_pvalues_cortical_pnn.csv", row.names = FALSE)


##### CellChat -----
# For WT

data.input <- wt_cortical_pnn[["RNA"]]$data # normalized data matrix
labels <- Idents(wt_cortical_pnn)
meta <- data.frame(labels = labels, row.names = names(labels))
cellchat_wt_cortical_pnn <- createCellChat(object = wt_cortical_pnn, group.by = "ident", assay = "RNA")


CellChatDB <- CellChatDB.mouse # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)

CellChatDB.use <- CellChatDB
cellchat_wt_cortical_pnn@DB <- CellChatDB.use

cellchat_wt_cortical_pnn <- subsetData(cellchat_wt_cortical_pnn)
cellchat_wt_cortical_pnn <- updateCellChat(cellchat_wt_cortical_pnn)
#future::plan("multisession", workers = 12) # do parallel
#future::plan("sequential")

cellchat_wt_cortical_pnn <- identifyOverExpressedGenes(cellchat_wt_cortical_pnn)
cellchat_wt_cortical_pnn <- identifyOverExpressedInteractions(cellchat_wt_cortical_pnn)
length(cellchat_wt_cortical_pnn@LR$LRsig)

ptm = Sys.time()
execution.time = Sys.time() - ptm
print(as.numeric(execution.time, units = "secs"))

cellchat_wt_cortical_pnn <- computeCommunProb(cellchat_wt_cortical_pnn, type = "triMean")
cellchat_wt_cortical_pnn <- filterCommunication(cellchat_wt_cortical_pnn, min.cells = 10)

cellchat_wt_cortical_pnn <- computeCommunProbPathway(cellchat_wt_cortical_pnn)

cellchat_wt_cortical_pnn <- aggregateNet(cellchat_wt_cortical_pnn)
execution.time = Sys.time() - ptm
print(as.numeric(execution.time, units = "secs"))

# For KO
data.input <- ko_cortical_pnn[["RNA"]]$data # normalized data matrix
labels <- Idents(ko_cortical_pnn)
meta <- data.frame(labels = labels, row.names = names(labels))
cellchat_ko_cortical_pnn <- createCellChat(object = ko_cortical_pnn, group.by = "ident", assay = "RNA")


CellChatDB <- CellChatDB.mouse # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)

CellChatDB.use <- CellChatDB
cellchat_ko_cortical_pnn@DB <- CellChatDB.use

cellchat_ko_cortical_pnn <- subsetData(cellchat_ko_cortical_pnn)
cellchat_ko_cortical_pnn <- updateCellChat(cellchat_ko_cortical_pnn)
#future::plan("multisession", workers = 12) # do parallel
#future::plan("sequential")

cellchat_ko_cortical_pnn <- identifyOverExpressedGenes(cellchat_ko_cortical_pnn)
cellchat_ko_cortical_pnn <- identifyOverExpressedInteractions(cellchat_ko_cortical_pnn)
length(cellchat_ko_cortical_pnn@LR$LRsig)

ptm = Sys.time()
execution.time = Sys.time() - ptm
print(as.numeric(execution.time, units = "secs"))

cellchat_ko_cortical_pnn <- computeCommunProb(cellchat_ko_cortical_pnn, type = "triMean")
cellchat_ko_cortical_pnn <- filterCommunication(cellchat_ko_cortical_pnn, min.cells = 10)

cellchat_ko_cortical_pnn <- computeCommunProbPathway(cellchat_ko_cortical_pnn)

cellchat_ko_cortical_pnn <- aggregateNet(cellchat_ko_cortical_pnn)
execution.time = Sys.time() - ptm
print(as.numeric(execution.time, units = "secs"))




##### Rename clusters ------
old_cellchat_names_pnn <- c(
  "Layer 2/3 IT neurons",
  "Layer 4 sensory neurons",
  "Layer 5a IT neurons",
  "Layer 5b PT neurons",
  "Layer 5/6 IT neurons",
  "Layer 6a corticothalamic neurons",
  "Layer 6b neurons",
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

new_names_pnn <- c(
  "L2/3 IT",
  "L4 Sensory",
  "L5a IT",
  "L5b PT",
  "L5/6 IT",
  "L6a CT",
  "L6b",
  "Deep ET",
  "CSN Type I",
  "CSN Type II",
  "PV+ Int",
  "SST+ Int",
  "VIP+ Int",
  "Astrocyte",
  "OPC",
  "NFOL",
  "PNOL",
  "MOL",
  "Microglia",
  "Endothelial",
  "Leptomeningeal FB",
  "Meningeal FB"
)

cellchat_wt_cortical_pnn <- updateClusterLabels(
  cellchat_wt_cortical_pnn,
  old.cluster.name = old_cellchat_names_pnn,
  new.cluster.name = new_names_pnn
)

cellchat_ko_cortical_pnn <- updateClusterLabels(
  cellchat_ko_cortical_pnn,
  old.cluster.name = old_cellchat_names_pnn,
  new.cluster.name = new_names_pnn
)

# Verify
levels(cellchat_wt_cortical_pnn@idents)
levels(cellchat_ko_cortical_pnn@idents)

##### Save CellChat object and re-load ------

saveRDS(cellchat_wt_cortical_pnn, file = "cellchat_wt_cortical_pnn_rev.rds")
saveRDS(cellchat_ko_cortical_pnn, file = "cellchat_ko_cortical_pnn_rev.rds")

cellchat_wt_cortical_pnn <- readRDS(file = "cellchat_wt_cortical_pnn_rev.rds")
cellchat_ko_cortical_pnn <- readRDS(file = "cellchat_ko_cortical_pnn_rev.rds")

##### Remove unwanted clusters from CellChat -------
cellchat_wt_cortical_pnn <- subsetCellChat(cellchat_wt_cortical_pnn, idents.use = c("CSN Type I", "CSN Type II", "Endothelial", "Leptomeningeal FB", "Meningeal FB"), invert = TRUE)
cellchat_ko_cortical_pnn <- subsetCellChat(cellchat_ko_cortical_pnn, idents.use = c("CSN Type I", "CSN Type II", "Endothelial", "Leptomeningeal FB", "Meningeal FB"), invert = TRUE)

#########
##### Merge objects and re-compute centrality --------
# For WT
cellchat_wt_cortical_pnn <- netAnalysis_computeCentrality(cellchat_wt_cortical_pnn, slot.name = "netP")
# For KO
cellchat_ko_cortical_pnn <- netAnalysis_computeCentrality(cellchat_ko_cortical_pnn, slot.name = "netP")

# List your CellChat objects
object.list_cortical_new <- list(WT = cellchat_wt_cortical_pnn, KO = cellchat_ko_cortical_pnn)

# Merge for comparison
cellchat_merged_cortical_pnn <- mergeCellChat(object.list_cortical_new, add.names = names(object.list_cortical_new))

### CHECK LEVELS
levels(cellchat_merged_cortical_pnn@idents$joint)
levels(cellchat_wt_cortical_pnn@idents)
levels(cellchat_ko_cortical_pnn@idents)

##########
### Interaction weight and counts ------------------

# List your CellChat objects
# Your original list of CellChat objects
object.list_cortical_pnn <- list(WT = cellchat_wt_cortical_pnn, KO = cellchat_ko_cortical_pnn)

### Interaction number and weight

par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat_merged_cortical_pnn, label.edge= F, weight.scale = T)
netVisual_diffInteraction(cellchat_merged_cortical_pnn, label.edge= F, weight.scale = T, measure = "weight")

gg1 <- netVisual_heatmap(cellchat_merged_cortical_pnn)
gg2 <- netVisual_heatmap(cellchat_merged_cortical_pnn, measure = "weight")
gg1 + gg2
png(filename = "interaction_weight_cortical_pnn.png", width = 8, height = 6, units = "in", res = 600)
gg2
dev.off()

## Quantify the matrix
# 1. Extract the raw count and weight matrices for each group
nets <- cellchat_merged_cortical_pnn@net       # list of two elements
mat_count_ctrl   <- nets[[1]]$count         # WT counts
mat_count_ko     <- nets[[2]]$count         # KO counts
mat_weight_ctrl  <- nets[[1]]$weight        # WT weights
mat_weight_ko    <- nets[[2]]$weight        # KO weights

# 2. Compute differential matrices (KO minus WT)
diff_count  <- mat_count_ko  - mat_count_ctrl
diff_weight <- mat_weight_ko - mat_weight_ctrl

# 3. (Optional) Assign row/column names if lost
rownames(diff_count)  <- rownames(mat_count_ctrl)
colnames(diff_count)  <- colnames(mat_count_ctrl)
rownames(diff_weight) <- rownames(mat_weight_ctrl)
colnames(diff_weight) <- colnames(mat_weight_ctrl)

# 4. Save to CSV for inspection
write.csv(diff_count,
          "netVisual_diff_counts_KO_vs_WT_pnn.csv",
          row.names = TRUE)
write.csv(diff_weight,
          "netVisual_diff_weights_KO_vs_WT_pnn.csv",
          row.names = TRUE)


### Information flow --------
gg1 <- rankNet(cellchat_merged_cortical_pnn, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat_merged_cortical_pnn, mode = "comparison", stacked = F, do.stat = TRUE)
gg1 + gg2

#### Similarity ------
cellchat_merged_cortical_pnn <- computeNetSimilarityPairwise(cellchat_merged_cortical_pnn, type = "functional")
cellchat_merged_cortical_pnn <- netEmbedding(cellchat_merged_cortical_pnn, type = "functional", umap.method = "uwot")
cellchat_merged_cortical_pnn <- netClustering(cellchat_merged_cortical_pnn, type = "functional", do.parallel = FALSE)
png(filename = "similarity_functional_cortical_pnn.png", width = 6, height = 6, units = "in", res = 600)
rankSimilarity(cellchat_merged_cortical_pnn, type = "functional")
dev.off()

cellchat_merged_cortical_pnn <- computeNetSimilarityPairwise(cellchat_merged_cortical_pnn, type = "structural")
cellchat_merged_cortical_pnn <- netEmbedding(cellchat_merged_cortical_pnn, type = "structural", umap.method = "uwot")
cellchat_merged_cortical_pnn <- netClustering(cellchat_merged_cortical_pnn, type = "structural", do.parallel = FALSE)
png(filename = "similarity_structural_cortical_pnn.png", width = 6, height = 6, units = "in", res = 600)
rankSimilarity(cellchat_merged_cortical_pnn, type = "structural")
dev.off()


##########
#### Direction-independent bubble plot with filter (set to incoming) -------
thresh        <- 0.05
direction     <- "incoming"  # "outgoing" (PV -> oligos) or "incoming" (oligos -> PV)
sources_names <- c("PV+ Int")
cluster_names <- c("OPC",  "NFOL", "PNOL", "MOL")  # fixed x-axis order

# ---- String helpers (head = before first dot; tail = after first dot)
head_token <- function(x) sub("\\..*$", "", x)
tail_token <- function(x) sub("^[^.]+\\.", "", x)
has_dot    <- function(x) grepl("\\.", x)

# Pretty LR formatter:
#  "A_B"     -> "A + B"
#  "A_B_C"   -> "A + (B + C)"
#  length>=3 generalized to first + (rest joined by +)
lrp_pretty_plus <- function(full_label) {
  # Use the substring after the first dot if present
  tail <- sub("^[^.]+\\.", "", full_label)
  if (identical(tail, full_label)) tail <- full_label
  parts <- strsplit(tail, "_", fixed = TRUE)[[1]]
  if (length(parts) <= 1L || anyNA(parts)) {
    tail
  } else if (length(parts) == 2L) {
    paste(parts[1], parts[2], sep = " + ")
  } else {
    paste0(parts[1], " + (", paste(parts[-1], collapse = " + "), ")")
  }
}

# ---- Significance masking (mimics subsetCommunication p-value filter)
mask_prob_by_pval <- function(cellchat_net, thresh = 0.05) {
  prob <- cellchat_net$prob
  pval <- cellchat_net$pval
  stopifnot(identical(dim(prob), dim(pval)))
  bad <- !(pval < thresh) | is.na(pval)
  prob_masked <- prob
  prob_masked[bad] <- 0
  list(prob = prob_masked, dnames = dimnames(prob_masked))
}

# ---- Extract selected sender/receiver blocks across all pathways to long format
extract_long <- function(cellchat_net, send_idx, recv_idx, thresh = 0.05) {
  m <- mask_prob_by_pval(cellchat_net, thresh)
  paths <- m$dnames[[3]]
  purrr::map_dfr(seq_along(paths), function(i) {
    mat <- m$prob[send_idx, recv_idx, i, drop = FALSE]
    tibble::as_tibble(mat, rownames = "Source") |>
      tidyr::pivot_longer(-Source, names_to = "Target", values_to = "prob") |>
      dplyr::mutate(Pathway = paths[i])
  })
}

# ---- Build indices by matching head tokens (sender × receiver × pair tensor)
dimn_WT <- dimnames(cellchat_merged_cortical_pnn@net$WT$prob)
rnames  <- dimn_WT[[1]]          # senders
cnames  <- dimn_WT[[2]]          # receivers
rhead   <- head_token(rnames)    # sender heads (cell-type labels)
chead   <- head_token(cnames)    # receiver heads (cell-type labels)

if (direction == "outgoing") {
  # PV sends to oligos: rows head == PV; cols head in cluster_names
  send_idx <- which(rhead %in% sources_names)
  recv_idx <- which(chead %in% cluster_names)
  x_from   <- "Target"
  x_label  <- "Target Cell Type"
} else {
  # Oligos send to PV: rows head in cluster_names; cols head == PV
  send_idx <- which(rhead %in% cluster_names)
  recv_idx <- which(chead %in% sources_names)
  x_from   <- "Source"
  x_label  <- "Source Cell Type"
}

stopifnot(length(send_idx) > 0, length(recv_idx) > 0)

# ---- Extract and compute Δ = KO − WT
wt_long <- extract_long(cellchat_merged_cortical_pnn@net$WT, send_idx, recv_idx, thresh)
ko_long <- extract_long(cellchat_merged_cortical_pnn@net$KO, send_idx, recv_idx, thresh)

delta_df <- dplyr::left_join(
  wt_long, ko_long,
  by = c("Source","Target","Pathway"),
  suffix = c("_WT","_KO")
) |>
  dplyr::mutate(delta = prob_KO - prob_WT)

# ---- Derive columns for plotting
# Cluster: head from the x side (Target for outgoing, Source for incoming)
# LRP: pretty-printed tail from the same side string (falls back to full Pathway if needed)
delta_df <- delta_df |>
  dplyr::mutate(
    Source_head = head_token(Source),
    Target_head = head_token(Target),
    PairStr = if (x_from == "Target") Target else Source,
    CellType = if (x_from == "Target") Target_head else Source_head,
    LRP = ifelse(has_dot(PairStr),
                 vapply(PairStr, lrp_pretty_plus, character(1)),
                 vapply(Pathway,  lrp_pretty_plus, character(1)))
  )

# ---- Lock x to desired columns (keeps exactly one column per cluster in order)
delta_df <- delta_df |>
  dplyr::filter(CellType %in% cluster_names)

# ---- Biological plausibility filter (direction-aware)
implausible_prefixes <- if (direction == "outgoing") {
  # PV is sender: PV cells are GABAergic, not glutamatergic
  c("Glutamate-Glu")
} else {
  # Oligo is sender: not glutamatergic, not GABAergic, not monoaminergic
  c("Glutamate-Glu", "GABA-A", "GABA-B", "SerotoninDopamin")
}

delta_df <- delta_df |>
  dplyr::filter(!Reduce(`|`, lapply(implausible_prefixes, function(p) startsWith(Pathway, p))))

# Filter APP-SORL1 and LIPA in both directions (intracellular trafficking, not intercellular signaling)
delta_df <- delta_df |>
  dplyr::filter(Pathway != "APP_SORL1")

delta_df <- delta_df |> dplyr::filter(Pathway != "Cholesterol-Cholesterol-LIPA_RORA")

delta_df$CellType <- factor(delta_df$CellType, levels = cluster_names)

# ---- Top 10 per cluster by |Δ|, split by sign
top_inc_per_type <- delta_df |>
  dplyr::filter(delta > 0) |>
  dplyr::group_by(CellType) |>
  dplyr::slice_max(order_by = abs(delta), n = 10, with_ties = FALSE) |>
  dplyr::ungroup()

top_dec_per_type <- delta_df |>
  dplyr::filter(delta < 0) |>
  dplyr::group_by(CellType) |>
  dplyr::slice_max(order_by = abs(delta), n = 10, with_ties = FALSE) |>
  dplyr::ungroup()

# ---- Outline color by sign
add_outline_cols <- function(df) {
  dplyr::mutate(df, outline_col = dplyr::case_when(
    delta > 0 ~ "#b2182b",
    delta < 0 ~ "#2166ac",
    TRUE ~ NA_character_
  ))
}
top_inc_per_type <- add_outline_cols(top_inc_per_type)
top_dec_per_type <- add_outline_cols(top_dec_per_type)

# ---- Palettes
increase_pal <- c("#fff7bc","#fec44f","#fb6a4a","#de2d26","#b2182b")
decrease_pal <- c("#fff7bc","#c3eec6","#9ecae1","#3182bd","#2166ac")

# ---- Plotters (separate panels: no need for ggnewscale inside each)
plot_sign <- function(df, legend_title, x_label, pal, fill_var) {
  ggplot2::ggplot(df, ggplot2::aes(x = CellType, y = LRP)) +
    ggplot2::geom_point(ggplot2::aes(fill = !!fill_var, colour = outline_col),
                        shape = 21, size = 6, stroke = 0.3, alpha = 0.95) +
    ggplot2::scale_fill_gradientn(colors = pal, name = legend_title, na.value = "transparent") +
    ggplot2::scale_color_identity(guide = "none") +
    ggplot2::scale_x_discrete(limits = cluster_names, drop = FALSE) +
    ggplot2::labs(x = x_label, y = "Ligand–Receptor Pair") +
    ggplot2::theme_minimal(base_size = 14, base_family = "Arial") +
    ggplot2::theme(
      text = ggplot2::element_text(color = "black"),
      axis.text = ggplot2::element_text(color = "black"),
      axis.title = ggplot2::element_text(color = "black"),
      legend.text = ggplot2::element_text(color = "black"),
      legend.title = ggplot2::element_text(color = "black"),
      strip.text = ggplot2::element_text(color = "black"),
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1, color = "black"),
      axis.text.y = ggplot2::element_text(size = 10, color = "black")
    )
}

p_top_inc <- if (nrow(top_inc_per_type) > 0) {
  plot_sign(top_inc_per_type, "Δ > 0 (KO − WT)", x_label, increase_pal, rlang::quo(delta))
} else {
  ggplot2::ggplot() + ggplot2::theme_void() + ggplot2::ggtitle("No increases pass filters")
}

p_top_dec <- if (nrow(top_dec_per_type) > 0) {
  plot_sign(top_dec_per_type, "Δ < 0 (KO − WT)", x_label, decrease_pal, rlang::quo(-delta))
} else {
  ggplot2::ggplot() + ggplot2::theme_void() + ggplot2::ggtitle("No decreases pass filters")
}

# ---- Combine panels side-by-side
title_txt <- if (direction == "outgoing") {
  paste0("KO vs WT: Top 10 PV→Oligo (outgoing) per cluster (p < ", thresh, ")")
} else {
  paste0("KO vs WT: Top 10 Oligo→PV (incoming) per cluster (p < ", thresh, ")")
}

(p_top_inc | p_top_dec) +
  patchwork::plot_annotation(
    title    = title_txt,
    subtitle = "Left: increases (red); Right: decreases (blue)"
  )

# ---- Optional save
ggsave("differential_bubble_oligo_to_PV_pnn.pdf", width = 16, height = 8)

# Snapshot incoming results
top_inc_incoming <- top_inc_per_type
top_dec_incoming <- top_dec_per_type


#### Direction-independent bubble plot with filter (set to outgoing) -------
thresh        <- 0.05
direction     <- "outgoing"  # "outgoing" (PV -> oligos) or "incoming" (oligos -> PV)
sources_names <- c("PV+ Int")
cluster_names <- c("OPC",  "NFOL", "PNOL", "MOL")  # fixed x-axis order

# ---- String helpers (head = before first dot; tail = after first dot)
head_token <- function(x) sub("\\..*$", "", x)
tail_token <- function(x) sub("^[^.]+\\.", "", x)
has_dot    <- function(x) grepl("\\.", x)

# Pretty LR formatter:
#  "A_B"     -> "A + B"
#  "A_B_C"   -> "A + (B + C)"
#  length>=3 generalized to first + (rest joined by +)
lrp_pretty_plus <- function(full_label) {
  # Use the substring after the first dot if present
  tail <- sub("^[^.]+\\.", "", full_label)
  if (identical(tail, full_label)) tail <- full_label
  parts <- strsplit(tail, "_", fixed = TRUE)[[1]]
  if (length(parts) <= 1L || anyNA(parts)) {
    tail
  } else if (length(parts) == 2L) {
    paste(parts[1], parts[2], sep = " + ")
  } else {
    paste0(parts[1], " + (", paste(parts[-1], collapse = " + "), ")")
  }
}

# ---- Significance masking (mimics subsetCommunication p-value filter)
mask_prob_by_pval <- function(cellchat_net, thresh = 0.05) {
  prob <- cellchat_net$prob
  pval <- cellchat_net$pval
  stopifnot(identical(dim(prob), dim(pval)))
  bad <- !(pval < thresh) | is.na(pval)
  prob_masked <- prob
  prob_masked[bad] <- 0
  list(prob = prob_masked, dnames = dimnames(prob_masked))
}

# ---- Extract selected sender/receiver blocks across all pathways to long format
extract_long <- function(cellchat_net, send_idx, recv_idx, thresh = 0.05) {
  m <- mask_prob_by_pval(cellchat_net, thresh)
  paths <- m$dnames[[3]]
  purrr::map_dfr(seq_along(paths), function(i) {
    mat <- m$prob[send_idx, recv_idx, i, drop = FALSE]
    tibble::as_tibble(mat, rownames = "Source") |>
      tidyr::pivot_longer(-Source, names_to = "Target", values_to = "prob") |>
      dplyr::mutate(Pathway = paths[i])
  })
}

# ---- Build indices by matching head tokens (sender × receiver × pair tensor)
dimn_WT <- dimnames(cellchat_merged_cortical_pnn@net$WT$prob)
rnames  <- dimn_WT[[1]]          # senders
cnames  <- dimn_WT[[2]]          # receivers
rhead   <- head_token(rnames)    # sender heads (cell-type labels)
chead   <- head_token(cnames)    # receiver heads (cell-type labels)

if (direction == "outgoing") {
  # PV sends to oligos: rows head == PV; cols head in cluster_names
  send_idx <- which(rhead %in% sources_names)
  recv_idx <- which(chead %in% cluster_names)
  x_from   <- "Target"
  x_label  <- "Target Cell Type"
} else {
  # Oligos send to PV: rows head in cluster_names; cols head == PV
  send_idx <- which(rhead %in% cluster_names)
  recv_idx <- which(chead %in% sources_names)
  x_from   <- "Source"
  x_label  <- "Source Cell Type"
}

stopifnot(length(send_idx) > 0, length(recv_idx) > 0)

# ---- Extract and compute Δ = KO − WT
wt_long <- extract_long(cellchat_merged_cortical_pnn@net$WT, send_idx, recv_idx, thresh)
ko_long <- extract_long(cellchat_merged_cortical_pnn@net$KO, send_idx, recv_idx, thresh)


delta_df <- dplyr::left_join(
  wt_long, ko_long,
  by = c("Source","Target","Pathway"),
  suffix = c("_WT","_KO")
) |>
  dplyr::mutate(delta = prob_KO - prob_WT)

# ---- Derive columns for plotting
# Cluster: head from the x side (Target for outgoing, Source for incoming)
# LRP: pretty-printed tail from the same side string (falls back to full Pathway if needed)
delta_df <- delta_df |>
  dplyr::mutate(
    Source_head = head_token(Source),
    Target_head = head_token(Target),
    PairStr = if (x_from == "Target") Target else Source,
    CellType = if (x_from == "Target") Target_head else Source_head,
    LRP = ifelse(has_dot(PairStr),
                 vapply(PairStr, lrp_pretty_plus, character(1)),
                 vapply(Pathway,  lrp_pretty_plus, character(1)))
  )

# ---- Lock x to desired columns (keeps exactly one column per cluster in order)
delta_df <- delta_df |>
  dplyr::filter(CellType %in% cluster_names)

# ---- Biological plausibility filter (direction-aware)
implausible_prefixes <- if (direction == "outgoing") {
  # PV is sender: PV cells are GABAergic, not glutamatergic
  c("Glutamate-Glu")
} else {
  # Oligo is sender: not glutamatergic, not GABAergic, not monoaminergic
  c("Glutamate-Glu", "GABA-A", "GABA-B", "SerotoninDopamin")
}

delta_df <- delta_df |>
  dplyr::filter(!Reduce(`|`, lapply(implausible_prefixes, function(p) startsWith(Pathway, p))))

# Filter APP-SORL1 and LIPA in both directions (intracellular trafficking, not intercellular signaling)
delta_df <- delta_df |>
  dplyr::filter(Pathway != "APP_SORL1")

delta_df <- delta_df |> dplyr::filter(Pathway != "Cholesterol-Cholesterol-LIPA_RORA")

delta_df$CellType <- factor(delta_df$CellType, levels = cluster_names)

# ---- Top 10 per cluster by |Δ|, split by sign
top_inc_per_type <- delta_df |>
  dplyr::filter(delta > 0) |>
  dplyr::group_by(CellType) |>
  dplyr::slice_max(order_by = abs(delta), n = 10, with_ties = FALSE) |>
  dplyr::ungroup()

top_dec_per_type <- delta_df |>
  dplyr::filter(delta < 0) |>
  dplyr::group_by(CellType) |>
  dplyr::slice_max(order_by = abs(delta), n = 10, with_ties = FALSE) |>
  dplyr::ungroup()

# ---- Outline color by sign
add_outline_cols <- function(df) {
  dplyr::mutate(df, outline_col = dplyr::case_when(
    delta > 0 ~ "#b2182b",
    delta < 0 ~ "#2166ac",
    TRUE ~ NA_character_
  ))
}
top_inc_per_type <- add_outline_cols(top_inc_per_type)
top_dec_per_type <- add_outline_cols(top_dec_per_type)

# ---- Palettes
increase_pal <- c("#fff7bc","#fec44f","#fb6a4a","#de2d26","#b2182b")
decrease_pal <- c("#fff7bc","#c3eec6","#9ecae1","#3182bd","#2166ac")

# ---- Plotters (separate panels: no need for ggnewscale inside each)
plot_sign <- function(df, legend_title, x_label, pal, fill_var) {
  ggplot2::ggplot(df, ggplot2::aes(x = CellType, y = LRP)) +
    ggplot2::geom_point(ggplot2::aes(fill = !!fill_var, colour = outline_col),
                        shape = 21, size = 6, stroke = 0.3, alpha = 0.95) +
    ggplot2::scale_fill_gradientn(colors = pal, name = legend_title, na.value = "transparent") +
    ggplot2::scale_color_identity(guide = "none") +
    ggplot2::scale_x_discrete(limits = cluster_names, drop = FALSE) +
    ggplot2::labs(x = x_label, y = "Ligand–Receptor Pair") +
    ggplot2::theme_minimal(base_size = 14, base_family = "Arial") +
    ggplot2::theme(
      text = ggplot2::element_text(color = "black"),
      axis.text = ggplot2::element_text(color = "black"),
      axis.title = ggplot2::element_text(color = "black"),
      legend.text = ggplot2::element_text(color = "black"),
      legend.title = ggplot2::element_text(color = "black"),
      strip.text = ggplot2::element_text(color = "black"),
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1, color = "black"),
      axis.text.y = ggplot2::element_text(size = 10, color = "black")
    )
}

p_top_inc <- if (nrow(top_inc_per_type) > 0) {
  plot_sign(top_inc_per_type, "Δ > 0 (KO − WT)", x_label, increase_pal, rlang::quo(delta))
} else {
  ggplot2::ggplot() + ggplot2::theme_void() + ggplot2::ggtitle("No increases pass filters")
}

p_top_dec <- if (nrow(top_dec_per_type) > 0) {
  plot_sign(top_dec_per_type, "Δ < 0 (KO − WT)", x_label, decrease_pal, rlang::quo(-delta))
} else {
  ggplot2::ggplot() + ggplot2::theme_void() + ggplot2::ggtitle("No decreases pass filters")
}

# ---- Combine panels side-by-side
title_txt <- if (direction == "outgoing") {
  paste0("KO vs WT: Top 10 PV→Oligo (outgoing) per cluster (p < ", thresh, ")")
} else {
  paste0("KO vs WT: Top 10 Oligo→PV (incoming) per cluster (p < ", thresh, ")")
}

(p_top_inc | p_top_dec) +
  patchwork::plot_annotation(
    title    = title_txt,
    subtitle = "Left: increases (red); Right: decreases (blue)"
  )

# ---- Optional save
ggsave("differential_bubble_PV_to_oligo_pnn.pdf", width = 16, height = 8)

# Snapshot outgoing results
top_inc_outgoing <- top_inc_per_type
top_dec_outgoing <- top_dec_per_type



#### Save the matrix (both directions) ----
library(dplyr)
library(readr)

cols_keep <- c("Flow","Direction","rank","CellType","LRP","Pathway",
               "Source","Target","prob_WT","prob_KO","delta")

prepare_export <- function(top_inc, top_dec, flow_tag) {
  inc <- top_inc %>%
    dplyr::mutate(Direction = "Increase", Flow = flow_tag, abs_delta = abs(delta)) %>%
    dplyr::arrange(CellType, dplyr::desc(abs_delta)) %>%
    dplyr::group_by(CellType) %>%
    dplyr::mutate(rank = dplyr::row_number()) %>%
    dplyr::ungroup()
  
  dec <- top_dec %>%
    dplyr::mutate(Direction = "Decrease", Flow = flow_tag, abs_delta = abs(delta)) %>%
    dplyr::arrange(CellType, dplyr::desc(abs_delta)) %>%
    dplyr::group_by(CellType) %>%
    dplyr::mutate(rank = dplyr::row_number()) %>%
    dplyr::ungroup()
  
  dplyr::bind_rows(inc, dec) %>%
    dplyr::arrange(CellType, Direction, rank) %>%
    dplyr::select(dplyr::all_of(cols_keep))
}

all_export <- dplyr::bind_rows(
  prepare_export(top_inc_incoming, top_dec_incoming, "Oligo→PV"),
  prepare_export(top_inc_outgoing, top_dec_outgoing, "PV→Oligo")
)

readr::write_csv(all_export, "top10_LRP_all_both_directions.csv")



##########
### Complex heatmaps ----------
# Data pre-processing
object.list_cortical_pnn <- list(WT = cellchat_wt_cortical_pnn, KO = cellchat_ko_cortical_pnn)
length(object.list_cortical_pnn[[1]]@netP$pathways)
length(object.list_cortical_pnn[[2]]@netP$pathways)

# # Also check the union you pass to signaling:
# length(pathway.union_cortical_pnn)
# print(pathway.union_cortical_pnn)

# Re-check
object.list_cortical_pnn <- list(WT = cellchat_wt_cortical_pnn, KO = cellchat_ko_cortical_pnn)

for (nm in names(object.list_cortical_pnn)) {
  x <- object.list_cortical_pnn[[nm]]
  
  # If you have not done these on the subset, re-run:
  x <- subsetData(x)  # ensures data slots match current idents
  x <- identifyOverExpressedGenes(x)
  x <- identifyOverExpressedInteractions(x)
  # projectData(...) if you used it originally (e.g., human PPI) 
  # x <- projectData(x, PPI.human)  # uncomment if needed
  
  x <- computeCommunProb(x)                  # computes communication probability
  x <- filterCommunication(x)                # filters low-probability interactions
  x <- computeCommunProbPathway(x)           # builds pathway-level network in x@netP
  x <- aggregateNet(x)                       # aggregates net and netP
  x <- netAnalysis_computeCentrality(x, slot.name = "netP")  # centrality for pathways
  
  object.list_cortical_pnn[[nm]] <- x
}

# Side-by-side and one-plot

i = 1
pathway.union_cortical_pnn <- union(object.list_cortical_pnn[[i]]@netP$pathways, object.list_cortical_pnn[[i+1]]@netP$pathways)

# Outgoing
# combining all the identified signaling pathways from different datasets 

ht1_out = netAnalysis_signalingRole_heatmap(object.list_cortical_pnn[[i]], pattern = "outgoing", signaling = pathway.union_cortical_pnn, 
                                            title = names(object.list_cortical_pnn)[i], width = 8, height = 20, color.heatmap = "GnBu")
ht2_out = netAnalysis_signalingRole_heatmap(object.list_cortical_pnn[[i+1]], pattern = "outgoing", signaling = pathway.union_cortical_pnn, 
                                            title = names(object.list_cortical_pnn)[i+1], width = 8, height = 20, color.heatmap = "GnBu")
draw(ht1_out + ht2_out, ht_gap = unit(1, "cm"))

# Extract the matrix data from the heatmaps
mat1_out <- ht1_out@matrix
mat2_out <- ht2_out@matrix

# Calculate the difference matrix
mat_diff_out <- mat1_out - mat2_out
mat_diff_out[is.na(mat_diff_out)] <- 0
mat_diff_out[is.infinite(mat_diff_out)] <- 0
col_fun <- colorRamp2(
  c(min(mat_diff_out), 0, max(mat_diff_out)),
  c("#2166ac", "white", "#b2182b")
)

# Create a heatmap for the difference
ht_diff_out <- Heatmap(
  mat_diff_out,
  name = "Difference",
  col = col_fun,
  width = unit(19, "cm"),
  height = unit(30, "cm")
)

# Draw the difference heatmap
draw(ht_diff_out)

# Incoming
ht1_in = netAnalysis_signalingRole_heatmap(object.list_cortical_pnn[[i]], pattern = "incoming", signaling = pathway.union_cortical_pnn, title = names(object.list_cortical_pnn)[i], width = 8, height = 20, color.heatmap = "GnBu")
ht2_in = netAnalysis_signalingRole_heatmap(object.list_cortical_pnn[[i+1]], pattern = "incoming", signaling = pathway.union_cortical_pnn, title = names(object.list_cortical_pnn)[i+1], width = 8, height = 20, color.heatmap = "GnBu")
draw(ht1_in + ht2_in, ht_gap = unit(1, "cm"))

# Extract the matrix data from the heatmaps
mat1_in <- ht1_in@matrix
mat2_in <- ht2_in@matrix

# Calculate the difference matrix
mat_diff_in <- mat1_in - mat2_in
mat_diff_in[is.na(mat_diff_in)] <- 0
mat_diff_in[is.infinite(mat_diff_in)] <- 0
col_fun <- colorRamp2(
  c(min(mat_diff_in), 0, max(mat_diff_in)),
  c("#2166ac", "white", "#b2182b")
)

# Create a heatmap for the difference
ht_diff_in <- Heatmap(
  mat_diff_in,
  name = "Difference",
  col = col_fun,
  width = unit(19, "cm"),
  height = unit(30, "cm")
)

# Draw the difference heatmap
draw(ht_diff_in)

# Both outgoing and incoming
ht1_all = netAnalysis_signalingRole_heatmap(object.list_cortical_pnn[[i]], pattern = "all", signaling = pathway.union_cortical_pnn, title = names(object.list_cortical_pnn)[i], width = 8, height = 20, color.heatmap = "OrRd")
ht2_all = netAnalysis_signalingRole_heatmap(object.list_cortical_pnn[[i+1]], pattern = "all", signaling = pathway.union_cortical_pnn, title = names(object.list_cortical_pnn)[i+1], width = 8, height = 20, color.heatmap = "OrRd")
draw(ht1_all + ht2_all, ht_gap = unit(1, "cm"))#

# Extract the matrix data from the heatmaps
mat1_all <- ht1_all@matrix
mat2_all <- ht2_all@matrix

# Calculate the difference matrix
mat_diff_all <- mat1_all - mat2_all
mat_diff_all[is.na(mat_diff_all)] <- 0
mat_diff_all[is.infinite(mat_diff_all)] <- 0
col_fun <- colorRamp2(
  c(min(mat_diff_all), 0, max(mat_diff_all)),
  c("#2166ac", "white", "#b2182b")
)

# Create a heatmap for the difference
ht_diff_all <- Heatmap(
  mat_diff_all,
  name = "Difference",
  col = col_fun,
  width = unit(19, "cm"),
  height = unit(30, "cm")
)

# Draw the difference heatmap
draw(ht_diff_all)


# PNG at 300 dpi, 11x14 inches
png("diff_heatmap_all_cortical_pnn.png", width = 11, height = 14, units = "in", res = 300)
draw(ht_diff_all)
dev.off()

png("diff_heatmap_in_cortical_pnn.png", width = 11, height = 14, units = "in", res = 300)
draw(ht_diff_in)
dev.off()

png("diff_heatmap_out_cortical_pnn.png", width = 11, height = 14, units = "in", res = 300)
draw(ht_diff_out)
dev.off()



### Complex heatmaps reduced to PNN relevant stuff ----
### Focused pathways (order preserved)
target_pathways <- c(
  "TENASCIN","SEMA3","RELN","NRG","CSPG4","CNTN","L1CAM","EPHA","EPHB","LAMININ"
)

# Compute available pathways across datasets and focus on targets
available_union <- Reduce(union, lapply(object.list_cortical_pnn, function(x) x@netP$pathways))
pathway.focus <- target_pathways[target_pathways %in% available_union]
if (length(pathway.focus) == 0) stop("None of the target pathways were found in netP$pathways")

# Convenience alias
i <- 1

### OUTGOING (focused)
ht1_out <- netAnalysis_signalingRole_heatmap(
  object.list_cortical_pnn[[i]], pattern = "outgoing",
  signaling = pathway.focus, title = names(object.list_cortical_pnn)[i],
  width = 8, height = 20, color.heatmap = "GnBu"
)
ht2_out <- netAnalysis_signalingRole_heatmap(
  object.list_cortical_pnn[[i+1]], pattern = "outgoing",
  signaling = pathway.focus, title = names(object.list_cortical_pnn)[i+1],
  width = 8, height = 20, color.heatmap = "GnBu"
)
draw(ht1_out + ht2_out, ht_gap = unit(1, "cm"))

# Difference (outgoing)
mat1_out <- ht1_out@matrix
mat2_out <- ht2_out@matrix
mat_diff_out <- mat1_out - mat2_out
mat_diff_out[is.na(mat_diff_out) | is.infinite(mat_diff_out)] <- 0
col_fun_out <- colorRamp2(c(min(mat_diff_out), 0, max(mat_diff_out)), c("#2166ac","white","#b2182b"))
ht_diff_out <- Heatmap(mat_diff_out, name = "Difference", col = col_fun_out,
                       width = unit(10, "cm"), height = unit(7, "cm"))
draw(ht_diff_out)

### INCOMING (focused)
ht1_in <- netAnalysis_signalingRole_heatmap(
  object.list_cortical_pnn[[i]], pattern = "incoming",
  signaling = pathway.focus, title = names(object.list_cortical_pnn)[i],
  width = 8, height = 20, color.heatmap = "GnBu"
)
ht2_in <- netAnalysis_signalingRole_heatmap(
  object.list_cortical_pnn[[i+1]], pattern = "incoming",
  signaling = pathway.focus, title = names(object.list_cortical_pnn)[i+1],
  width = 8, height = 20, color.heatmap = "GnBu"
)
draw(ht1_in + ht2_in, ht_gap = unit(1, "cm"))

# Difference (incoming)
mat1_in <- ht1_in@matrix
mat2_in <- ht2_in@matrix
mat_diff_in <- mat1_in - mat2_in
mat_diff_in[is.na(mat_diff_in) | is.infinite(mat_diff_in)] <- 0
col_fun_in <- colorRamp2(c(min(mat_diff_in), 0, max(mat_diff_in)), c("#2166ac","white","#b2182b"))
ht_diff_in <- Heatmap(mat_diff_in, name = "Difference", col = col_fun_in,
                      width = unit(10, "cm"), height = unit(7, "cm"))
draw(ht_diff_in)

### ALL (focused)
ht1_all <- netAnalysis_signalingRole_heatmap(
  object.list_cortical_pnn[[i]], pattern = "all",
  signaling = pathway.focus, title = names(object.list_cortical_pnn)[i],
  width = 8, height = 20, color.heatmap = "OrRd"
)
ht2_all <- netAnalysis_signalingRole_heatmap(
  object.list_cortical_pnn[[i+1]], pattern = "all",
  signaling = pathway.focus, title = names(object.list_cortical_pnn)[i+1],
  width = 8, height = 20, color.heatmap = "OrRd"
)
draw(ht1_all + ht2_all, ht_gap = unit(1, "cm"))

# Difference (all)
mat1_all <- ht1_all@matrix
mat2_all <- ht2_all@matrix
mat_diff_all <- mat1_all - mat2_all
mat_diff_all[is.na(mat_diff_all) | is.infinite(mat_diff_all)] <- 0
col_fun_all <- colorRamp2(c(min(mat_diff_all), 0, max(mat_diff_all)), c("#2166ac","white","#b2182b"))
ht_diff_all <- Heatmap(mat_diff_all, name = "Difference", col = col_fun_all,
                       width = unit(10, "cm"), height = unit(7, "cm"))
draw(ht_diff_all)


png("diff_heatmap_all_cortical_pnn_PNN.png", width = 8, height = 8, units = "in", res = 300)
draw(ht_diff_all); dev.off()
png("diff_heatmap_in_cortical_pnn_PNN.png", width = 8, height = 8, units = "in", res = 300)
draw(ht_diff_in); dev.off()
png("diff_heatmap_out_cortical_pnn_PNN.png", width = 8, height = 8, units = "in", res = 300)
draw(ht_diff_out); dev.off()


###########
### Scatter plot of interaction strength -------------
num.link <- sapply(object.list_cortical_pnn, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list_cortical_pnn)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list_cortical_pnn[[i]], title = names(object.list_cortical_pnn)[i], weight.MinMax = weight.MinMax)
}
patchwork::wrap_plots(plots = gg)



### Scatter plot of cluster-based signaling changes ----------------------
gg1 <- netAnalysis_signalingChanges_scatter(cellchat_merged_cortical_pnn, idents.use = "PV+ Int", signaling.exclude = "MIF")
#gg2 <- netAnalysis_signalingChanges_scatter(cellchat_merged_cortical_pnn, idents.use = "SST+ Int", signaling.exclude = "MIF")
#gg3 <- netAnalysis_signalingChanges_scatter(cellchat_merged_cortical_pnn, idents.use = "VIP+ Int", signaling.exclude = "MIF")

gg4 <- netAnalysis_signalingChanges_scatter(cellchat_merged_cortical_pnn, idents.use = "OPC", signaling.exclude = "MIF")
gg5 <- netAnalysis_signalingChanges_scatter(cellchat_merged_cortical_pnn, idents.use = "New Oligo", signaling.exclude = "MIF")
gg6 <- netAnalysis_signalingChanges_scatter(cellchat_merged_cortical_pnn, idents.use = "Myelinating Oligo", signaling.exclude = "MIF")
gg7 <- netAnalysis_signalingChanges_scatter(cellchat_merged_cortical_pnn, idents.use = "Astrocyte", signaling.exclude = "MIF")
gg8 <- netAnalysis_signalingChanges_scatter(cellchat_merged_cortical_pnn, idents.use = "Perineuronal Oligo", signaling.exclude = "MIF")
print(gg1)
#patchwork::wrap_plots(plots = list(gg1,gg2,gg3))
patchwork::wrap_plots(plots = list(gg4,gg8, gg5,gg6))




##########
### Extra chord diagram --------------------------
par(mfrow = c(1, 2), xpd=TRUE)
# compare all the interactions sending from PVs to OPCs
for (i in 1:length(object.list_cortical_pnn)) {
  netVisual_chord_gene(object.list_cortical_pnn[[i]], sources.use = 9, targets.use = c(14), #lab.cex = 0.3, 
                       title.name = paste0("Signaling from PV+ interneurons - ", names(object.list_cortical_pnn)[i]), small.gap = 4.5, big.gap=8)
}

# show all the significant signaling pathways from PV INs to oligos
par(mfrow = c(1, 2), xpd=TRUE)
for (i in 1:length(object.list_cortical_pnn)) {
  netVisual_chord_gene(object.list_cortical_pnn[[i]], sources.use = c(9), targets.use = c(13, 14,15,16),slot.name = "netP", title.name = paste0("Signaling pathways sending from PV+ interneurons - ", names(object.list_cortical_pnn)[i]), legend.pos.x = 10, small.gap = 5, big.gap=8)
}





###########
#### Circle interaction map (fix attempt)---------
suppressPackageStartupMessages({
  library(CellChat)
  library(igraph)
})

par(mfrow = c(1, 1), xpd = NA)

edge_width_max <- 2
weight_cap     <- 0.2
arrow_width    <- edge_width_max / 2
arrow_size     <- 0.3 * edge_width_max
label_offset_in <- 0.6

# Requested order (cluster 9 first); pass either indices or names
sources.use <- c(9,13,14,15,16)
targets.use <- c(9,13,14,15,16)

# Inputs
object  <- cellchat_merged_cortical_pnn
comp    <- c(1, 2)
measure <- "weight"  # "count" or "weight"

# Access nets and ensure dimnames via idents
net1 <- object@net[[comp[1]]]
net2 <- object@net[[comp[2]]]
mat1 <- net1[[measure]]
mat2 <- net2[[measure]]

fix_dimnames <- function(m, idents_levels) {
  if (is.null(rownames(m)) && length(idents_levels) == nrow(m)) rownames(m) <- idents_levels
  if (is.null(colnames(m)) && length(idents_levels) == ncol(m)) colnames(m) <- idents_levels
  m
}
lvl1 <- if (is.list(object@idents)) levels(object@idents[[comp[1]]]) else levels(object@idents)
lvl2 <- if (is.list(object@idents)) levels(object@idents[[comp[2]]]) else levels(object@idents)
mat1 <- fix_dimnames(mat1, lvl1)
mat2 <- fix_dimnames(mat2, lvl2)

group_names_all <- rownames(mat1)
if (is.null(group_names_all)) stop("Group names are missing; ensure idents have levels or matrices have dimnames.")

# Map requested indices to names using the original names
map_to_names <- function(v, names_all) {
  if (is.numeric(v)) {
    if (any(v < 1 | v > length(names_all))) stop("Index out of range in sources.use/targets.use")
    names_all[v]
  } else {
    as.character(v)
  }
}
sources_names <- map_to_names(sources.use, group_names_all)
targets_names <- map_to_names(targets.use, group_names_all)

# Relevel idents so the first requested name becomes the first sector; reorder all 2D nets accordingly
first_label <- sources_names[1]
new_levels <- c(first_label, setdiff(group_names_all, first_label))

if (is.list(object@idents)) {
  for (k in seq_along(object@idents)) {
    cur_lev <- levels(object@idents[[k]])
    keep_lev <- new_levels[new_levels %in% cur_lev]
    object@idents[[k]] <- factor(object@idents[[k]], levels = keep_lev)
  }
} else {
  cur_lev <- levels(object@idents)
  keep_lev <- new_levels[new_levels %in% cur_lev]
  object@idents <- factor(object@idents, levels = keep_lev)
}

mat_names_2d <- c("count","sum","weight","count.merged","weight.merged")
for (i in seq_along(object@net)) {
  for (mn in intersect(mat_names_2d, names(object@net[[i]]))) {
    m <- object@net[[i]][[mn]]
    rn <- rownames(m); cn <- colnames(m)
    if (!is.null(rn) && !is.null(cn)) {
      keep <- new_levels[new_levels %in% rn & new_levels %in% cn]
      object@net[[i]][[mn]] <- m[keep, keep, drop = FALSE]
    }
  }
}

# Refresh matrices and valid group names after reordering
net1 <- object@net[[comp[1]]]
net2 <- object@net[[comp[2]]]
mat1 <- net1[[measure]]
mat2 <- net2[[measure]]
group_names_all <- rownames(mat1)

# Keep only requested names that still exist
sources_names <- sources_names[sources_names %in% group_names_all]
targets_names <- targets_names[targets_names %in% group_names_all]

# Plot with default CellChat colors and asp=1 for perfect circle
netVisual_diffInteraction(
  object             = object,
  comparison         = comp,
  measure            = measure,
  color.edge         = c("#b2182b", "#2166ac"),
  weight.scale       = FALSE,
  edge.weight.max    = weight_cap,
  edge.width.max     = edge_width_max,
  sources.use        = sources_names,
  targets.use        = targets_names,
  remove.isolate     = TRUE,
  vertex.label.cex   = 1e-6,
  vertex.label.color = NA,
  title.name         = "",
  margin             = 0.8,
  alpha.edge         = 1,
  shape              = "circle",
  arrow.width        = arrow_width,
  arrow.size         = arrow_size
)

# Function to wrap multi-word labels into two lines
wrap_label <- function(label, max_width = 15) {
  # Count number of words
  words <- strsplit(label, " ")[[1]]
  if (length(words) <= 1) {
    return(label)
  }
  # Use strwrap to split into lines
  wrapped <- strwrap(label, width = max_width)
  # Combine with newline
  paste(wrapped, collapse = "\n")
}

# Label placement aligned to plotted order among active nodes
node_names_sel <- unique(c(sources_names, setdiff(targets_names, sources_names)))
sub1 <- mat1[node_names_sel, node_names_sel, drop = FALSE]
sub2 <- mat2[node_names_sel, node_names_sel, drop = FALSE]
sub1z <- ifelse(is.na(sub1), 0, sub1)
sub2z <- ifelse(is.na(sub2), 0, sub2)

deg_any <- (rowSums(sub1z != 0) + colSums(sub1z != 0) +
              rowSums(sub2z != 0) + colSums(sub2z != 0)) > 0
active_names <- node_names_sel[deg_any]

order_active_names <- unique(c(intersect(sources_names, active_names),
                               setdiff(intersect(targets_names, active_names),
                                       intersect(sources_names, active_names))))
n_active <- length(order_active_names)

if (n_active > 0) {
  g_dummy <- igraph::make_empty_graph(n = n_active, directed = TRUE)
  coords  <- igraph::layout_in_circle(g_dummy, order = seq_len(n_active))
  x <- coords[,1]; y <- coords[,2]
  theta <- atan2(y, x); ux <- cos(theta); uy <- sin(theta)
  
  dx_in <- grconvertX(1, "inches", "user") - grconvertX(0, "inches", "user")
  dy_in <- grconvertY(1, "inches", "user") - grconvertY(0, "inches", "user")
  vx <- ux * dx_in; vy <- uy * dy_in
  vnorm <- sqrt(vx^2 + vy^2)
  x_lab <- x + label_offset_in * vx / vnorm
  y_lab <- y + label_offset_in * vy / vnorm
  
  deg <- theta * 180 / pi
  pos_right  <- (deg >= -45 & deg < 45)
  pos_top    <- (deg >= 45  & deg < 135)
  pos_left   <- (deg >= 135 | deg <= -135)
  pos_bottom <- (deg > -135 & deg < -45)
  
  adjx <- rep(0.5, n_active); adjy <- rep(0.5, n_active)
  adjx[pos_right]  <- 0;   adjy[pos_right]  <- 0.5
  adjx[pos_left]   <- 0.9; adjy[pos_left]   <- 0.5
  adjx[pos_top]    <- 0.5; adjy[pos_top]    <- 0
  adjx[pos_bottom] <- 0.5; adjy[pos_bottom] <- 1
  
  # Wrap labels for multi-word labels
  wrapped_labels <- sapply(order_active_names, wrap_label, max_width = 15)
  
  for (i in seq_len(n_active)) {
    text(
      x = x_lab[i], y = y_lab[i], labels = wrapped_labels[i],
      adj = c(adjx[i], adjy[i]), srt = 0,
      cex = 1.3, col = "black", xpd = NA
    )
  }
}






##########
# Differential chord diagram PV to Oligo (SAVES ALSO LOCALLY)----
library(dplyr)
library(circlize)
library(ComplexHeatmap)
library(RColorBrewer)
library(scales)
library(grid)

# --- Data extraction and preprocessing ---
comm_list <- subsetCommunication(cellchat_merged_cortical_pnn)
comm1 <- comm_list[[1]]
comm2 <- comm_list[[2]]
comm_col <- if ("prob" %in% colnames(comm1)) "prob" else if ("value" %in% colnames(comm1)) "value" else colnames(comm1)[3]
signaling_col <- if ("interaction_name" %in% colnames(comm1)) "interaction_name" else if ("pathway_name" %in% colnames(comm1)) "pathway_name" else NA

comm1_sub <- comm1[, c("source", "target", signaling_col, comm_col)]
comm2_sub <- comm2[, c("source", "target", signaling_col, comm_col)]
names(comm1_sub)[4] <- "val1"
names(comm2_sub)[4] <- "val2"
diff_data <- full_join(comm1_sub, comm2_sub,
                       by = c("source", "target", signaling_col)) %>%
  mutate(val1 = replace_na(val1, 0),
         val2 = replace_na(val2, 0),
         diff = val2 - val1)

# --- Biologically implausible ligand families to exclude ---
exclude_ligand_families <- c(
  "GLUTAMATE",
  "CHOLESTEROL",
  "DEHYDROEPIANDROSTERONE",
  "PPIA",
  "TUB",
  "SEROTONINDOPAMIN"
)

# Helper: extract ligand core for filtering (same logic as sectors_df)
extract_ligand_core_filter <- function(interaction) {
  if (is.na(interaction)) return(NA_character_)
  ligand <- sub("_.*", "", interaction)
  ligand_core <- sub("[0-9]+$", "", ligand)
  toupper(ligand_core)
}

map_ligand_group_filter <- function(core) {
  if (is.na(core)) return("UNKNOWN")
  core_up <- toupper(core)
  if (grepl("^LAM[A-C]", core_up)) return("LAMININ")
  if (grepl("^LRRC",    core_up)) return("LRRC")
  if (grepl("^SEMA",    core_up)) return("SEMA")
  fallback <- sub("([A-Z]+).*", "\\1", core_up)
  if (nzchar(fallback)) return(fallback)
  return("OTHER")
}

# Tag each interaction with its ligand family and exclude
diff_data <- diff_data %>%
  mutate(
    ligand_core_tmp = vapply(get(signaling_col), extract_ligand_core_filter, character(1)),
    ligand_fam_tmp  = vapply(ligand_core_tmp,    map_ligand_group_filter,    character(1))
  ) %>%
  filter(!ligand_fam_tmp %in% exclude_ligand_families) %>%
  dplyr::select(-ligand_core_tmp, -ligand_fam_tmp)

# --- Subset cell populations ---
sources.use <- c("PV+ Int")
targets.use <- c("OPC", "NFOL", "PNOL","MOL")

df <- diff_data %>%
  filter(source %in% sources.use, target %in% targets.use, diff != 0) %>%
  mutate(pathway = get(signaling_col),
         sector_from = paste("FROM", source, pathway, sep = "|"),
         sector_to   = paste("TO", target, pathway, sep = "|"))


# --- Unique sectors (nodes) ---
sectors_in_df <- c(
  df %>% transmute(sector = sector_from) %>% pull(sector) %>% unique(),
  df %>% transmute(sector = sector_to)   %>% pull(sector) %>% unique()
) %>% unique()

# --- Parse sectors and derive ligand core ---
parse_sector <- function(sec) {
  x <- strsplit(sec, "\\|")[[1]]
  data.frame(
    sector = sec,
    side = x[1],
    celltype = x[2],
    pathway = x[3],
    stringsAsFactors = FALSE
  )
}
sectors_df <- do.call(rbind, lapply(sectors_in_df, parse_sector))

extract_ligand_core <- function(interaction) {
  if (is.na(interaction)) return(NA_character_)
  ligand <- sub("_.*", "", interaction)     # before underscore
  ligand_core <- sub("[0-9]+$", "", ligand) # drop trailing digits
  toupper(ligand_core)
}
sectors_df <- sectors_df %>%
  mutate(ligand_core = vapply(pathway, extract_ligand_core, character(1)))

# --- Group similar ligands into families ---
map_ligand_group <- function(core) {
  if (is.na(core)) return("UNKNOWN")
  core_up <- toupper(core)
  if (grepl("^LAM[A-C]", core_up)) return("LAMININ")  # LAMA/LAMB/LAMC
  if (grepl("^LRRC",    core_up)) return("LRRC")      # LRRC/LRRC4B...
  if (grepl("^SEMA",    core_up)) return("SEMA")      # SEMA4D/SEMA5A/SEMA6A...
  if (grepl("^PCDH", core_up)) return("PCDH")         # Protocadherins
  if (core_up == "L1CAM") return("L1CAM")  # ensure L1CAM stays intact
  fallback <- sub("([A-Z]+).*", "\\1", core_up)       # alphabetic stem
  if (nzchar(fallback)) return(fallback)
  return("OTHER")
}
sectors_df <- sectors_df %>%
  mutate(ligand_group = vapply(ligand_core, map_ligand_group, character(1)))

# --- Weights for ordering (sum of |diff| touching a sector) ---
sector_weights <- bind_rows(
  df %>% group_by(sector = sector_from) %>% summarise(w = sum(abs(diff)), .groups = "drop"),
  df %>% group_by(sector = sector_to)   %>% summarise(w = sum(abs(diff)), .groups = "drop")
) %>% group_by(sector) %>% summarise(w = sum(w), .groups = "drop")

sectors_df <- sectors_df %>% left_join(sector_weights, by = "sector")
sectors_df$w[is.na(sectors_df$w)] <- 0

# --- Per-cell type ligand-family weights for ordering ---
group_weights_ct <- sectors_df %>%
  group_by(celltype, ligand_group) %>%
  summarise(group_w = sum(w), .groups = "drop")

# --- Order sectors: biggest ligand-family first within each cell type ---
celltype_order <- unique(sectors_df$celltype)
ordered_sectors <- character(0)

for (ct in celltype_order) {
  gtab <- group_weights_ct %>% filter(celltype == ct) %>% arrange(desc(group_w))
  groups_ct <- gtab$ligand_group
  secs_ct <- sectors_df %>% filter(celltype == ct)
  for (g in groups_ct) {
    secs_ct_g <- secs_ct %>% filter(ligand_group == g)
    from_first <- secs_ct_g %>% filter(side == "FROM") %>% arrange(pathway, sector)
    to_next    <- secs_ct_g %>% filter(side == "TO")   %>% arrange(pathway, sector)
    ordered_sectors <- c(ordered_sectors, from_first$sector, to_next$sector)
  }
}
ordered_sectors <- intersect(unique(ordered_sectors), sectors_in_df)

# --- Colors ---
# Cell type colors (user defaults)
celltype_colors_user <- c(
  "PV+ Int" = "#DD4124",
  "OPC" = "#4BACC6",
  "PNOL" = "#ea8a33",
  "NFOL" = "#E8A0BF",
  "MOL" = "#9B2335"
)
missing_ct <- setdiff(celltype_order, names(celltype_colors_user))
if (length(missing_ct) > 0) {
  extra_cols <- colorRampPalette(brewer.pal(8, "Set1"))(length(missing_ct))
  names(extra_cols) <- missing_ct
  celltype_colors <- c(celltype_colors_user, extra_cols)
} else {
  celltype_colors <- celltype_colors_user
}
celltype_colors <- celltype_colors[celltype_order]

# --- Very distinct ligand-family colors (avoid black and very dark) ---
# Convert hex to Lab using base grDevices
hex_to_lab <- function(hex_vec) {
  rgb01 <- t(grDevices::col2rgb(hex_vec) / 255)
  grDevices::convertColor(rgb01, from = "sRGB", to = "Lab")
}
exclude_dark <- function(cols, L_min = 20) {
  lab <- hex_to_lab(cols)
  cols[!is.na(lab[, 1]) & lab[, 1] >= L_min]
}

# Okabe–Ito without black/near-black
okabe_ito_nonblack <- function(maxn = 9) {
  pal <- tryCatch(grDevices::palette.colors(maxn, palette = "Okabe-Ito"),
                  error = function(e) NULL)
  if (is.null(pal)) return(character(0))
  pal <- pal[!tolower(pal) %in% c("#000000", "black")]
  pal <- exclude_dark(pal, L_min = 20)
  unique(pal)
}

# Dense HCL candidates (vivid, within sRGB gamut)
hcl_candidates <- function(m_per_ring = 720) {
  h <- seq(0, 360 - 360/m_per_ring, length.out = m_per_ring)
  c1 <- grDevices::hcl(h = h, c = 70, l = 65)
  c2 <- grDevices::hcl(h = (h + 180/m_per_ring) %% 360, c = 70, l = 55)
  unique(c(c1, c2))
}

# Greedy farthest-point selection in Lab with optional seeds
select_farthest_lab_seeded <- function(candidates_hex, n, seed_hex = character(0)) {
  candidates_hex <- unique(candidates_hex)
  candidates_hex <- candidates_hex[!is.na(candidates_hex)]
  selected <- unique(seed_hex)
  selected <- selected[selected %in% candidates_hex]
  remaining <- setdiff(candidates_hex, selected)
  if (length(selected) == 0 && length(remaining) > 0) {
    lab <- hex_to_lab(remaining)
    dmat <- as.matrix(dist(lab))
    start <- remaining[which.max(rowMeans(dmat))]
    selected <- c(selected, start)
    remaining <- setdiff(remaining, start)
  }
  while (length(selected) < n && length(remaining) > 0) {
    lab_sel <- hex_to_lab(selected)
    lab_rem <- hex_to_lab(remaining)
    d <- as.matrix(dist(rbind(lab_sel, lab_rem)))
    ns <- nrow(lab_sel); nr <- nrow(lab_rem)
    d_sr <- d[seq_len(ns), ns + seq_len(nr), drop = FALSE]
    mind <- apply(d_sr, 2, min)
    next_idx <- which.max(mind)
    selected <- c(selected, remaining[next_idx])
    remaining <- remaining[-next_idx]
  }
  selected[seq_len(min(n, length(selected)))]
}

# Reorder colors to maximize adjacency contrast
order_max_contrast <- function(hex_vec) {
  if (length(hex_vec) <= 2) return(hex_vec)
  lab <- hex_to_lab(hex_vec)
  dmat <- as.matrix(dist(lab))
  start <- which.max(rowMeans(dmat))
  order_idx <- start
  remaining <- setdiff(seq_along(hex_vec), start)
  current <- start
  while (length(remaining) > 0) {
    nxt <- remaining[which.max(dmat[current, remaining])]
    order_idx <- c(order_idx, nxt)
    remaining <- setdiff(remaining, nxt)
    current <- nxt
  }
  hex_vec[order_idx]
}

# Build distinct palette: avoid black, maximize pairwise Lab distance
make_distinct_ligand_palette <- function(n) {
  seed <- if (n <= 9) head(okabe_ito_nonblack(9), n) else character(0)
  cand <- unique(c(okabe_ito_nonblack(9), hcl_candidates(720)))
  cand <- exclude_dark(cand, L_min = 20)
  sel  <- select_farthest_lab_seeded(cand, n, seed_hex = seed)
  order_max_contrast(sel)
}

# Identify ligand groups and assign colors with adjacency contrast along ring
ligand_groups_all <- unique(sectors_df$ligand_group)
n_lg <- length(ligand_groups_all)
palette_distinct <- make_distinct_ligand_palette(n_lg)

# Determine the appearance order of ligand groups along the ring using FROM sectors
ligand_group_per_sector_all <- setNames(sectors_df$ligand_group, sectors_df$sector)
side_per_sector_all         <- setNames(sectors_df$side,         sectors_df$sector)
lg_ring   <- ligand_group_per_sector_all[ordered_sectors]
side_ring <- side_per_sector_all[ordered_sectors]
ring_groups <- unique(lg_ring[side_ring == "FROM"])
ring_groups <- ring_groups[!is.na(ring_groups)]

# Map colors to ring groups to maximize adjacent contrast; remaining groups follow
assign_groups <- c(ring_groups, setdiff(ligand_groups_all, ring_groups))
ligand_group_colors <- setNames(palette_distinct[seq_along(assign_groups)], assign_groups)

# --- Chord node colors (for sector grids; bands come from highlights)
chord_colors <- setNames(
  colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(length(ordered_sectors)),
  ordered_sectors
)

# --- Edge appearance ---
color.edge <- c("#b2182b", "#2166ac")
link_col <- ifelse(df$diff >= 0, color.edge[1], color.edge[2])
max_width <- 10
link_lwd <- if (nrow(df) > 0) rescale(abs(df$diff), to = c(1, max_width)) else numeric()

# --- Gaps based on new order (small within-ct, larger between-ct) ---
sector_celltype_ordered <- setNames(sectors_df$celltype, sectors_df$sector)[ordered_sectors]
gap_vec <- ifelse(
  c(TRUE, diff(as.numeric(factor(sector_celltype_ordered, levels = celltype_order))) == 0),
  0.2, 4
)

# Name gaps to sector order (add this line)
names(gap_vec) <- ordered_sectors

# --- Plotting with spacers and selective inner ring (sources only) ---

out_pdf <- "differential_chord_PV_oligos.pdf"
pdf(out_pdf, width = 13, height = 9)

circos.clear()
circos.par(
  start.degree = 90,
  gap.after = gap_vec,
  cell.padding = c(0, 0, 0, 0),
  track.margin = c(0, 0),
  points.overflow.warning = FALSE,
  unit.circle.segments = 500
)

if (nrow(df) > 0 && length(ordered_sectors) > 1) {
  chordDiagram(
    x = df[, c("sector_from", "sector_to", "diff")],
    grid.col = chord_colors[ordered_sectors],
    directional = 1,
    direction.type = c("arrows"),
    link.arr.type = "big.arrow",
    link.arr.length = 0.20,
    transparency = 0.6,
    link.lwd = link_lwd,
    link.border = NA,
    col = link_col,
    order = ordered_sectors,
    annotationTrack = NULL,
    reduce = 0,
    preAllocateTracks = list(
      list(track.height = 0.09, bg.border = NA),         # 1: outer cell-type ring
      list(track.height = circlize::mm_h(2), bg.border = NA),      # 2: spacer outer-inner
      list(track.height = 0.08, bg.border = NA),         # 3: inner ligand-family ring
      list(track.height = circlize::mm_h(1.5), bg.border = NA)     # 4: spacer inner-links
    )
  )
  
  # --- Outer cell-type band: continuous per cell type (track 1) ---
  for (ct in celltype_order) {
    secs_ct <- ordered_sectors[sector_celltype_ordered == ct]
    if (length(secs_ct) > 0) {
      circlize::highlight.sector(
        sector.index = secs_ct,
        track.index = 1,
        col = celltype_colors[ct],
        border = NA,
        padding = c(0, 0, 0, 0)
      )
    }
  }
  
  # --- Inner ligand-family band: draw ONLY for sources (FROM), not targets (track 3) ---
  ligand_group_per_sector <- setNames(sectors_df$ligand_group, sectors_df$sector)
  side_per_sector         <- setNames(sectors_df$side,         sectors_df$sector)
  ligand_group_per_sector <- ligand_group_per_sector[ordered_sectors]
  side_per_sector         <- side_per_sector[ordered_sectors]
  
  for (lg in names(ligand_group_colors)) {
    secs_lg_from <- names(ligand_group_per_sector)[
      ligand_group_per_sector == lg & side_per_sector == "FROM"
    ]
    if (length(secs_lg_from) > 0) {
      circlize::highlight.sector(
        sector.index = secs_lg_from,
        track.index  = 3,
        col          = ligand_group_colors[lg],
        border       = NA,
        padding      = c(0, 0, 0, 0)
      )
    }
  }
  
  # --- Legends (packed and movable) ---
  edge_legend <- ComplexHeatmap::Legend(
    labels = ComplexHeatmap::gt_render(c(
      "Increased in <i>Dnmt1</i>-KO",
      "Decreased in <i>Dnmt1</i>-KO"
    )),
    legend_gp = grid::gpar(fill = color.edge),
    labels_gp = grid::gpar(fontsize = 10),
    title = "Edge"
  )
  ct_legend <- ComplexHeatmap::Legend(
    labels = names(celltype_colors),
    legend_gp = grid::gpar(fill = unname(celltype_colors), col = unname(celltype_colors)),
    title = "Cell Population"
  )
  ligand_legend <- ComplexHeatmap::Legend(
    labels = names(ligand_group_colors),
    legend_gp = grid::gpar(fill = unname(ligand_group_colors), col = unname(ligand_group_colors)),
    title = "Ligand Family"
  )
  pd <- ComplexHeatmap::packLegend(edge_legend, ct_legend, ligand_legend, direction = "vertical")
  ComplexHeatmap::draw(pd, x = grid::unit(0.98, "npc"), y = grid::unit(0.95, "npc"), just = c("right", "top"))
} else {
  message("Not enough edges or sectors to plot a chord diagram. Check your filter and input.")
}

dev.off()
message("Saved chord diagram to: ", normalizePath(out_pdf))

# --- Optional: print top pathway changes
top_up   <- df %>% arrange(-diff) %>% head(10)
top_down <- df %>% arrange(diff)  %>% head(10)
cat("\nMost increased signaling pairs (source, target, pathway, diff):\n")
print(top_up[, c("source", "target", signaling_col, "diff")])
cat("\nMost decreased signaling pairs:\n")
print(top_down[, c("source", "target", signaling_col, "diff")])

save(ligand_group_colors, file = "ligand_group_colors_PV_to_oligo.RData")



# Differential chord diagram Oligo to PV (SAVES ALSO LOCALLY)----
library(dplyr)
library(circlize)
library(ComplexHeatmap)
library(RColorBrewer)
library(scales)
library(grid)

# --- Data extraction and preprocessing ---
comm_list <- subsetCommunication(cellchat_merged_cortical_pnn)
comm1 <- comm_list[[1]]
comm2 <- comm_list[[2]]
comm_col <- if ("prob" %in% colnames(comm1)) "prob" else if ("value" %in% colnames(comm1)) "value" else colnames(comm1)[3]
signaling_col <- if ("interaction_name" %in% colnames(comm1)) "interaction_name" else if ("pathway_name" %in% colnames(comm1)) "pathway_name" else NA

comm1_sub <- comm1[, c("source", "target", signaling_col, comm_col)]
comm2_sub <- comm2[, c("source", "target", signaling_col, comm_col)]
names(comm1_sub)[4] <- "val1"
names(comm2_sub)[4] <- "val2"
diff_data <- full_join(comm1_sub, comm2_sub,
                       by = c("source", "target", signaling_col)) %>%
  mutate(val1 = replace_na(val1, 0),
         val2 = replace_na(val2, 0),
         diff = val2 - val1)


# --- Biologically implausible ligand families to exclude ---
exclude_ligand_families <- c(
  "GLUTAMATE",
  "CHOLESTEROL",
  "DEHYDROEPIANDROSTERONE",
  "PPIA",
  "TUB",
  "SEROTONINDOPAMIN"
)

# Helper: extract ligand core for filtering (same logic as sectors_df)
extract_ligand_core_filter <- function(interaction) {
  if (is.na(interaction)) return(NA_character_)
  ligand <- sub("_.*", "", interaction)
  ligand_core <- sub("[0-9]+$", "", ligand)
  toupper(ligand_core)
}

map_ligand_group_filter <- function(core) {
  if (is.na(core)) return("UNKNOWN")
  core_up <- toupper(core)
  if (grepl("^LAM[A-C]", core_up)) return("LAMININ")
  if (grepl("^LRRC",    core_up)) return("LRRC")
  if (grepl("^SEMA",    core_up)) return("SEMA")
  fallback <- sub("([A-Z]+).*", "\\1", core_up)
  if (nzchar(fallback)) return(fallback)
  return("OTHER")
}

# Tag each interaction with its ligand family and exclude
diff_data <- diff_data %>%
  mutate(
    ligand_core_tmp = vapply(get(signaling_col), extract_ligand_core_filter, character(1)),
    ligand_fam_tmp  = vapply(ligand_core_tmp,    map_ligand_group_filter,    character(1))
  ) %>%
  filter(!ligand_fam_tmp %in% exclude_ligand_families) %>%
  dplyr::select(-ligand_core_tmp, -ligand_fam_tmp)

# --- Subset cell populations (Flipped: targets become sources, sources become targets) ---
sources.use <- c("OPC", "NFOL", "PNOL","MOL")
targets.use <- c("PV+ Int")

df <- diff_data %>%
  filter(source %in% sources.use, target %in% targets.use, diff != 0) %>%
  mutate(pathway = get(signaling_col),
         sector_from = paste("FROM", source, pathway, sep = "|"),
         sector_to   = paste("TO", target, pathway, sep = "|"))


# --- Unique sectors (nodes) ---
sectors_in_df <- c(
  df %>% transmute(sector = sector_from) %>% pull(sector) %>% unique(),
  df %>% transmute(sector = sector_to)   %>% pull(sector) %>% unique()
) %>% unique()

# --- Parse sectors and ligand grouping ---
parse_sector <- function(sec) {
  x <- strsplit(sec, "\\|")[[1]]
  data.frame(
    sector = sec,
    side = x[1],
    celltype = x[2],
    pathway = x[3],
    stringsAsFactors = FALSE
  )
}
sectors_df <- do.call(rbind, lapply(sectors_in_df, parse_sector))

extract_ligand_core <- function(interaction) {
  if (is.na(interaction)) return(NA_character_)
  ligand <- sub("_.*", "", interaction)
  ligand_core <- sub("[0-9]+$", "", ligand)
  toupper(ligand_core)
}
sectors_df <- sectors_df %>%
  mutate(ligand_core = vapply(pathway, extract_ligand_core, character(1)))

# --- Updated: keep L1CAM intact ---
map_ligand_group <- function(core) {
  if (is.na(core)) return("UNKNOWN")
  core_up <- toupper(core)
  if (grepl("^LAM[A-C]", core_up)) return("LAMININ")
  if (grepl("^LRRC", core_up)) return("LRRC")
  if (grepl("^SEMA", core_up)) return("SEMA")
  if (core_up == "L1CAM") return("L1CAM")  # ensure L1CAM stays intact
  if (grepl("^PCDH", core_up)) return("PCDH")         # Protocadherins
  fallback <- sub("([A-Z]{2,}).*", "\\1", core_up)  # require ≥2 letters
  if (nzchar(fallback)) return(fallback)
  return("OTHER")
}
sectors_df <- sectors_df %>%
  mutate(ligand_group = vapply(ligand_core, map_ligand_group, character(1)))

# --- Weight computation and ordering preserved ---
sector_weights <- bind_rows(
  df %>% group_by(sector = sector_from) %>% summarise(w = sum(abs(diff)), .groups = "drop"),
  df %>% group_by(sector = sector_to)   %>% summarise(w = sum(abs(diff)), .groups = "drop")
) %>% group_by(sector) %>% summarise(w = sum(w), .groups = "drop")

sectors_df <- sectors_df %>% left_join(sector_weights, by = "sector")
sectors_df$w[is.na(sectors_df$w)] <- 0

group_weights_ct <- sectors_df %>%
  group_by(celltype, ligand_group) %>%
  summarise(group_w = sum(w), .groups = "drop")

celltype_order <- unique(sectors_df$celltype)
ordered_sectors <- character(0)
for (ct in celltype_order) {
  gtab <- group_weights_ct %>% filter(celltype == ct) %>% arrange(desc(group_w))
  groups_ct <- gtab$ligand_group
  secs_ct <- sectors_df %>% filter(celltype == ct)
  for (g in groups_ct) {
    secs_ct_g <- secs_ct %>% filter(ligand_group == g)
    from_first <- secs_ct_g %>% filter(side == "FROM") %>% arrange(pathway, sector)
    to_next    <- secs_ct_g %>% filter(side == "TO")   %>% arrange(pathway, sector)
    ordered_sectors <- c(ordered_sectors, from_first$sector, to_next$sector)
  }
}
ordered_sectors <- intersect(unique(ordered_sectors), sectors_in_df)

# --- Keep same colors for populations and ligand families ---
celltype_colors_user <- c(
  "PV+ Int" = "#DD4124",
  "OPC" = "#4BACC6",
  "PNOL" = "#ea8a33",
  "NFOL" = "#E8A0BF",
  "MOL" = "#9B2335"
)
missing_ct <- setdiff(celltype_order, names(celltype_colors_user))
if (length(missing_ct) > 0) {
  extra_cols <- colorRampPalette(brewer.pal(8, "Set1"))(length(missing_ct))
  names(extra_cols) <- missing_ct
  celltype_colors <- c(celltype_colors_user, extra_cols)
} else {
  celltype_colors <- celltype_colors_user
}
celltype_colors <- celltype_colors[celltype_order]

# --- Reuse ligand color map from PV -> Oligo if available ---
if (file.exists("ligand_group_colors_PV_to_oligo.RData")) {
  load("ligand_group_colors_PV_to_oligo.RData")
} else {
  warning("Using default new ligand color map; PV_to_oligo map not found.")
  ligand_groups_all <- unique(sectors_df$ligand_group)
  n_lg <- length(ligand_groups_all)
  palette_distinct <- colorRampPalette(brewer.pal(12, "Paired"))(n_lg)
  ligand_group_colors <- setNames(palette_distinct[seq_along(ligand_groups_all)], ligand_groups_all)
}

# Add any new ligand families (not present in the saved map)
missing <- setdiff(unique(sectors_df$ligand_group), names(ligand_group_colors))
if (length(missing) > 0) {
  add_cols <- colorRampPalette(brewer.pal(8, "Set2"))(length(missing))
  names(add_cols) <- missing
  ligand_group_colors <- c(ligand_group_colors, add_cols)
}

# --- Chord Plot ---
chord_colors <- setNames(colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(length(ordered_sectors)), ordered_sectors)
color.edge <- c("#b2182b", "#2166ac")
link_col <- ifelse(df$diff >= 0, color.edge[1], color.edge[2])
max_width <- 10
link_lwd <- if (nrow(df) > 0) rescale(abs(df$diff), to = c(1, max_width)) else numeric()
sector_celltype_ordered <- setNames(sectors_df$celltype, sectors_df$sector)[ordered_sectors]
gap_vec <- ifelse(c(TRUE, diff(as.numeric(factor(sector_celltype_ordered, levels = celltype_order))) == 0), 0.2, 4)
names(gap_vec) <- ordered_sectors

out_pdf <- "differential_chord_Oligos_to_PV.pdf"
pdf(out_pdf, width = 13, height = 9)

circos.clear()
circos.par(start.degree = 90, gap.after = gap_vec, cell.padding = c(0, 0, 0, 0),
           track.margin = c(0, 0), points.overflow.warning = FALSE,
           unit.circle.segments = 500)

if (nrow(df) > 0 && length(ordered_sectors) > 1) {
  chordDiagram(
    x = df[, c("sector_from", "sector_to", "diff")],
    grid.col = chord_colors[ordered_sectors],
    directional = 1,
    direction.type = c("arrows"),
    link.arr.type = "big.arrow",
    link.arr.length = 0.20,
    transparency = 0.6,
    link.lwd = link_lwd,
    link.border = NA,
    col = link_col,
    order = ordered_sectors,
    annotationTrack = NULL,
    reduce = 0,
    preAllocateTracks = list(
      list(track.height = 0.09, bg.border = NA),
      list(track.height = circlize::mm_h(2), bg.border = NA),
      list(track.height = 0.08, bg.border = NA),
      list(track.height = circlize::mm_h(1.5), bg.border = NA)
    )
  )
  
  for (ct in celltype_order) {
    secs_ct <- ordered_sectors[sector_celltype_ordered == ct]
    if (length(secs_ct) > 0) {
      circlize::highlight.sector(sector.index = secs_ct, track.index = 1,
                                 col = celltype_colors[ct], border = NA, padding = c(0, 0, 0, 0))
    }
  }
  
  ligand_group_per_sector <- setNames(sectors_df$ligand_group, sectors_df$sector)
  side_per_sector         <- setNames(sectors_df$side,         sectors_df$sector)
  ligand_group_per_sector <- ligand_group_per_sector[ordered_sectors]
  side_per_sector         <- side_per_sector[ordered_sectors]
  
  for (lg in names(ligand_group_colors)) {
    secs_lg_from <- names(ligand_group_per_sector)[ligand_group_per_sector == lg & side_per_sector == "FROM"]
    if (length(secs_lg_from) > 0) {
      circlize::highlight.sector(sector.index = secs_lg_from, track.index = 3,
                                 col = ligand_group_colors[lg], border = NA, padding = c(0, 0, 0, 0))
    }
  }
  
  edge_legend <- ComplexHeatmap::Legend(
    labels = ComplexHeatmap::gt_render(c("Increased in <i>Dnmt1</i>-KO", "Decreased in <i>Dnmt1</i>-KO")),
    legend_gp = grid::gpar(fill = color.edge), labels_gp = grid::gpar(fontsize = 10), title = "Edge"
  )
  ct_legend <- ComplexHeatmap::Legend(
    labels = names(celltype_colors),
    legend_gp = grid::gpar(fill = unname(celltype_colors), col = unname(celltype_colors)),
    title = "Cell Population"
  )
  ligand_legend <- ComplexHeatmap::Legend(
    labels = names(ligand_group_colors),
    legend_gp = grid::gpar(fill = unname(ligand_group_colors), col = unname(ligand_group_colors)),
    title = "Ligand Family"
  )
  pd <- ComplexHeatmap::packLegend(edge_legend, ct_legend, ligand_legend, direction = "vertical")
  ComplexHeatmap::draw(pd, x = grid::unit(0.98, "npc"),
                       y = grid::unit(0.95, "npc"), just = c("right", "top"))
} else {
  message("Not enough edges or sectors to plot a chord diagram. Check your filter and input.")
}

dev.off()
message("Saved chord diagram to: ", normalizePath(out_pdf))

top_up   <- df %>% arrange(-diff) %>% head(10)
top_down <- df %>% arrange(diff)  %>% head(10)
cat("\nMost increased signaling pairs (source, target, pathway, diff):\n")
print(top_up[, c("source", "target", signaling_col, "diff")])
cat("\nMost decreased signaling pairs:\n")
print(top_down[, c("source", "target", signaling_col, "diff")])





#####################
#### Violin plot for expression of genes across clusters (z score) ------
library(Seurat)
library(ggplot2)
library(dplyr)
library(colorspace)

seu <- pv_cortical_filtered_pnn

exclude_clusters <- c("Endothelial cells", "Leptomeningeal cells", "Meningeal fibroblasts")
cluster_names <- c(
  "Layer 2/3 IT neurons",
  "Layer 4 sensory neurons", 
  "Layer 5a IT neurons",
  "Layer 5b PT neurons",
  "Layer 5/6 IT neurons",
  "Layer 6 corticothalamic neurons",
  "Deep-layer extratelencephalic neurons",
  "Corticospinal neurons (Type I)",
  "Corticospinal neurons (Type II)",
  "Atypical excitatory neurons",
  "PV+ interneurons",
  "SST+ interneurons",
  "VIP+ interneurons",
  "Astrocytes",
  "Oligodendrocyte precursor cells",
  "Perineuronal oligodendrocytes",
  "Newly formed oligodendrocytes",
  "Myelinating oligodendrocytes",
  "Microglia"
)
ko_colors <- c(
  "#6B5B95", "#45B8AC", "#955251", "#4E84C4", "#B565A7", "#88B04B", "#C3447A",
  "#009B77", "#EFC050", "#7B6888", "#7FCDCD", "#DD4124", "#5B5EA6", "#E07A5F",
  "#4BACC6", "#ea8a33", "#E8A0BF", "#9B2335", "#C17BAE"
)
names(ko_colors) <- cluster_names
ctrl_colors <- sapply(ko_colors, colorspace::lighten, amount=0.4)
names(ctrl_colors) <- cluster_names
all_colors <- c(
  setNames(ctrl_colors, paste0(cluster_names, "_Ctrl")),
  setNames(ko_colors, paste0(cluster_names, "_KO"))
)

# Recode Genotype (assumes values are correct)
seu$Genotype_short <- recode(
  seu$Genotype,
  "PV-Cre/tdTom" = "Ctrl",
  "PV-Cre/tdTom/Dnmt1 loxP2" = "KO"
)

seu$ClusterName <- factor(
  seu$seurat_clusters,
  levels = 0:(length(cluster_names) - 1),
  labels = cluster_names
)

# Build dataframe and scale (z-score) the Dnmt1 expression
plotdf <- FetchData(seu, vars = c("ClusterName", "Genotype_short", "Dnmt1")) %>%
  filter(!is.na(ClusterName), !is.na(Genotype_short)) %>%
  filter(!ClusterName %in% exclude_clusters)
plotdf$pair <- paste0(plotdf$ClusterName, "_", plotdf$Genotype_short)

# Calculate z score across all cells for Dnmt1 (standard score)
plotdf$Dnmt1_zscore <- as.numeric(scale(plotdf$Dnmt1))

# Plot z-scored Dnmt1 values
ggplot(plotdf, aes(x=ClusterName, y=Dnmt1_zscore, fill=pair)) +
  geom_violin(position = position_dodge(width = 0.8), width = 0.75, scale = 'width', trim = TRUE) +
  scale_fill_manual(values = all_colors) +
  geom_jitter(
    position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
    size = 0.4, alpha = 0.6, color = "black", shape = 16
  ) +
  xlab("") +
  ylab("Dnmt1 Z-score (across all cells)") +
  ggtitle("Z-scored Dnmt1 expression across clusters and genotypes") +
  theme(
    axis.text.x = element_text(angle=75, hjust=1, vjust=1, color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", size = 0.8),
    panel.background = element_blank(),
    plot.background = element_blank(),
    legend.position="none"
  )




#### Violin plot for expression of genes across clusters (z score, ctrl-ko together)----
library(Seurat)
library(ggplot2)
library(dplyr)
library(colorspace)

seu <- pv_cortical_filtered_pnn

exclude_clusters <- c("Endothelial cells", "Leptomeningeal cells", "Meningeal fibroblasts")
cluster_names <- c(
  "Layer 2/3 IT neurons",
  "Layer 4 sensory neurons", 
  "Layer 5a IT neurons",
  "Layer 5b PT neurons",
  "Layer 5/6 IT neurons",
  "Layer 6 corticothalamic neurons",
  "Deep-layer extratelencephalic neurons",
  "Corticospinal neurons (Type I)",
  "Corticospinal neurons (Type II)",
  "Atypical excitatory neurons",
  "PV+ interneurons",
  "SST+ interneurons",
  "VIP+ interneurons",
  "Astrocytes",
  "Oligodendrocyte precursor cells",
  "Perineuronal oligodendrocytes",
  "Newly formed oligodendrocytes",
  "Myelinating oligodendrocytes",
  "Microglia"
)

cluster_colors <- c(
  "#6B5B95", "#45B8AC", "#955251", "#4E84C4", "#B565A7", "#88B04B", "#C3447A",
  "#009B77", "#EFC050", "#7B6888", "#7FCDCD", "#DD4124", "#5B5EA6", "#E07A5F",
  "#4BACC6", "#ea8a33", "#E8A0BF", "#9B2335", "#C17BAE"
)
names(cluster_colors) <- cluster_names

seu$ClusterName <- factor(
  seu$seurat_clusters,
  levels = 0:(length(cluster_names) - 1),
  labels = cluster_names
)

# Dataframe with only clusters, removing genotype dependency
plotdf <- FetchData(seu, vars = c("ClusterName", "Dnmt1")) %>%
  filter(!is.na(ClusterName)) %>%
  filter(!ClusterName %in% exclude_clusters)
plotdf$Dnmt1_zscore <- as.numeric(scale(plotdf$Dnmt1))

ggplot(plotdf, aes(x = ClusterName, y = Dnmt1_zscore, fill = ClusterName)) +
  geom_violin(width = 0.75, scale = 'width', trim = TRUE) +
  scale_fill_manual(values = cluster_colors) +
  geom_jitter(
    width = 0.15,
    size = 0.4, alpha = 0.6, color = "black", shape = 16
  ) +
  xlab("") +
  ylab("Dnmt1 Z-score (across all cells)") +
  ggtitle("Z-scored Dnmt1 expression across clusters") +
  theme(
    axis.text.x = element_text(angle = 75, hjust = 1, vjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", size = 0.8),
    panel.background = element_blank(),
    plot.background = element_blank(),
    legend.position = "none"
  )


#### Violin plot for expression of genes across clusters (z score, only ctrl)----
library(Seurat)
library(ggplot2)
library(dplyr)
library(colorspace)

seu <- pv_cortical_filtered_pnn

exclude_clusters <- c(
  "Layer 2/3 IT neurons",
  "Layer 4 sensory neurons", 
  "Layer 5a IT neurons",
  "Layer 5b PT neurons",
  "Layer 5/6 IT neurons",
  "Layer 6 corticothalamic neurons",
  "Deep-layer extratelencephalic neurons",
  "Corticospinal neurons (Type I)",
  "Corticospinal neurons (Type II)",
  "Atypical excitatory neurons",
  "SST+ interneurons",
  "VIP+ interneurons",
  "Perineuronal oligodendrocytes",
  "Microglia",
  "Endothelial cells", 
  "Leptomeningeal cells", 
  "Meningeal fibroblasts"
)

# exclude_clusters <- c(  "Endothelial cells", "Leptomeningeal cells", "Meningeal fibroblasts")

# All possible cluster names
cluster_names <- c(
  "Layer 2/3 IT neurons",
  "Layer 4 sensory neurons", 
  "Layer 5a IT neurons",
  "Layer 5b PT neurons",
  "Layer 5/6 IT neurons",
  "Layer 6 corticothalamic neurons",
  "Deep-layer extratelencephalic neurons",
  "Corticospinal neurons (Type I)",
  "Corticospinal neurons (Type II)",
  "Atypical excitatory neurons",
  "PV+ interneurons",
  "SST+ interneurons",
  "VIP+ interneurons",
  "Astrocytes",
  "Oligodendrocyte precursor cells",
  "Perineuronal oligodendrocytes",
  "Newly formed oligodendrocytes",
  "Myelinating oligodendrocytes",
  "Microglia"
)

cluster_colors <- c(
  "#6B5B95", "#45B8AC", "#955251", "#4E84C4", "#B565A7", "#88B04B", "#C3447A",
  "#009B77", "#EFC050", "#7B6888", "#7FCDCD", "#DD4124", "#5B5EA6", "#E07A5F",
  "#4BACC6", "#ea8a33", "#E8A0BF", "#9B2335", "#C17BAE"
)
names(cluster_colors) <- cluster_names

# Set up factor for ClusterName and recode Genotype
seu$ClusterName <- factor(
  seu$seurat_clusters,
  levels = 0:(length(cluster_names) - 1),
  labels = cluster_names
)

seu$Genotype_short <- recode(
  seu$Genotype,
  "PV-Cre/tdTom" = "Ctrl",
  "PV-Cre/tdTom/Dnmt1 loxP2" = "KO"
)

# Build and filter dataframe: only Ctrl and not in exclude list
plotdf <- FetchData(seu, vars = c("ClusterName", "Genotype_short", "Dnmt1")) %>%
  filter(!is.na(ClusterName)) %>%
  filter(!ClusterName %in% exclude_clusters) %>%
  filter(Genotype_short == "Ctrl")
plotdf$Dnmt1_zscore <- as.numeric(scale(plotdf$Dnmt1))

# Only keep remaining clusters for color vector
kept_clusters <- setdiff(cluster_names, exclude_clusters)
kept_colors <- cluster_colors[kept_clusters]

ggplot(plotdf, aes(x = ClusterName, y = Dnmt1_zscore, fill = ClusterName)) +
  geom_violin(width = 0.6, scale = 'width', trim = TRUE, linewidth = 0.3) +
  scale_fill_manual(values = kept_colors) +
  geom_jitter(
    width = 0.15,
    size = 0.4, alpha = 0.6, color = "black", shape = 16
  ) +
  xlab("") +
  ylab("Z-score (across all cells)") +
  ggtitle("Dnmt1") +
  theme(
    axis.text.x = element_text(angle = 75, hjust = 1, vjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", size = 0.8),
    panel.background = element_blank(),
    plot.background = element_blank(),
    legend.position = "none"
  )




#### UMAP and gene expression -------
library(cowplot)

p2 <- FeaturePlot(
  object = pv_cortical_filtered_pnn,
  features = "Dnmt3a",
  reduction = "umap",
  cols = c("ghostwhite", "#b2182b"),
  pt.size = 0.1
)
p3 <- FeaturePlot(
  object = pv_cortical_filtered_pnn,
  features = "Ffar3",
  reduction = "umap",
  cols = c("ghostwhite", "#b2182b"),
  pt.size = 0.1
)
p4 <- FeaturePlot(
  object = pv_cortical_filtered_pnn,
  features = "Slc16a1",
  reduction = "umap",
  cols = c("ghostwhite", "#b2182b"),
  pt.size = 0.1
)
p5 <- FeaturePlot(
  object = pv_cortical_filtered_pnn,
  features = "Slc16a7",
  reduction = "umap",
  cols = c("ghostwhite", "#b2182b"),
  pt.size = 0.1
)

# Align and arrange plots
plot_grid(p2, p3, p4, p5, ncol = 2, align = "hv", axis = "tblr")
final_plot <- plot_grid(p2, p3, p4, p5, ncol = 2, align = "hv", axis = "tblr")
ggsave("aligned_plots.pdf", final_plot, width = 9, height = 6.7)





############################### CORTICAL DATA ########################### ------------
##### Set everything up and load the cortical data --------
setwd("~/Downloads/Seurat_CellChat/revision/cellchat_new")
pv_cortical_filtered <- readRDS("~/Downloads/Seurat_CellChat/pv_cortical_subclustered_final.rds")

##### Quick reorder ------
pv_cortical_filtered <- RenameIdents(
  pv_cortical_filtered,
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
  "Myelinating oligodendrocytes",
  "Microglia",
  "Endothelial cells",
  "Leptomeningeal cells",
  "Meningeal fibroblasts"
)

pv_cortical_filtered@active.ident <- factor(pv_cortical_filtered@active.ident, levels = new_order)
levels(Idents(pv_cortical_filtered))

#################
##### QC METRICS TABLE (add after pv_cortical_filtered is defined) ----
library(dplyr)

# Extract metadata
meta <- pv_cortical_filtered@meta.data

# Assumes genotype is stored in a column called "Genotype"
qc_summary <- meta %>%
  group_by(Genotype) %>%
  summarise(
    n_nuclei        = n(),
    median_nCount   = median(nCount_RNA, na.rm = TRUE),
    mean_nCount     = mean(nCount_RNA, na.rm = TRUE),
    median_nFeature = median(nFeature_RNA, na.rm = TRUE),
    mean_nFeature   = mean(nFeature_RNA, na.rm = TRUE),
    .groups = "drop"
  )

meta$active_ident <- as.character(Idents(pv_cortical_filtered))

# Per-cluster cell counts split by genotype
cluster_counts <- meta %>%
  group_by(Genotype, active_ident) %>%   # use whichever column holds final cluster labels
  summarise(n_nuclei = n(), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = Genotype, values_from = n_nuclei, values_fill = 0)

write.csv(qc_summary,      "QC_metrics_per_genotype.csv",    row.names = FALSE)
write.csv(cluster_counts,  "QC_nuclei_per_cluster.csv",      row.names = FALSE)

print(qc_summary)
print(cluster_counts)

##### DEG: all clusters, one CSV each, saved to ./DEG/ ----

# Create output folder if it doesn't exist
dir.create("DEG", showWarnings = FALSE)

# Get all cluster names from the object
all_clusters_source <- levels(Idents(pv_cortical_filtered_pnn))

cat("Found", length(all_clusters_source), "clusters. Starting export...\n\n")

for (clust in all_clusters_source) {
  
  cells <- WhichCells(pv_cortical_filtered_pnn, idents = clust)
  
  # Skip clusters with too few cells to test
  if (length(cells) < 10) {
    cat("Skipping cluster (< 10 cells):", clust, "\n")
    next
  }
  
  cat("Processing:", clust, paste0("(", length(cells), " cells)"), "\n")
  
  subobj <- subset(pv_cortical_filtered_pnn, cells = cells)
  Idents(subobj) <- "Genotype"
  
  # Check both genotypes are represented
  geno_counts <- table(Idents(subobj))
  if (length(geno_counts) < 2 || any(geno_counts < 3)) {
    cat("  Skipping — one genotype has fewer than 3 cells.\n")
    next
  }
  
  de_full <- tryCatch({
    FindMarkers(
      object          = subobj,
      ident.1         = "PV-Cre/tdTom/Dnmt1 loxP2",
      ident.2         = "PV-Cre/tdTom",
      min.pct         = 0.1,
      logfc.threshold = 0,
      test.use        = "wilcox"
    )
  }, error = function(e) {
    cat("  ERROR in FindMarkers:", conditionMessage(e), "\n")
    return(NULL)
  })
  
  if (is.null(de_full) || nrow(de_full) == 0) {
    cat("  No results returned, skipping.\n")
    next
  }
  
  de_full$gene    <- rownames(de_full)
  de_full$cluster <- clust
  de_full         <- de_full[!grepl("^ENSMUSG", de_full$gene), ]
  de_full         <- de_full[order(de_full$p_val_adj), ]
  
  # Reorder columns for readability
  de_full <- de_full[, c("cluster", "gene", "avg_log2FC", "p_val", "p_val_adj",
                         "pct.1", "pct.2")]
  
  # Sanitize cluster name for filename
  safe_name <- gsub("[^A-Za-z0-9_-]", "_", clust)
  safe_name <- gsub("_{2,}", "_", safe_name)
  safe_name <- gsub("_$", "", safe_name)
  
  out_path <- file.path("DEG", paste0("SourceData_DEG_", safe_name, ".csv"))
  write.csv(de_full, file = out_path, row.names = FALSE)
  
  cat("  Saved:", out_path, paste0("(", nrow(de_full), " genes)\n"))
}

cat("\nDone. All files saved to ./DEG/\n")



#################
# Subset by Genotype -----------------------------------------------------------
# Check exact genotype names
table(pv_cortical_filtered$Genotype)

# Create subsets
wt_cortical <- subset(pv_cortical_filtered, subset = Genotype == "PV-Cre/tdTom")
ko_cortical <- subset(pv_cortical_filtered, subset = Genotype == "PV-Cre/tdTom/Dnmt1 loxP2")

# Validate subsets and check cell counts
dim(wt_cortical)  
dim(ko_cortical)

# For wt object
table(Idents(wt_cortical))

# For ko object
table(Idents(ko_cortical))

# For wt
wt_cortical_cluster_counts <- as.data.frame(table(Idents(wt_cortical)))
colnames(wt_cortical_cluster_counts) <- c("Cluster", "CellNumber")
print(wt_cortical_cluster_counts)

# For ko
ko_cortical_cluster_counts <- as.data.frame(table(Idents(ko_cortical)))
colnames(ko_cortical_cluster_counts) <- c("Cluster", "CellNumber")
print(ko_cortical_cluster_counts)

# For wt
ggplot(wt_cortical_cluster_counts, aes(x=Cluster, y=CellNumber)) +
  geom_bar(stat="identity") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Cell Numbers per Cluster (WT)")


# For ko
ggplot(ko_cortical_cluster_counts, aes(x=Cluster, y=CellNumber)) +
  geom_bar(stat="identity") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Cell Numbers per Cluster (KO)")


# Proportional to total WT and KO cells---------------------

# 1. Get cluster counts for WT and KO
wt_cortical_cluster_counts <- as.data.frame(table(Idents(wt_cortical)))
colnames(wt_cortical_cluster_counts) <- c("Cluster", "CellNumber_WT")

ko_cortical_cluster_counts <- as.data.frame(table(Idents(ko_cortical)))
colnames(ko_cortical_cluster_counts) <- c("Cluster", "CellNumber_KO")

# 2. Get total number of cells in each dataset
total_wt_cortical <- ncol(wt_cortical)  # or dim(wt)[1]
total_ko_cortical <- ncol(ko_cortical)  # or dim(ko)[1]

# 3. Calculate proportion of each cluster in WT and KO
wt_cortical_cluster_counts$Proportion_WT <- wt_cortical_cluster_counts$CellNumber_WT / total_wt_cortical
ko_cortical_cluster_counts$Proportion_KO <- ko_cortical_cluster_counts$CellNumber_KO / total_ko_cortical

# 4. Merge the tables
merged_cortical_counts <- merge(wt_cortical_cluster_counts, ko_cortical_cluster_counts, by = "Cluster", all = TRUE)
merged_cortical_counts[is.na(merged_cortical_counts)] <- 0  # Replace NAs with 0

# 5. Calculate relative difference in proportions (as percent)
merged_cortical_counts$RelativeDifference <- (merged_cortical_counts$Proportion_KO - merged_cortical_counts$Proportion_WT) / merged_cortical_counts$Proportion_WT * 100

# 6. Handle Inf/NaN (if Proportion_WT is zero)
merged_cortical_counts$RelativeDifference[!is.finite(merged_cortical_counts$RelativeDifference)] <- NA

merged_cortical_counts$Proportion_WT <- round(merged_cortical_counts$Proportion_WT * 100, 2)
merged_cortical_counts$Proportion_KO <- round(merged_cortical_counts$Proportion_KO * 100, 2)
merged_cortical_counts$RelativeDifference <- round(merged_cortical_counts$RelativeDifference, 2)
write.csv(merged_cortical_counts, file = "cluster_relative_difference_WT_KO_cortical.csv", row.names = FALSE)

# Assuming merged_counts from your previous code
merged_cortical_counts$p_value <- NA

for (i in 1:nrow(merged_cortical_counts)) {
  k_WT <- merged_cortical_counts$CellNumber_WT[i]
  k_KO <- merged_cortical_counts$CellNumber_KO[i]
  # Only test if both groups have at least one cell
  if ((k_WT + k_KO) > 0) {
    test <- prop.test(
      x = c(k_WT, k_KO),
      n = c(total_wt_cortical, total_ko_cortical),
      alternative = "two.sided"
    )
    merged_cortical_counts$p_value[i] <- test$p.value
  }
}
merged_cortical_counts$adj_p_value <- p.adjust(merged_cortical_counts$p_value, method = "fdr")
write.csv(merged_cortical_counts, file = "cluster_relative_difference_with_pvalues_cortical.csv", row.names = FALSE)


# Side by side proportion analysis for all cell types ------------------------


# Get cluster counts and calculate proportions
wt_cortical_prop <- as.data.frame(table(Idents(wt_cortical)))
colnames(wt_cortical_prop) <- c("Cluster", "Count_WT")
wt_cortical_prop$Proportion_WT <- wt_cortical_prop$Count_WT / ncol(wt_cortical) * 100

ko_cortical_prop <- as.data.frame(table(Idents(ko_cortical)))
colnames(ko_cortical_prop) <- c("Cluster", "Count_KO")
ko_cortical_prop$Proportion_KO <- ko_cortical_prop$Count_KO / ncol(ko_cortical) * 100

# Merge data
merged_cortical_prop <- merge(wt_cortical_prop, ko_cortical_prop, by = "Cluster", all = TRUE)
merged_cortical_prop[is.na(merged_cortical_prop)] <- 0

# Convert to long format for plotting
plot_cortical_data <- merged_cortical_prop %>% 
  pivot_longer(
    cols = c(Proportion_WT, Proportion_KO),
    names_to = "Condition",
    values_to = "Proportion"
  ) %>% 
  mutate(Condition = factor(
    Condition, 
    levels = c("Proportion_WT", "Proportion_KO"),
    labels = c("WT", "KO")
  ))

# Build plot: remove grid, keep axis lines, customize text
p <- ggplot(plot_cortical_data, aes(x = Cluster, y = Proportion, fill = Condition)) +
  geom_col(position = position_dodge(0.9), width = 0.8) +
  scale_fill_manual(values = c("WT" = "#c6c6c6", "KO" = "#7FCDCD")) +
  labs(
    title = "Number of cells per cluster",
    x = "Cluster",
    y = "Proportion of Total Cells (%)"
  ) +
  # theme_classic removes gridlines and shows axis lines
  theme_classic(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 70, hjust = 1),
    axis.title.x = element_text(margin = margin(t = 20)),
    legend.position = "right",
    # Explicitly ensure grid is gone and axis lines are present
    panel.grid = element_blank(),
    axis.line = element_line(color = "black")
  )

print(p)

# Save 
ggsave(
  filename = "cortical_cluster_proportions.pdf",
  plot = p,
  width = 8,
  height = 5,
  units = "in",
  device = "pdf"
)



#################
# Select gene dot blot marker with adjustable color palette ------
# Libraries
library(Seurat)
library(dplyr)
library(ggplot2)

# Curated mouse-case gene panels (MGI style: first letter uppercase, rest lowercase)
# 20-gene panel (replace Gfap -> Scnn1a)
genes_20_mouse <- c(
  # Pan-neuronal, excitatory, inhibitory
  "Snap25","Slc17a7","Gad1",
  # Excitatory / laminar
  "Cux2","Rorb","Fezf2","Tle4","Scnn1a", "Drd1", "Slc17a6", "Tox", "Penk",
  # Interneurons
  "Pvalb","Sst","Vip", "Reln",
  # Oligodendrocyte lineage
  "Pdgfra","Enpp6","Mbp",
  # Astrocytes
  "Aqp4",
  # Microglia
  "Tmem119",
  # Endothelial
  "Pecam1",
  # Meningeal
  "Foxc1","Dcn","Col1a1"
)

# 33-gene panel (30 original + 3 new)
genes_30_mouse <- c(
  # Pan-neuronal, excitatory, inhibitory
  "Snap25","Slc17a7","Gad1",
  # Excitatory / laminar
  "Cux2","Satb2","Rorb","Scnn1a","Pcp4","Fezf2","Bcl11b","Tle4","Foxp2",
  # Interneurons
  "Pvalb","Sst","Vip","Npy", "Reln",
  # Oligodendrocyte lineage
  "Pdgfra","Cspg4","Enpp6","Mbp","Mog","Plp1",
  # Astrocytes
  "Aqp4","Gfap","Slc1a2",
  # Microglia
  "P2ry12","Tmem119","Cx3cr1",
  # Endothelial
  "Pecam1","Vwf","Flt1",
  # Meningeal
  "Foxc1","Col1a1"
)

# Robust case-insensitive matching to features present in the Seurat object
all_features <- rownames(pv_cortical_filtered[["RNA"]])

match_to_features <- function(gvec, features) {
  lut <- setNames(features, toupper(features))         # map UPPER -> original
  matched <- lut[toupper(gvec)]
  unique(unname(matched[!is.na(matched)]))
}

features_20 <- match_to_features(genes_20_mouse, all_features)
features_30 <- match_to_features(genes_30_mouse, all_features)

# Optional: report any genes not found
missing_report <- function(requested, matched) {
  req_up <- toupper(requested)
  feat_up <- toupper(all_features)
  missing <- requested[!(req_up %in% feat_up)]
  if (length(missing) > 0) message("Not found: ", paste(missing, collapse = ", "))
}
missing_report(genes_20_mouse, features_20)
missing_report(genes_30_mouse, features_30)

# Define a blue gradient centered on #2166ac
# Light -> medium -> dark blues (ColorBrewer-like), including #2166ac
pal_blues <- colorRampPalette(c(
  "#f7fbff", "#deebf7", "#c6dbef", "#9ecae1",
  "#6baed6", "#4292c6", "#2171b5", "#2166ac", "#084594"
))(100)

# Keep your clipping consistent across DotPlot and the scale
cmin <- -1.5
cmax <- 2.5

# DotPlot: 20-gene panel
p20 <- DotPlot(
  pv_cortical_filtered,
  features = features_20,
  dot.scale = 5.5,
  col.min = cmin,
  col.max = cmax
) +
  theme(axis.text.x = element_text(angle = 45, #vjust = 0.5, 
                                   hjust = 1, face = "italic")) +
  labs(title = "Marker panel") +
  scale_colour_gradientn(
    colours = pal_blues,
    limits = c(cmin, cmax),
    oob = scales::squish,
    name = "Avg exp (scaled)"
  )

# DotPlot: 30-gene panel
p30 <- DotPlot(
  pv_cortical_filtered,
  features = features_30,
  dot.scale = 5.5,
  col.min = cmin,
  col.max = cmax
) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, face = "italic")) +
  labs(title = "Marker panel") +
  scale_colour_gradientn(
    colours = pal_blues,
    limits = c(cmin, cmax),
    oob = scales::squish,
    name = "Avg exp (scaled)"
  )

# Render or save
print(p20)
#print(p30)

ggsave("dotplot_20.pdf", p20, width = 12, height = 7)
#ggsave("dotplot_30.png", p30, width = 12, height = 6, dpi = 600)


#################
# Oligo population -----
# Load required libraries
library(Seurat)
library(ggplot2)
library(dplyr)
library(scales)

# Assumes objects: wt_cortical, ko_cortical, total_wt_cortical, total_ko_cortical

# Cells in each group
wt_nmo <- WhichCells(wt_cortical, idents = "Newly formed oligodendrocytes")
wt_myo <- WhichCells(wt_cortical, idents = "Myelinating oligodendrocytes")
ko_nmo <- WhichCells(ko_cortical, idents = "Newly formed oligodendrocytes")
ko_myo <- WhichCells(ko_cortical, idents = "Myelinating oligodendrocytes")

# Union for mature
wt_mature <- union(wt_nmo, wt_myo)
ko_mature <- union(ko_nmo, ko_myo)

# Counts
n_wt_nmo <- length(wt_nmo)
n_wt_myo <- length(wt_myo)
n_wt_mature <- length(wt_mature)

n_ko_nmo <- length(ko_nmo)
n_ko_myo <- length(ko_myo)
n_ko_mature <- length(ko_mature)

# Proportions
prop_wt_nmo <- n_wt_nmo / total_wt_cortical
prop_wt_myo <- n_wt_myo / total_wt_cortical
prop_wt_mature <- n_wt_mature / total_wt_cortical

prop_ko_nmo <- n_ko_nmo / total_ko_cortical
prop_ko_myo <- n_ko_myo / total_ko_cortical
prop_ko_mature <- n_ko_mature / total_ko_cortical

# Proportion tests
test_nmo <- prop.test(c(n_wt_nmo, n_ko_nmo), c(total_wt_cortical, total_ko_cortical))
test_myo <- prop.test(c(n_wt_myo, n_ko_myo), c(total_wt_cortical, total_ko_cortical))
test_mature <- prop.test(c(n_wt_mature, n_ko_mature), c(total_wt_cortical, total_ko_cortical))

# Data frame for plotting
plot_data <- data.frame(
  Group = rep(c("NFOL", "MOL", "Mature oligodendrocytes"), each = 2),
  Condition = factor(rep(c("WT", "KO"), 3), levels = c("WT", "KO")),
  Proportion = c(prop_wt_nmo, prop_ko_nmo, prop_wt_myo, prop_ko_myo, prop_wt_mature, prop_ko_mature),
  N_Pos = c(n_wt_nmo, n_ko_nmo, n_wt_myo, n_ko_myo, n_wt_mature, n_ko_mature)
)

# Significance
plot_data$P_Value <- rep(c(test_nmo$p.value, test_myo$p.value, test_mature$p.value), each = 2)
plot_data$Significance <- cut(plot_data$P_Value, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf), 
                              labels = c("***", "**", "*", "ns"))

# Explicit plotting order
plot_data$Group <- factor(plot_data$Group, levels = c("NFOL", "MOL", "Mature oligodendrocytes"))

# Annotation positions
annotation_data <- plot_data %>% 
  group_by(Group) %>% 
  summarise(
    y_pos = max(Proportion) + 0.002, 
    p_value = P_Value[1],
    significance = Significance[1], .groups = 'drop'
  )


# Pdgfc expressing cells - Perineuronal oligodendrocytes normalized to all cells --------------------
# Load necessary libraries
library(Seurat)
library(ggplot2)
library(dplyr)
library(scales)

# Assumes 'wt_cortical', 'ko_cortical', 'total_wt_cortical', 'total_ko_cortical' exist

gene_name <- "Pdgfc"

# Expression matrices
expr_mat_wt <- GetAssayData(wt_cortical, layer = "data")
expr_mat_ko <- GetAssayData(ko_cortical, layer = "data")

# Counts of Pdgfc+ cells in WT
n_wt_OPC <- sum(expr_mat_wt[gene_name, WhichCells(wt_cortical, idents = "Oligodendrocyte precursor cells")] > 0)
n_wt_NMO <- sum(expr_mat_wt[gene_name, WhichCells(wt_cortical, idents = "Newly formed oligodendrocytes")] > 0)
n_wt_MYO <- sum(expr_mat_wt[gene_name, WhichCells(wt_cortical, idents = "Myelinating oligodendrocytes")] > 0)

# Counts of Pdgfc+ cells in KO
n_ko_OPC <- sum(expr_mat_ko[gene_name, WhichCells(ko_cortical, idents = "Oligodendrocyte precursor cells")] > 0)
n_ko_NMO <- sum(expr_mat_ko[gene_name, WhichCells(ko_cortical, idents = "Newly formed oligodendrocytes")] > 0)
n_ko_MYO <- sum(expr_mat_ko[gene_name, WhichCells(ko_cortical, idents = "Myelinating oligodendrocytes")] > 0)

# Proportions normalized to all cells
prop_wt_OPC <- n_wt_OPC / total_wt_cortical
prop_wt_NMO <- n_wt_NMO / total_wt_cortical
prop_wt_MYO <- n_wt_MYO / total_wt_cortical

prop_ko_OPC <- n_ko_OPC / total_ko_cortical
prop_ko_NMO <- n_ko_NMO / total_ko_cortical
prop_ko_MYO <- n_ko_MYO / total_ko_cortical

# Plotting data
plot_data <- data.frame(
  Cell_Type = rep(c("OPC", "NFOL", "MOL"), each = 2),
  Condition = rep(c("WT", "KO"), 3),
  Proportion = c(prop_wt_OPC, prop_ko_OPC,
                 prop_wt_NMO, prop_ko_NMO, 
                 prop_wt_MYO, prop_ko_MYO)
)

# Enforce display order
plot_data$Cell_Type <- factor(plot_data$Cell_Type, levels = c("OPC", "NFOL", "MOL"))
plot_data$Condition <- factor(plot_data$Condition, levels = c("WT", "KO"))

# Tests
opc_test <- prop.test(c(n_wt_OPC, n_ko_OPC), c(total_wt_cortical, total_ko_cortical))
nmo_test <- prop.test(c(n_wt_NMO, n_ko_NMO), c(total_wt_cortical, total_ko_cortical))
myo_test <- prop.test(c(n_wt_MYO, n_ko_MYO), c(total_wt_cortical, total_ko_cortical))

p_values <- c(opc_test$p.value, nmo_test$p.value, myo_test$p.value)
p_values_adj <- p.adjust(p_values, method = "fdr")

annotation_data <- plot_data |>
  dplyr::group_by(Cell_Type) |>
  dplyr::summarise(y_pos = max(Proportion) + 0.002, .groups = "drop") |>
  dplyr::mutate(
    p_value = p_values_adj,
    significance = dplyr::case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE            ~ "ns"
    )
  )

# Plot with small space between paired bars
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
    aes(x = Cell_Type, y = y_pos + 0.001, label = paste("p =", format.pval(p_value, digits = 2))),
    inherit.aes = FALSE, hjust = 0.5, vjust = 1, size = 4
  ) +
  # geom_text(
  #   data = annotation_data, 
  #   aes(x = Cell_Type, y = y_pos + 0.015, label = significance),
  #   inherit.aes = FALSE, hjust = 0.5, vjust = 0, size = 6, fontface = "bold"
  # ) +
  labs(
    x = "Cell Type", 
    y = "Proportion of Pdgfc+ cells (of total)", 
    fill = "Genotype"
  ) +
  scale_x_discrete(labels = c(
    "OPC",
    "NFOL", 
    "MOL"
  )) +
  scale_y_continuous(labels = label_percent(accuracy = 0.1),
                     expand = expansion(mult = c(0, 0.05)),
                     limits = c(0, 0.016)) +
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
  "Pdgfc-population_OL.pdf",
  plot = p1,
  width = 6,
  height = 6,
  dpi = 600,
  bg = "transparent"
)

################
# ---- DNMT1 EXPRESSION IN PV CLUSTER: violin plot with p-value and log2FC ----

library(Seurat)
library(ggplot2)
library(dplyr)

# Subset to PV+ interneurons only
pv_cells <- subset(pv_cortical_filtered_pnn, idents = "PV+ interneurons")
pv_cells$Genotype_short <- recode(pv_cells$Genotype,
                                  "PV-Cre/tdTom"              = "Ctrl",
                                  "PV-Cre/tdTom/Dnmt1 loxP2" = "KO"
)
Idents(pv_cells) <- "Genotype_short"

# Run FindMarkers for Dnmt1 only, no thresholds
dnmt1_stats <- FindMarkers(
  object          = pv_cells,
  ident.1         = "KO",
  ident.2         = "Ctrl",
  features        = "Dnmt1",
  min.pct         = 0,
  logfc.threshold = 0,
  test.use        = "wilcox"
)

# Extract values for annotation
log2fc  <- round(dnmt1_stats$avg_log2FC, 3)
pval    <- dnmt1_stats$p_val_adj
pval_fmt <- ifelse(pval < 0.001,
                   formatC(pval, format = "e", digits = 2),
                   round(pval, 4))

cat("Dnmt1 in PV+ interneurons — log2FC (KO/Ctrl):", log2fc,
    "| adj. p-value:", pval_fmt, "\n")

# Build annotation label for the plot
annot_label <- paste0(
  "\u0394log2FC = ", log2fc,
  "\nadj. p = ", pval_fmt
)

# Get y position for annotation (just above the highest violin)
plot_data <- FetchData(pv_cells, vars = c("Dnmt1", "Genotype_short"))
y_max     <- max(plot_data$Dnmt1, na.rm = TRUE)
y_annot   <- y_max * 1.05

# Violin plot
p_dnmt1 <- VlnPlot(
  pv_cells,
  features     = "Dnmt1",
  cols         = c("Ctrl" = "#c6c6c6", "KO" = "#7FCDCD"),
  pt.size      = 0.8,
  log          = FALSE
) +
  # Significance bar
  annotate("segment",
           x = 1, xend = 2,
           y = y_annot * 0.98, yend = y_annot * 0.98,
           linewidth = 0.5, color = "black") +
  
  # Annotation text
  annotate("text",
           x     = 1.5,
           y     = y_annot * 1.08,
           label = annot_label,
           size  = 3.5,
           hjust = 0.5,
           vjust = 0,
           color = "black") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) +
  labs(
    title    = "Dnmt1 expression in PV+ interneurons",
    x        = NULL,
    y        = "Normalized expression"
  ) +
  theme(
    legend.position  = "none",
    plot.title       = element_text(hjust = 0.5, size = 13, face = "bold"),
    plot.subtitle    = element_text(hjust = 0.5, size = 10, color = "grey40"),
    axis.text.x      = element_text(size = 12),
    panel.grid       = element_blank(),
    axis.line        = element_line(color = "black", linewidth = 0.5)
  )

print(p_dnmt1)

# Save
ggsave(
  filename ="Dnmt1_VlnPlot_PV_Ctrl_vs_KO.pdf",
  plot     = p_dnmt1,
  width    = 4, 
  height   = 5.5,
  bg       = "white"
)

# Also save the raw stats as CSV
write.csv(dnmt1_stats,
          file      = "Dnmt1_KO_vs_WT_PVcluster_stats.csv",
          row.names = TRUE)

cat("Saved plot and stats to ./revision/\n")

# ---- DNMT1 KNOCKOUT VALIDATION: proportion of expressing nuclei ----

pv_cells <- subset(pv_cortical_filtered_pnn, idents = "PV+ interneurons")
pv_cells$Genotype_short <- recode(pv_cells$Genotype,
                                  "PV-Cre/tdTom"              = "Ctrl",
                                  "PV-Cre/tdTom/Dnmt1 loxP2" = "KO"
)

expr_mat <- GetAssayData(pv_cells, layer = "data")

ctrl_cells <- WhichCells(pv_cells, expression = Genotype_short == "Ctrl")
ko_cells   <- WhichCells(pv_cells, expression = Genotype_short == "KO")

# Count expressing nuclei (counts > 0)
n_ctrl_expressing <- sum(expr_mat["Dnmt1", ctrl_cells] > 0)
n_ko_expressing   <- sum(expr_mat["Dnmt1", ko_cells]   > 0)
n_ctrl_total      <- length(ctrl_cells)
n_ko_total        <- length(ko_cells)

prop_ctrl <- n_ctrl_expressing / n_ctrl_total
prop_ko   <- n_ko_expressing   / n_ko_total

cat("Ctrl: ", n_ctrl_expressing, "/", n_ctrl_total,
    "(", round(prop_ctrl * 100, 1), "% expressing)\n")
cat("KO:   ", n_ko_expressing,   "/", n_ko_total,
    "(", round(prop_ko   * 100, 1), "% expressing)\n")

# Fisher's exact test (better than prop.test for sparse counts)
cont_table <- matrix(
  c(n_ctrl_expressing, n_ctrl_total - n_ctrl_expressing,
    n_ko_expressing,   n_ko_total   - n_ko_expressing),
  nrow = 2,
  dimnames = list(
    Genotype   = c("Ctrl", "KO"),
    Expressing = c("Yes",  "No")
  )
)
fisher_result <- fisher.test(cont_table)
pval_fmt <- ifelse(fisher_result$p.value < 0.001,
                   formatC(fisher_result$p.value, format = "e", digits = 2),
                   round(fisher_result$p.value, 4))

cat("Fisher's exact test p-value:", pval_fmt, "\n")
cat("Odds ratio:", round(fisher_result$estimate, 3), "\n")

# PLOT: bar chart of % expressing nuclei
plot_df <- data.frame(
  Genotype   = factor(c("Ctrl", "KO"), levels = c("Ctrl", "KO")),
  Proportion = c(prop_ctrl, prop_ko) * 100,
  Total      = c(n_ctrl_total, n_ko_total),
  Expressing = c(n_ctrl_expressing, n_ko_expressing)
)

# annotation position
y_annot <- max(plot_df$Proportion) * 1.15

annot_label <- paste0(
  "Fisher's p = ", pval_fmt,
  "\nlog\u2082FC = ", round(log2(prop_ko / prop_ctrl), 3)
)

p_dnmt1 <- ggplot(plot_df, aes(x = Genotype, y = Proportion, fill = Genotype)) +
  geom_bar(stat = "identity", width = 0.55, alpha = 0.85) +
  geom_text(aes(label = paste0(Expressing, "/", Total)),
            vjust = -0.5, size = 3.5, color = "black") +
  scale_fill_manual(values = c("Ctrl" = "#c6c6c6", "KO" = "#7FCDCD")) +
  # bracket
  annotate("segment",
           x = 1, xend = 2,
           y = y_annot * 0.95, yend = y_annot * 0.95,
           linewidth = 0.5, color = "black") +
  annotate("text",
           x = 1.5, y = y_annot * 1.02,
           label = annot_label,
           size = 3.5, hjust = 0.5, color = "black") +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.35)),
    labels = function(x) paste0(x, "%")
  ) +
  labs(
    title    = "Dnmt1-expressing nuclei in PV+ interneurons",
    subtitle = "Proportion of nuclei with detectable Dnmt1 transcripts",
    x        = NULL,
    y        = "% Dnmt1+ nuclei"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    plot.title      = element_text(hjust = 0.5, size = 13, face = "bold"),
    plot.subtitle   = element_text(hjust = 0.5, size = 10, color = "grey40"),
    axis.text.x     = element_text(size = 12),
    axis.line       = element_line(color = "black", linewidth = 0.5)
  )

print(p_dnmt1)

ggsave(
  filename = "Dnmt1_PV_expressing_nuclei_proportion.pdf",
  p_dnmt1, width = 4, height = 5.5, dpi = 300, bg = "white"
)

# Save stats
stats_out <- data.frame(
  metric            = c("n_Ctrl_total", "n_KO_total",
                        "n_Ctrl_expressing", "n_KO_expressing",
                        "pct_Ctrl_expressing", "pct_KO_expressing",
                        "log2FC_proportion", "Fisher_pvalue"),
  value             = c(n_ctrl_total, n_ko_total,
                        n_ctrl_expressing, n_ko_expressing,
                        round(prop_ctrl * 100, 2), round(prop_ko * 100, 2),
                        round(log2(prop_ko / prop_ctrl), 3),
                        fisher_result$p.value)
)
write.csv(stats_out,
          "Dnmt1_PV_proportion_stats.csv",
          row.names = FALSE)


# that were upregulated in KO PV cells (endocytosis/vesicle genes)
known_targets <- c("Dnm1", "Clta", "Cltb", "Ap2a1", "Syt1")  # example - use your actual genes

pv_cells <- subset(pv_cortical_filtered_pnn, idents = "PV+ interneurons")
Idents(pv_cells) <- "Genotype"

target_de <- FindMarkers(
  object          = pv_cells,
  ident.1         = "PV-Cre/tdTom/Dnmt1 loxP2",
  ident.2         = "PV-Cre/tdTom",
  features        = known_targets,
  min.pct         = 0,
  logfc.threshold = 0,
  test.use        = "wilcox"
)
target_de$gene <- rownames(target_de)
print(target_de)
write.csv(target_de,
          "KnownDnmt1Targets_PV_KO_vs_Ctrl.csv",
          row.names = FALSE)

################
# ---- PV SUBTYPE DOTPLOTS (matching existing color palette) ----
library(Seurat)
library(ggplot2)
library(dplyr)

# ---- Your exact palette ----
pal_blues <- colorRampPalette(c(
  "#f7fbff", "#deebf7", "#c6dbef", "#9ecae1",
  "#6baed6", "#4292c6", "#2171b5", "#2166ac", "#084594"
))(100)

cmin <- -1.5
cmax <-  2.5

# ---- Subset and recluster PV interneurons ----
pv_only <- subset(pv_cortical_filtered_pnn, idents = "PV+ interneurons")

pv_only <- NormalizeData(pv_only)
pv_only <- FindVariableFeatures(pv_only, nfeatures = 2000)
pv_only <- ScaleData(pv_only)
pv_only <- RunPCA(pv_only, npcs = 20)
pv_only <- FindNeighbors(pv_only, dims = 1:10)
pv_only <- FindClusters(pv_only, resolution = 0.3)
pv_only <- RunUMAP(pv_only, dims = 1:10)

cat("PV sub-clusters found:\n")
print(table(Idents(pv_only)))

# ---- Case-insensitive gene matching (your exact helper) ----
all_features_pv <- rownames(pv_only[["RNA"]])

match_to_features <- function(gvec, features) {
  lut     <- setNames(features, toupper(features))
  matched <- lut[toupper(gvec)]
  unique(unname(matched[!is.na(matched)]))
}

missing_report <- function(requested, features) {
  missing <- requested[!(toupper(requested) %in% toupper(features))]
  if (length(missing) > 0) message("Not found: ", paste(missing, collapse = ", "))
}

# ---- Find data-driven top markers per sub-cluster ----
pv_subtype_markers <- FindAllMarkers(
  pv_only,
  min.pct         = 0.1,
  logfc.threshold = 0.25,
  only.pos        = TRUE,
  test.use        = "wilcox"
)

top5_pv <- pv_subtype_markers %>%
  filter(p_val_adj < 0.05) %>%
  filter(!grepl("^ENSMUSG", gene)) %>%
  filter(!grepl("^Gm",      gene)) %>%
  filter(!grepl("Rik$",     gene)) %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 5) %>%
  ungroup()

top5_genes_pv <- unique(top5_pv$gene)

write.csv(pv_subtype_markers,
          "PV_subtype_marker_genes.csv",
          row.names = FALSE)

cat("Top 5 marker genes per PV sub-cluster:\n")
print(top5_pv %>% dplyr::select(cluster, gene, avg_log2FC, p_val_adj))

features_datadriven <- match_to_features(top5_genes_pv, all_features_pv)
missing_report(top5_genes_pv, all_features_pv)

# ---- DOTPLOT 1: data-driven top marker genes ----
p_dot_markers <- DotPlot(
  pv_only,
  features  = features_datadriven,
  dot.scale = 5,
  col.min   = cmin,
  col.max   = cmax
) +
  scale_colour_gradientn(
    colours = pal_blues,
    limits  = c(cmin, cmax),
    oob     = scales::squish,
    name    = "Avg exp (scaled)"
  ) +
  labs(
    title = "PV interneuron sub-cluster: data-driven marker genes",
    x     = NULL,
    y     = "Sub-cluster"
  ) +
  theme_classic() +
  theme(
    axis.text.x     = element_text(angle = 90, vjust = 0.5, hjust = 1,
                                   face = "italic", size = 9),
    axis.text.y     = element_text(size = 10),
    plot.title      = element_text(hjust = 0.5, face = "bold", size = 13),
    legend.position = "right",
    panel.grid      = element_blank(),
    axis.line       = element_line(color = "black", linewidth = 0.4)
  )

print(p_dot_markers)

ggsave(
  filename = "DotPlot_PV_subtype_datadriven_markers.pdf",
  p_dot_markers,
  width  = max(8, length(features_datadriven) * 0.38 + 3),
  height = 5,
  dpi    = 300,
  bg     = "white"
)

# ---- DOTPLOT 2: canonical known PV subtype markers ----
# Chandelier (axo-axonic): Unc5b, Adarb2, Pthlh, Nkx2-1
# Basket cells:            Cntn1, Syt2, Mef2c, Kcnc1
# General PV identity:     Pvalb, Erbb4, Gad1, Gad2
# PNN / maturity:          Hapln1, Bcan, Tnr
# PV-SST overlap:          Sst, Cort

genes_pv_canonical <- c(
  # General PV identity
  "Pvalb", "Erbb4", "Gad1", "Gad2",
  # Chandelier markers
  "Unc5b", "Adarb2", "Pthlh", "Nkx2-1",
  # Basket cell markers
  "Cntn1", "Syt2", "Mef2c", "Kcnc1",
  # PNN / maturity markers
  "Hapln1", "Bcan", "Tnr",
  # PV-SST overlap
  "Sst", "Cort"
)

features_canonical <- match_to_features(genes_pv_canonical, all_features_pv)
missing_report(genes_pv_canonical, all_features_pv)

cat("Canonical markers found:", paste(features_canonical, collapse = ", "), "\n")

p_dot_known <- DotPlot(
  pv_only,
  features  = features_canonical,
  dot.scale = 6,
  col.min   = cmin,
  col.max   = cmax
) +
  scale_colour_gradientn(
    colours = pal_blues,
    limits  = c(cmin, cmax),
    oob     = scales::squish,
    name    = "Avg exp (scaled)"
  ) +
  # Vertical dividers between gene groups
  # Positions: after General (4), Chandelier (8), Basket (12), PNN (15)
  # Adjust if genes are missing from your object
  geom_vline(
    xintercept = c(4.5, 8.5, 12.5, 15.5),
    linetype   = "dashed",
    color      = "grey70",
    linewidth  = 0.3
  ) +
  labs(
    title    = "Canonical PV interneuron subtype markers",
    subtitle = "General identity | Chandelier | Basket | PNN/maturity | PV-SST",
    x        = NULL,
    y        = "Sub-cluster"
  ) +
  theme_classic() +
  theme(
    axis.text.x     = element_text(angle = 90, vjust = 0.5, hjust = 1,
                                   face = "italic", size = 10),
    axis.text.y     = element_text(size = 10),
    plot.title      = element_text(hjust = 0.5, face = "bold", size = 13),
    plot.subtitle   = element_text(hjust = 0.5, size = 9, color = "grey40"),
    legend.position = "right",
    panel.grid      = element_blank(),
    axis.line       = element_line(color = "black", linewidth = 0.4)
  )

print(p_dot_known)

ggsave(
  filename = "DotPlot_PV_subtype_canonical_markers.pdf",
  p_dot_known,
  width  = length(features_canonical) * 0.5 + 3,
  height = 5,
  dpi    = 300,
  bg     = "white"
)

# ---- UMAP for context ----
p_umap_pv <- DimPlot(pv_only, reduction = "umap", label = FALSE, repel = TRUE) +
  ggtitle("PV interneuron sub-clusters") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave(
  filename = "UMAP_PV_subclusters.pdf",
  p_umap_pv,
  width = 6, height = 5, dpi = 300, bg = "white"
)

###############
# PV SUBCLUSTER ANALYSIS — PART 1: SETUP + PROPORTION ANALYSIS-------

library(Seurat)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(clusterProfiler)
library(org.Mm.eg.db)
library(ReactomePA)
library(DOSE)
library(cowplot)
library(RColorBrewer)

select <- dplyr::select
filter <- dplyr::filter

# Directory structure
base_dir       <- "pv_deg"
deg_dir        <- file.path(base_dir, "DEG_results")
plot_top10_dir <- file.path(deg_dir,  "volcano_top10_labeled")
plot_sig_dir   <- file.path(deg_dir,  "volcano_significant_labeled")
plot_clean_dir <- file.path(deg_dir,  "volcano_clean")
go_dir         <- file.path(base_dir, "GO_results")

for (d in c(base_dir, deg_dir, plot_top10_dir,
            plot_sig_dir, plot_clean_dir, go_dir)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

# ---- Assumes pv_only already exists from subclustering block
# If not, re-run:
# pv_only <- subset(pv_cortical_filtered_pnn, idents = "PV+ interneurons")
# pv_only <- NormalizeData(pv_only) %>% FindVariableFeatures() %>%
#            ScaleData() %>% RunPCA(npcs=20) %>%
#            FindNeighbors(dims=1:10) %>% FindClusters(resolution=0.3) %>%
#            RunUMAP(dims=1:10)

pv_only$Genotype_short <- recode(pv_only$Genotype,
                                 "PV-Cre/tdTom"              = "Ctrl",
                                 "PV-Cre/tdTom/Dnmt1 loxP2" = "KO"
)

wt_pv <- subset(pv_only, subset = Genotype_short == "Ctrl")
ko_pv <- subset(pv_only, subset = Genotype_short == "KO")

total_wt_pv <- ncol(wt_pv)
total_ko_pv <- ncol(ko_pv)

# ---- 1a. Proportions within PV population only 
wt_counts_pv <- as.data.frame(table(Idents(wt_pv)))
ko_counts_pv <- as.data.frame(table(Idents(ko_pv)))
colnames(wt_counts_pv) <- c("Cluster", "CellNumber_WT")
colnames(ko_counts_pv) <- c("Cluster", "CellNumber_KO")

merged_pv <- merge(wt_counts_pv, ko_counts_pv, by = "Cluster", all = TRUE)
merged_pv[is.na(merged_pv)] <- 0

merged_pv$Proportion_WT_withinPV <- merged_pv$CellNumber_WT / total_wt_pv
merged_pv$Proportion_KO_withinPV <- merged_pv$CellNumber_KO / total_ko_pv
merged_pv$RelativeDiff_withinPV  <- (
  (merged_pv$Proportion_KO_withinPV - merged_pv$Proportion_WT_withinPV) /
    merged_pv$Proportion_WT_withinPV * 100
)
merged_pv$RelativeDiff_withinPV[!is.finite(merged_pv$RelativeDiff_withinPV)] <- NA

merged_pv$p_value_withinPV <- NA
for (i in seq_len(nrow(merged_pv))) {
  k_wt <- merged_pv$CellNumber_WT[i]
  k_ko <- merged_pv$CellNumber_KO[i]
  if ((k_wt + k_ko) > 0) {
    merged_pv$p_value_withinPV[i] <- prop.test(
      x = c(k_wt, k_ko),
      n = c(total_wt_pv, total_ko_pv),
      alternative = "two.sided"
    )$p.value
  }
}
merged_pv$adj_p_value_withinPV <- p.adjust(merged_pv$p_value_withinPV, method = "fdr")

# ---- 1b. Proportions normalized to total cells in full dataset
pv_cortical_filtered_pnn$Genotype_short <- recode(pv_cortical_filtered_pnn$Genotype,
                                                  "PV-Cre/tdTom"              = "Ctrl",
                                                  "PV-Cre/tdTom/Dnmt1 loxP2" = "KO"
)

total_wt_all <- ncol(subset(pv_cortical_filtered_pnn,
                            subset = Genotype_short == "Ctrl"))
total_ko_all <- ncol(subset(pv_cortical_filtered_pnn,
                            subset = Genotype_short == "KO"))

merged_pv$Proportion_WT_ofTotal <- merged_pv$CellNumber_WT / total_wt_all
merged_pv$Proportion_KO_ofTotal <- merged_pv$CellNumber_KO / total_ko_all
merged_pv$RelativeDiff_ofTotal  <- (
  (merged_pv$Proportion_KO_ofTotal - merged_pv$Proportion_WT_ofTotal) /
    merged_pv$Proportion_WT_ofTotal * 100
)
merged_pv$RelativeDiff_ofTotal[!is.finite(merged_pv$RelativeDiff_ofTotal)] <- NA

merged_pv$p_value_ofTotal <- NA
for (i in seq_len(nrow(merged_pv))) {
  k_wt <- merged_pv$CellNumber_WT[i]
  k_ko <- merged_pv$CellNumber_KO[i]
  if ((k_wt + k_ko) > 0) {
    merged_pv$p_value_ofTotal[i] <- prop.test(
      x = c(k_wt, k_ko),
      n = c(total_wt_all, total_ko_all),
      alternative = "two.sided"
    )$p.value
  }
}
merged_pv$adj_p_value_ofTotal <- p.adjust(merged_pv$p_value_ofTotal, method = "fdr")

# Round percentages
merged_pv$Proportion_WT_withinPV <- round(merged_pv$Proportion_WT_withinPV * 100, 2)
merged_pv$Proportion_KO_withinPV <- round(merged_pv$Proportion_KO_withinPV * 100, 2)
merged_pv$Proportion_WT_ofTotal  <- round(merged_pv$Proportion_WT_ofTotal  * 100, 4)
merged_pv$Proportion_KO_ofTotal  <- round(merged_pv$Proportion_KO_ofTotal  * 100, 4)

write.csv(merged_pv,
          file.path(base_dir, "PV_subcluster_proportion_analysis.csv"),
          row.names = FALSE)

cat("Part 1 done. Proportion table saved to ./pv_deg/\n")
print(merged_pv)


# PV SUBCLUSTER ANALYSIS — PART 2: DEG + VOLCANO PLOTS --------

wt_label  <- "PV-Cre/tdTom"
ko_label  <- "PV-Cre/tdTom/Dnmt1 loxP2"
logfc_thr <- 0.26
padj_thr  <- 0.05

sanitize_id <- function(x) {
  x <- gsub("[^A-Za-z0-9_]", "_", x)
  x <- gsub("_{2,}", "_", x)
  gsub("_$", "", x)
}

build_volcano <- function(de, clust, label_mode = c("top10", "signif", "none")) {
  label_mode <- match.arg(label_mode)
  de$label   <- ""
  
  if (label_mode == "top10") {
    is_up   <- (de$avg_log2FC >= logfc_thr) & (de$p_val_adj < padj_thr)
    is_down <- (de$avg_log2FC <= -logfc_thr) & (de$p_val_adj < padj_thr)
    if (any(is_up & is_down, na.rm = TRUE)) is_down[is_up] <- FALSE
    up_df   <- head(de[is_up,   , drop = FALSE], 10)
    down_df <- head(de[is_down, , drop = FALSE], 10)
    if (nrow(up_df)   > 0) de$label[match(up_df$gene,   de$gene)] <- up_df$gene
    if (nrow(down_df) > 0) de$label[match(down_df$gene, de$gene)] <- down_df$gene
  } else if (label_mode == "signif") {
    de$label[de$diffexpressed != "NO"] <- de$gene[de$diffexpressed != "NO"]
  }
  
  p <- ggplot(de, aes(x = avg_log2FC, y = neg_log10_padj, color = diffexpressed)) +
    geom_point(alpha = 0.6, size = 1.5, show.legend = TRUE) +
    scale_color_manual(
      values = c("UP" = "#b2182b", "DOWN" = "#2166ac", "NO" = "grey70"),
      labels = c("DOWN" = "Downregulated", "NO" = "Not significant", "UP" = "Upregulated")
    ) +
    geom_vline(xintercept = c(-logfc_thr, logfc_thr),
               linetype = "dashed", color = "grey60", linewidth = 0.2) +
    geom_hline(yintercept = -log10(padj_thr),
               linetype = "dashed", color = "grey60", linewidth = 0.2) +
    labs(
      title = paste0("PV sub-cluster ", clust, ": KO vs WT"),
      x     = "Log2 Fold Change",
      y     = "-Log10(Adjusted P-value)",
      color = "Differential Expression"
    ) +
    theme_classic() +
    theme(
      plot.title      = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.title      = element_text(size = 12),
      axis.text       = element_text(size = 10),
      legend.position = "right",
      axis.line       = element_line(color = "black", linewidth = 0.5)
    ) +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))
  
  if (label_mode == "top10" && any(nzchar(de$label))) {
    set.seed(123)
    p <- p + ggrepel::geom_text_repel(
      aes(label = label),
      max.overlaps = 20, size = 3,
      box.padding = 0.5, point.padding = 0.3,
      show.legend = FALSE,
      segment.color = "grey30", segment.size = 0.3,
      min.segment.length = 0
    )
  } else if (label_mode == "signif" && any(nzchar(de$label))) {
    set.seed(123)
    p <- p + ggrepel::geom_text_repel(
      aes(label = label),
      max.overlaps = Inf, size = 2.5,
      box.padding = 0.1, point.padding = 0.1,
      show.legend = FALSE,
      segment.color = NA, segment.size = 0
    )
  }
  p
}

# ---- Iterate over all PV sub-clusters 
pv_clusters <- levels(Idents(pv_only))

for (clust in pv_clusters) {
  
  cells  <- WhichCells(pv_only, idents = clust)
  subobj <- subset(pv_only, cells = cells)
  Idents(subobj) <- "Genotype"
  
  geno_counts <- table(Idents(subobj))
  if (length(geno_counts) < 2 || any(geno_counts < 3)) {
    cat("Skipping", clust, "— insufficient cells in one genotype.\n")
    next
  }
  
  de <- tryCatch(
    FindMarkers(subobj,
                ident.1         = ko_label,
                ident.2         = wt_label,
                min.pct         = 0.1,
                logfc.threshold = logfc_thr,
                test.use        = "wilcox"),
    error = function(e) { message("FindMarkers error: ", e$message); NULL }
  )
  if (is.null(de) || nrow(de) == 0) {
    cat("No DEGs for cluster:", clust, "\n")
    next
  }
  
  de$gene <- rownames(de)
  if (!"avg_log2FC" %in% colnames(de) && "avg_logFC" %in% colnames(de)) {
    de$avg_log2FC <- de$avg_logFC
  }
  de <- de[!grepl("^ENSMUSG", de$gene), , drop = FALSE]
  de <- de[order(de$p_val_adj), , drop = FALSE]
  
  de$diffexpressed <- "NO"
  de$diffexpressed[de$avg_log2FC >= logfc_thr & de$p_val_adj < padj_thr] <- "UP"
  de$diffexpressed[de$avg_log2FC <= -logfc_thr & de$p_val_adj < padj_thr] <- "DOWN"
  de$diffexpressed <- factor(de$diffexpressed, levels = c("UP", "DOWN", "NO"))
  
  de$neg_log10_padj <- -log10(de$p_val_adj)
  finite_vals <- de$neg_log10_padj[is.finite(de$neg_log10_padj)]
  if (length(finite_vals) == 0) finite_vals <- 0
  de$neg_log10_padj[!is.finite(de$neg_log10_padj)] <- max(finite_vals) + 10
  
  safe_clust <- sanitize_id(clust)
  
  write.csv(de,
            file.path(deg_dir,
                      paste0("PV_subcluster_", safe_clust, "_KO_vs_WT.csv")),
            row.names = FALSE)
  message("DEG saved: PV sub-cluster ", clust)
  
  ggsave(file.path(plot_top10_dir,
                   paste0("PV_subcluster_", safe_clust, "_volcano_top10.pdf")),
         build_volcano(de, clust, "top10"),
         width = 10, height = 8, device = "pdf")
  
  ggsave(file.path(plot_sig_dir,
                   paste0("PV_subcluster_", safe_clust, "_volcano_sig_labeled.pdf")),
         build_volcano(de, clust, "signif"),
         width = 10, height = 8, device = "pdf")
  
  ggsave(file.path(plot_clean_dir,
                   paste0("PV_subcluster_", safe_clust, "_volcano_clean.pdf")),
         build_volcano(de, clust, "none"),
         width = 10, height = 8, device = "pdf")
}

cat("Part 2 done. DEG results saved to ./pv_deg/DEG_results/\n")


# PV SUBCLUSTER ANALYSIS — PART 3: GENE ONTOLOGY -------


p_adj_cut       <- 0.05
padj_cut_plot   <- 0.05
top_n           <- 20
wrap_char_limit <- 45
plot_width      <- 9
plot_height     <- 6
plot_dpi        <- 600
panel_rel_width <- 0.6
palette_choice  <- "brewer"
brewer_pal      <- "Reds"
brewer_direction <- 1

enrich_funcs <- list(
  GO_BP = function(ids) enrichGO(ids, OrgDb=org.Mm.eg.db, keyType="ENTREZID",
                                 ont="BP", pAdjustMethod="BH",
                                 pvalueCutoff=p_adj_cut, qvalueCutoff=p_adj_cut,
                                 readable=TRUE),
  GO_MF = function(ids) enrichGO(ids, OrgDb=org.Mm.eg.db, keyType="ENTREZID",
                                 ont="MF", pAdjustMethod="BH",
                                 pvalueCutoff=p_adj_cut, qvalueCutoff=p_adj_cut,
                                 readable=TRUE),
  GO_CC = function(ids) enrichGO(ids, OrgDb=org.Mm.eg.db, keyType="ENTREZID",
                                 ont="CC", pAdjustMethod="BH",
                                 pvalueCutoff=p_adj_cut, qvalueCutoff=p_adj_cut,
                                 readable=TRUE),
  KEGG = function(ids) enrichKEGG(ids, organism="mmu", pAdjustMethod="BH",
                                  pvalueCutoff=p_adj_cut, qvalueCutoff=p_adj_cut),
  Reactome = function(ids) enrichPathway(ids, organism="mouse", pAdjustMethod="BH",
                                         pvalueCutoff=p_adj_cut, qvalueCutoff=p_adj_cut,
                                         readable=TRUE),
  DO = function(ids) enrichDO(ids, ont="DO", pAdjustMethod="BH",
                              pvalueCutoff=p_adj_cut, qvalueCutoff=p_adj_cut,
                              readable=TRUE)
)

compute_fold_enrichment <- function(gr_chr, br_chr) {
  gr <- do.call(rbind, strsplit(gr_chr, "/", fixed = TRUE))
  br <- do.call(rbind, strsplit(br_chr, "/", fixed = TRUE))
  k <- as.numeric(gr[,1]); n <- as.numeric(gr[,2])
  M <- as.numeric(br[,1]); N <- as.numeric(br[,2])
  (k / n) / (M / N)
}

wrap_label_by_chars <- function(x, limit = wrap_char_limit) {
  vapply(x, function(s) {
    if (nchar(s) <= limit) return(s)
    left <- gregexpr("\\s", substr(s, 1, limit))[[1]]
    if (length(left) && any(left > 0)) {
      cut <- max(left[left > 0])
      return(paste0(substr(s, 1, cut-1), "\n",
                    trimws(substr(s, cut+1, nchar(s)))))
    }
    right <- regexpr("\\s", substr(s, limit+1, nchar(s)))
    if (right[1] != -1) {
      cut <- limit + right[1]
      return(paste0(substr(s, 1, cut-1), "\n",
                    trimws(substr(s, cut+1, nchar(s)))))
    }
    paste0(substr(s, 1, limit), "\n", substr(s, limit+1, nchar(s)))
  }, character(1))
}

size_breaks_integer <- function(vals, n = 4) {
  uc <- sort(unique(as.integer(round(vals))))
  if (length(uc) <= n) return(uc)
  br <- as.integer(round(pretty(range(uc), n = n)))
  unique(br[br >= min(uc) & br <= max(uc)])
}

# ---- Step 1: Run enrichment and save CSVs
csv_files <- list.files(deg_dir,
                        pattern = "PV_subcluster_.*_KO_vs_WT\\.csv$",
                        full.names = TRUE)

for (csv_path in csv_files) {
  
  cluster_id <- sub("PV_subcluster_(.*)_KO_vs_WT\\.csv$", "\\1",
                    basename(csv_path))
  de_df      <- read.csv(csv_path, stringsAsFactors = FALSE)
  genes_sym  <- dplyr::filter(de_df, p_val_adj <= p_adj_cut)$gene
  
  if (length(genes_sym) < 5) {
    message("Cluster ", cluster_id, ": fewer than 5 sig. genes, skipping GO.")
    next
  }
  
  tryCatch({
    gene_map   <- bitr(genes_sym, fromType="SYMBOL",
                       toType="ENTREZID", OrgDb=org.Mm.eg.db)
    entrez_ids <- unique(gene_map$ENTREZID)
    
    if (length(entrez_ids) < 3) {
      message("Cluster ", cluster_id, ": fewer than 3 Entrez IDs, skipping GO.")
      next
    }
    
    for (fname in names(enrich_funcs)) {
      tryCatch({
        ego    <- enrich_funcs[[fname]](entrez_ids)
        if (is.null(ego) || nrow(ego@result) == 0) {
          message("Cluster ", cluster_id, " - ", fname, ": no terms found.")
          next
        }
        res_df <- as.data.frame(ego@result)
        eps    <- .Machine$double.xmin
        res_df$NegLog10Padj <- -log10(pmax(as.numeric(res_df$p.adjust), eps))
        
        keep_cols <- intersect(
          c("ID","Description","GeneRatio","BgRatio","pvalue",
            "p.adjust","qvalue","geneID","Count","NegLog10Padj"),
          colnames(res_df)
        )
        out_df <- res_df[order(res_df$p.adjust), keep_cols]
        out_csv <- file.path(go_dir,
                             sprintf("PV_subcluster_%s_%s.csv",
                                     cluster_id, fname))
        write.csv(out_df, out_csv, row.names = FALSE)
        message("Saved: ", basename(out_csv))
      }, error = function(e) {
        message("Error in ", fname, " for cluster ", cluster_id, ": ", e$message)
      })
    }
  }, error = function(e) {
    message("Gene mapping error for cluster ", cluster_id, ": ", e$message)
  })
}

# ---- Step 2: Plot GO dot plots from CSVs 
go_csvs <- list.files(go_dir,
                      pattern = "PV_subcluster_.*_GO_(BP|MF|CC)\\.csv$",
                      full.names = TRUE)

for (csv_path in go_csvs) {
  df <- read.csv(csv_path, stringsAsFactors = FALSE)
  
  df$FoldEnrichment <- compute_fold_enrichment(df$GeneRatio, df$BgRatio)
  
  if (!"GeneCount" %in% names(df)) {
    df$GeneCount <- if ("Count" %in% names(df)) df$Count else
      vapply(strsplit(df$geneID, "/", fixed=TRUE), length, integer(1))
  }
  df$GeneCount <- as.integer(round(df$GeneCount))
  
  if (!"NegLog10Padj" %in% names(df)) {
    df$NegLog10Padj <- -log10(pmax(as.numeric(df$p.adjust), .Machine$double.xmin))
  }
  
  plot_df <- df[is.finite(df$FoldEnrichment) &
                  !is.na(df$p.adjust) &
                  df$p.adjust < padj_cut_plot, ]
  plot_df <- plot_df[order(plot_df$p.adjust), ]
  
  if (nrow(plot_df) == 0) {
    message("No terms to plot in ", basename(csv_path))
    next
  }
  if (nrow(plot_df) > top_n) plot_df <- head(plot_df, top_n)
  
  plot_df$Label          <- wrap_label_by_chars(plot_df$Description)
  plot_df$GeneCountCapped <- plot_df$GeneCount
  size_br                 <- size_breaks_integer(plot_df$GeneCount)
  
  base_plot <- ggplot(plot_df,
                      aes(x = FoldEnrichment,
                          y = reorder(Label, FoldEnrichment))) +
    geom_point(aes(size = GeneCountCapped, color = NegLog10Padj)) +
    scale_size_continuous(name   = "# Genes",
                          breaks = size_br,
                          labels = as.character,
                          range  = c(2.5, 7.0)) +
    scale_x_continuous(name   = "Fold enrichment",
                       expand = expansion(mult = c(0.08, 0.05))) +
    scale_y_discrete(name   = NULL,
                     expand = expansion(add = 0.6)) +
    scale_color_distiller(name      = "-log10(p.adj)",
                          palette   = brewer_pal,
                          direction = brewer_direction) +
    labs(title = sprintf("%s (top %d, padj < %.2f)",
                         sub("\\.csv$", "", basename(csv_path)),
                         nrow(plot_df), padj_cut_plot)) +
    theme_classic(base_size = 12) +
    theme(
      panel.grid  = element_blank(),
      axis.line   = element_line(color = "black"),
      plot.title  = element_text(hjust = 0.5),
      axis.text.y = element_text(size = 10)
    )
  
  legend_grob   <- cowplot::get_legend(base_plot +
                                         theme(legend.position = "right"))
  combined_plot <- cowplot::plot_grid(
    base_plot + theme(legend.position = "none"),
    legend_grob,
    ncol = 2, align = "h",
    rel_widths = c(panel_rel_width, 1 - panel_rel_width)
  )
  
  out_pdf <- file.path(go_dir,
                       paste0(sub("\\.csv$", "", basename(csv_path)),
                              sprintf("_top%d_dotplot.pdf", top_n)))
  ggsave(out_pdf, combined_plot,
         device = pdf, width = plot_width,
         height = plot_height, units = "in")
  message("Saved: ", basename(out_png))
}

cat("Part 3 done. GO results saved to ./pv_deg/GO_results/\n")

###############
#### OLIGO UMAP: per-cluster WT/KO split panels ----

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

# ---- Overall cluster-colored UMAP (unchanged)
p_oligo_clusters <- DimPlot(oligo_only, reduction = "umap",
                            label = FALSE, repel = TRUE, pt.size  = 0.3) +
  scale_colour_manual(values = oligo_colors) +
  ggtitle("Oligodendrocyte lineage sub-types") +
  theme(legend.position = "right")

# ---- Overall genotype UMAP
p_oligo_genotype <- DimPlot(oligo_only, reduction = "umap",
                            group.by = "Genotype_short",
                            cols     = c("Ctrl" = "#c6c6c6", "KO" = "#7FCDCD"),
                            pt.size  = 0.3) +
  ggtitle("Oligodendrocyte UMAP: Ctrl vs KO") +
  theme(legend.position = "right")

#### Per-cluster WT/KO panels ----
# Strategy: for each cluster, build a metadata column that is
# "Ctrl" / "KO" for cells IN that cluster, and "Other" for the rest.
# Other cells are plotted in light grey as background context.

cluster_plot_list <- list()

for (cl in oligo_clusters) {
  
  # Build a per-plot color label column
  oligo_only$plot_group <- ifelse(
    as.character(Idents(oligo_only)) == cl,
    as.character(oligo_only$Genotype_short),
    "Other"
  )
  
  # Order so "Other" is drawn first (background), then Ctrl, then KO on top
  oligo_only$plot_group <- factor(oligo_only$plot_group,
                                  levels = c("Other", "Ctrl", "KO"))
  
  p <- DimPlot(oligo_only,
               reduction = "umap",
               group.by  = "plot_group",
               order     = c("KO", "Ctrl", "Other"),   # KO drawn on top
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
    guides(color = guide_legend(
      override.aes = list(size = 3, alpha = 1)
    ))
  
  cluster_plot_list[[cl]] <- p
}

# ---- Save each panel individually 
short_names <- c(
  "Oligodendrocyte precursor cells" = "OPC",
  "Perineuronal oligodendrocytes"   = "PNOL",
  "Newly mature oligodendrocytes"   = "NFOL",
  "Myelinating oligodendrocytes"    = "MOL"
)

for (cl in oligo_clusters) {
  pdf(paste0("OligoOnly_UMAP_", short_names[cl], "_WT_KO.pdf"),
      width = 6, height = 5)
  print(cluster_plot_list[[cl]])
  dev.off()
}

# ---- Save all four panels in one combined figure 
library(patchwork)

p_combined_percluster <- (
  cluster_plot_list[[1]] + cluster_plot_list[[2]] +
    cluster_plot_list[[3]] + cluster_plot_list[[4]]
) +
  plot_layout(ncol = 2) +
  plot_annotation(
    title    = "Oligodendrocyte lineage: Ctrl vs KO per sub-type",
    theme    = theme(plot.title = element_text(hjust = 0.5,
                                               face = "bold", size = 14))
  )

pdf("OligoOnly_UMAP_percluster_WT_KO.pdf", width = 11, height = 8)
print(p_combined_percluster)
dev.off()

# ---- Save original combined overview
pdf("OligoOnly_UMAP_combined.pdf", width = 14, height = 5)
print(p_oligo_clusters + p_oligo_genotype)
dev.off()

cat("All oligo UMAP plots saved to ./revision/\n")


################
#### FULL DATASET WT/KO UMAP — split into two panels ----
pv_cortical_filtered_pnn$Genotype_short <- recode(pv_cortical_filtered_pnn$Genotype,
                                                  "PV-Cre/tdTom"              = "Ctrl",
                                                  "PV-Cre/tdTom/Dnmt1 loxP2" = "KO"
)

p_wt_ko_split <- DimPlot(pv_cortical_filtered_pnn,
                         reduction = "umap",
                         group.by  = "Genotype_short",
                         split.by  = "Genotype_short",
                         cols      = c("Ctrl" = "#c6c6c6", "KO" = "#7FCDCD"),
                         pt.size   = 0.1,
                         alpha     = 1) +
  ggtitle(NULL) +
  theme(
    legend.position = "none",
    strip.text      = element_text(size = 13, face = "bold")  # panel titles "Ctrl" / "KO"
  )

pdf("FullDataset_UMAP_WT_KO_split.pdf",
    width = 12, height = 6)
p_wt_ko_split
dev.off()


###############
#### Extract top 50 upregulated and downregulated genes per cluster from cortical (no pnn separation) -----
all_markers_cortical <- FindAllMarkers(pv_cortical_filtered,
                                       min.pct = 0.25,          # Gene must be detected in ≥25% cells
                                       logfc.threshold = 0.25,  # Minimum fold change threshold
                                       only.pos = FALSE,         # Keep both upregulated and downregulated markers
                                       test.use = "wilcox"      # Default Wilcoxon Rank Sum test
)

top50_up <- all_markers_cortical %>%
  filter(p_val_adj < 0.05, avg_log2FC > 0) %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 50) %>%
  mutate(direction = "upregulated") %>%
  ungroup()

# Extract top 50 downregulated genes per cluster
top50_down <- all_markers_cortical %>%
  filter(p_val_adj < 0.05, avg_log2FC < 0) %>%
  group_by(cluster) %>%
  slice_min(order_by = avg_log2FC, n = 50) %>%
  mutate(direction = "downregulated") %>%
  ungroup()

# Combine both into a single dataframe
top50_combined <- bind_rows(top50_up, top50_down) %>%
  arrange(cluster, direction, desc(avg_log2FC))

# Write to CSV
write.csv(top50_combined, "top50_markers_cortical.csv", row.names = FALSE)

message("Done! File saved as top50_markers_cortical_pnn.csv")

###############
#### DEG analysis: Newly Mature Oligos vs Myelinating Oligos ------

library(Seurat)
library(dplyr)
library(readr)


seurat_obj          <- pv_cortical_filtered
cluster_col         <- "sub.cluster"
cluster_new         <- "Newly mature oligodendrocytes"
cluster_myelinating <- "Myelinating oligodendrocytes"
output_file         <- "DEGs_NewlyMatureOligo_vs_MyelinatingOligo.csv"

# Set identity to sub.cluster
Idents(seurat_obj) <- cluster_col

# Verify the cluster labels exist
stopifnot(
  cluster_new %in% levels(Idents(seurat_obj)),
  cluster_myelinating %in% levels(Idents(seurat_obj))
)

# Run FindMarkers: Newly Mature Oligo (ident.1) vs Myelinating Oligo (ident.2)
# Positive FC  = higher in Newly Mature Oligo
# Negative FC  = higher in Myelinating Oligo
degs <- FindMarkers(
  object          = seurat_obj,
  ident.1         = cluster_new,
  ident.2         = cluster_myelinating,
  test.use        = "wilcox",
  logfc.threshold = 0,       # return all genes, filter afterwards
  min.pct         = 0.1,     # detected in ≥10% of cells in either cluster
  only.pos        = FALSE
)

# Add gene names as column, sort by avg_log2FC
degs <- degs %>%
  tibble::rownames_to_column("gene") %>%
  arrange(desc(avg_log2FC))

# Add direction column
degs <- degs %>%
  mutate(direction = case_when(
    p_val_adj < 0.05 & avg_log2FC > 0  ~ "upregulated_in_NewlyMatureOligo",
    p_val_adj < 0.05 & avg_log2FC < 0  ~ "upregulated_in_MyelinatingOligo",
    TRUE                                ~ "not_significant"
  ))

# Print summary
cat("Total genes tested:                  ", nrow(degs), "\n")
cat("Significant DEGs (padj < 0.05):      ", sum(degs$p_val_adj < 0.05), "\n")
cat("Higher in Newly Mature Oligo:        ", sum(degs$direction == "upregulated_in_NewlyMatureOligo"), "\n")
cat("Higher in Myelinating Oligo:         ", sum(degs$direction == "upregulated_in_MyelinatingOligo"), "\n")

# Write to CSV
write_csv(degs, output_file)
cat("Saved to:", output_file, "\n")

##############
#### Cell chat for wildtype --------------------------------------------------------

data.input <- wt_cortical[["RNA"]]$data # normalized data matrix
labels <- Idents(wt_cortical)
meta <- data.frame(labels = labels, row.names = names(labels))
cellchat_wt_cortical <- createCellChat(object = wt_cortical, group.by = "ident", assay = "RNA")


CellChatDB <- CellChatDB.mouse # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)

CellChatDB.use <- CellChatDB
cellchat_wt_cortical@DB <- CellChatDB.use

cellchat_wt_cortical <- subsetData(cellchat_wt_cortical)
cellchat_wt_cortical <- updateCellChat(cellchat_wt_cortical)
#future::plan("multisession", workers = 12) # do parallel
#future::plan("sequential")

cellchat_wt_cortical <- identifyOverExpressedGenes(cellchat_wt_cortical)
cellchat_wt_cortical <- identifyOverExpressedInteractions(cellchat_wt_cortical)
length(cellchat_wt_cortical@LR$LRsig)

ptm = Sys.time()
execution.time = Sys.time() - ptm
print(as.numeric(execution.time, units = "secs"))

cellchat_wt_cortical <- computeCommunProb(cellchat_wt_cortical, type = "triMean")
cellchat_wt_cortical <- filterCommunication(cellchat_wt_cortical, min.cells = 10)

cellchat_wt_cortical <- computeCommunProbPathway(cellchat_wt_cortical)

cellchat_wt_cortical <- aggregateNet(cellchat_wt_cortical)
execution.time = Sys.time() - ptm
print(as.numeric(execution.time, units = "secs"))


#### Cell chat for knockout --------------------------------------------------------

data.input <- ko_cortical[["RNA"]]$data # normalized data matrix
labels <- Idents(ko_cortical)
meta <- data.frame(labels = labels, row.names = names(labels))
cellchat_ko_cortical <- createCellChat(object = ko_cortical, group.by = "ident", assay = "RNA")


CellChatDB <- CellChatDB.mouse # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)

CellChatDB.use <- CellChatDB
cellchat_ko_cortical@DB <- CellChatDB.use

cellchat_ko_cortical <- subsetData(cellchat_ko_cortical)
cellchat_ko_cortical <- updateCellChat(cellchat_ko_cortical)
#future::plan("multisession", workers = 12) # do parallel
#future::plan("sequential")

cellchat_ko_cortical <- identifyOverExpressedGenes(cellchat_ko_cortical)
cellchat_ko_cortical <- identifyOverExpressedInteractions(cellchat_ko_cortical)
length(cellchat_ko_cortical@LR$LRsig)

ptm = Sys.time()
execution.time = Sys.time() - ptm
print(as.numeric(execution.time, units = "secs"))

cellchat_ko_cortical <- computeCommunProb(cellchat_ko_cortical, type = "triMean")
cellchat_ko_cortical <- filterCommunication(cellchat_ko_cortical, min.cells = 10)

cellchat_ko_cortical <- computeCommunProbPathway(cellchat_ko_cortical)

cellchat_ko_cortical <- aggregateNet(cellchat_ko_cortical)
execution.time = Sys.time() - ptm
print(as.numeric(execution.time, units = "secs"))


#### Change labels of the clusters to abbreviations -------
old_cellchat_names <- c(
  "Layer 2/3 IT neurons",
  "Layer 4 sensory neurons",
  "Layer 5a IT neurons",
  "Layer 5b PT neurons",
  "Layer 5/6 IT neurons",
  "Layer 6a corticothalamic neurons",
  "Layer 6b neurons",
  "Deep-layer extratelencephalic neurons",
  "Corticospinal neurons (Type I)",
  "Corticospinal neurons (Type II)",
  "PV+ interneurons",
  "SST+ interneurons",
  "VIP+ interneurons",
  "Astrocytes",
  "Oligodendrocyte precursor cells",
  "Newly formed oligodendrocytes",
  "Myelinating oligodendrocytes",
  "Microglia",
  "Endothelial cells",
  "Leptomeningeal cells",
  "Meningeal fibroblasts"
)

new_names <- c(
  "L2/3 IT",
  "L4 Sensory",
  "L5a IT",
  "L5b PT",
  "L5/6 IT",
  "L6a CT",
  "L6b",
  "Deep ET",
  "CSN Type I",
  "CSN Type II",
  "PV+ Int",
  "SST+ Int",
  "VIP+ Int",
  "Astrocyte",
  "OPC",
  "NFOL",
  "MOL",
  "Microglia",
  "Endothelial",
  "Leptomeningeal FB",
  "Meningeal FB"
)

cellchat_wt_cortical <- updateClusterLabels(
  cellchat_wt_cortical,
  old.cluster.name = old_cellchat_names,
  new.cluster.name = new_names
)

cellchat_ko_cortical <- updateClusterLabels(
  cellchat_ko_cortical,
  old.cluster.name = old_cellchat_names,
  new.cluster.name = new_names
)

# Verify
levels(cellchat_wt_cortical@idents)
levels(cellchat_ko_cortical@idents)

# Save
saveRDS(cellchat_wt_cortical, "cellchat_wt_cortical_only.rds")
saveRDS(cellchat_ko_cortical, "cellchat_ko_cortical_only.rds")

cellchat_wt_cortical <- readRDS("cellchat_wt_cortical_only.rds")
cellchat_ko_cortical <- readRDS("cellchat_ko_cortical_only.rds")

#### Remove unwanted clusters from CellChat -------
cellchat_wt_cortical <- subsetCellChat(cellchat_wt_cortical, idents.use = c("CSN Type I", "CSN Type II", "Endothelial", "Leptomeningeal FB", "Meningeal FB"), invert = TRUE)
cellchat_ko_cortical <- subsetCellChat(cellchat_ko_cortical, idents.use = c("CSN Type I", "CSN Type II", "Endothelial", "Leptomeningeal FB", "Meningeal FB"), invert = TRUE)

#### Merge objects and re-compute centrality --------
# For WT
cellchat_wt_cortical <- netAnalysis_computeCentrality(cellchat_wt_cortical, slot.name = "netP")
# For KO
cellchat_ko_cortical <- netAnalysis_computeCentrality(cellchat_ko_cortical, slot.name = "netP")

# List your CellChat objects
object.list_cortical_new <- list(WT = cellchat_wt_cortical, KO = cellchat_ko_cortical)

# Merge for comparison
cellchat_merged_cortical <- mergeCellChat(object.list_cortical_new, add.names = names(object.list_cortical_new))

### CHECK LEVELS
levels(cellchat_merged_cortical@idents$joint)
levels(cellchat_wt_cortical@idents)
levels(cellchat_ko_cortical@idents)

##############
#### Interaction weight and counts ------------------

# List your CellChat objects
# Your original list of CellChat objects
object.list_cortical_pnn <- list(WT = cellchat_wt_cortical, KO = cellchat_ko_cortical)

### Interaction number and weight

par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat_merged_cortical, label.edge= F, weight.scale = T)
netVisual_diffInteraction(cellchat_merged_cortical, label.edge= F, weight.scale = T, measure = "weight")

gg1 <- netVisual_heatmap(cellchat_merged_cortical)
gg2 <- netVisual_heatmap(cellchat_merged_cortical, measure = "weight")
gg1 + gg2
pdf("interaction_weight_cortical.pdf", width = 8, height = 6)
gg2
dev.off()

## Quantify the matrix
# 1. Extract the raw count and weight matrices for each group
nets <- cellchat_merged_cortical@net       # list of two elements
mat_count_ctrl   <- nets[[1]]$count         # WT counts
mat_count_ko     <- nets[[2]]$count         # KO counts
mat_weight_ctrl  <- nets[[1]]$weight        # WT weights
mat_weight_ko    <- nets[[2]]$weight        # KO weights

# 2. Compute differential matrices (KO minus WT)
diff_count  <- mat_count_ko  - mat_count_ctrl
diff_weight <- mat_weight_ko - mat_weight_ctrl

# 3. (Optional) Assign row/column names if lost
rownames(diff_count)  <- rownames(mat_count_ctrl)
colnames(diff_count)  <- colnames(mat_count_ctrl)
rownames(diff_weight) <- rownames(mat_weight_ctrl)
colnames(diff_weight) <- colnames(mat_weight_ctrl)

# 4. Save to CSV for inspection
write.csv(diff_count,
          "netVisual_diff_counts_KO_vs_WT.csv",
          row.names = TRUE)
write.csv(diff_weight,
          "netVisual_diff_weights_KO_vs_WT.csv",
          row.names = TRUE)


#### Circle interaction map (PV and excitatory)---------
suppressPackageStartupMessages({
  library(CellChat)
  library(igraph)
})

par(mfrow = c(1, 1), xpd = NA)

edge_width_max <- 2
weight_cap     <- 0.2
arrow_width    <- edge_width_max / 2
arrow_size     <- 0.3 * edge_width_max
label_offset_in <- 0.43

# Requested order (cluster 9 first); pass either indices or names
sources.use <- c(9,1,2,3,4,5,6,7,8,10,11)
targets.use <- c(9,1,2,3,4,5,6,7,8,10,11)

# Inputs
object  <- cellchat_merged_cortical
comp    <- c(1, 2)
measure <- "weight"  # "count" or "weight"

# Access nets and ensure dimnames via idents
net1 <- object@net[[comp[1]]]
net2 <- object@net[[comp[2]]]
mat1 <- net1[[measure]]
mat2 <- net2[[measure]]

fix_dimnames <- function(m, idents_levels) {
  if (is.null(rownames(m)) && length(idents_levels) == nrow(m)) rownames(m) <- idents_levels
  if (is.null(colnames(m)) && length(idents_levels) == ncol(m)) colnames(m) <- idents_levels
  m
}
lvl1 <- if (is.list(object@idents)) levels(object@idents[[comp[1]]]) else levels(object@idents)
lvl2 <- if (is.list(object@idents)) levels(object@idents[[comp[2]]]) else levels(object@idents)
mat1 <- fix_dimnames(mat1, lvl1)
mat2 <- fix_dimnames(mat2, lvl2)

group_names_all <- rownames(mat1)
if (is.null(group_names_all)) stop("Group names are missing; ensure idents have levels or matrices have dimnames.")

# Map requested indices to names using the original names
map_to_names <- function(v, names_all) {
  if (is.numeric(v)) {
    if (any(v < 1 | v > length(names_all))) stop("Index out of range in sources.use/targets.use")
    names_all[v]
  } else {
    as.character(v)
  }
}
sources_names <- map_to_names(sources.use, group_names_all)
targets_names <- map_to_names(targets.use, group_names_all)

# Relevel idents so the first requested name becomes the first sector; reorder all 2D nets accordingly
first_label <- sources_names[1]
new_levels <- c(first_label, setdiff(group_names_all, first_label))

if (is.list(object@idents)) {
  for (k in seq_along(object@idents)) {
    cur_lev <- levels(object@idents[[k]])
    keep_lev <- new_levels[new_levels %in% cur_lev]
    object@idents[[k]] <- factor(object@idents[[k]], levels = keep_lev)
  }
} else {
  cur_lev <- levels(object@idents)
  keep_lev <- new_levels[new_levels %in% cur_lev]
  object@idents <- factor(object@idents, levels = keep_lev)
}

mat_names_2d <- c("count","sum","weight","count.merged","weight.merged")
for (i in seq_along(object@net)) {
  for (mn in intersect(mat_names_2d, names(object@net[[i]]))) {
    m <- object@net[[i]][[mn]]
    rn <- rownames(m); cn <- colnames(m)
    if (!is.null(rn) && !is.null(cn)) {
      keep <- new_levels[new_levels %in% rn & new_levels %in% cn]
      object@net[[i]][[mn]] <- m[keep, keep, drop = FALSE]
    }
  }
}

# Refresh matrices and valid group names after reordering
net1 <- object@net[[comp[1]]]
net2 <- object@net[[comp[2]]]
mat1 <- net1[[measure]]
mat2 <- net2[[measure]]
group_names_all <- rownames(mat1)

# Keep only requested names that still exist
sources_names <- sources_names[sources_names %in% group_names_all]
targets_names <- targets_names[targets_names %in% group_names_all]

# Plot with default CellChat colors and asp=1 for perfect circle
netVisual_diffInteraction(
  object             = object,
  comparison         = comp,
  measure            = measure,
  color.edge         = c("#b2182b", "#2166ac"),
  weight.scale       = FALSE,
  edge.weight.max    = weight_cap,
  edge.width.max     = edge_width_max,
  sources.use        = sources_names,
  targets.use        = targets_names,
  remove.isolate     = TRUE,
  vertex.label.cex   = 1e-6,
  vertex.label.color = NA,
  title.name         = "",
  margin             = 0.8,
  alpha.edge         = 1,
  shape              = "circle",
  arrow.width        = arrow_width,
  arrow.size         = arrow_size
)

# Function to wrap multi-word labels into two lines
wrap_label <- function(label, max_width = 15) {
  # Count number of words
  words <- strsplit(label, " ")[[1]]
  if (length(words) <= 1) {
    return(label)
  }
  # Use strwrap to split into lines
  wrapped <- strwrap(label, width = max_width)
  # Combine with newline
  paste(wrapped, collapse = "\n")
}

# Label placement aligned to plotted order among active nodes
node_names_sel <- unique(c(sources_names, setdiff(targets_names, sources_names)))
sub1 <- mat1[node_names_sel, node_names_sel, drop = FALSE]
sub2 <- mat2[node_names_sel, node_names_sel, drop = FALSE]
sub1z <- ifelse(is.na(sub1), 0, sub1)
sub2z <- ifelse(is.na(sub2), 0, sub2)

deg_any <- (rowSums(sub1z != 0) + colSums(sub1z != 0) +
              rowSums(sub2z != 0) + colSums(sub2z != 0)) > 0
active_names <- node_names_sel[deg_any]

order_active_names <- unique(c(intersect(sources_names, active_names),
                               setdiff(intersect(targets_names, active_names),
                                       intersect(sources_names, active_names))))
n_active <- length(order_active_names)

if (n_active > 0) {
  g_dummy <- igraph::make_empty_graph(n = n_active, directed = TRUE)
  coords  <- igraph::layout_in_circle(g_dummy, order = seq_len(n_active))
  x <- coords[,1]; y <- coords[,2]
  theta <- atan2(y, x); ux <- cos(theta); uy <- sin(theta)
  
  dx_in <- grconvertX(1, "inches", "user") - grconvertX(0, "inches", "user")
  dy_in <- grconvertY(1, "inches", "user") - grconvertY(0, "inches", "user")
  vx <- ux * dx_in; vy <- uy * dy_in
  vnorm <- sqrt(vx^2 + vy^2)
  x_lab <- x + label_offset_in * vx / vnorm
  y_lab <- y + label_offset_in * vy / vnorm
  
  deg <- theta * 180 / pi
  pos_right  <- (deg >= -45 & deg < 45)
  pos_top    <- (deg >= 45  & deg < 135)
  pos_left   <- (deg >= 135 | deg <= -135)
  pos_bottom <- (deg > -135 & deg < -45)
  
  adjx <- rep(0.5, n_active); adjy <- rep(0.5, n_active)
  adjx[pos_right]  <- 0;   adjy[pos_right]  <- 0.5
  adjx[pos_left]   <- 0.9; adjy[pos_left]   <- 0.5
  adjx[pos_top]    <- 0.5; adjy[pos_top]    <- 0
  adjx[pos_bottom] <- 0.5; adjy[pos_bottom] <- 1
  
  # Wrap labels for multi-word labels
  wrapped_labels <- sapply(order_active_names, wrap_label, max_width = 15)
  
  for (i in seq_len(n_active)) {
    text(
      x = x_lab[i], y = y_lab[i], labels = wrapped_labels[i],
      adj = c(adjx[i], adjy[i]), srt = 0,
      cex = 1.3, col = "black", xpd = NA
    )
  }
}


#### Plot UMAP ---------------

cluster_ids <- levels(Idents(pv_cortical_filtered))
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
    "#9B2335",  # 17 faded burgundy
    "#C17BAE",  # 18 dusty orchid
    "#DECF3F",  # 19 olive gold
    "#789262",  # 20 sage green
    "#BC243C"   # 21 dark raspberry
  ),
  cluster_ids
)
DimPlot(pv_cortical_filtered, reduction = "umap",repel = TRUE, pt.size = 0.01) + 
  scale_colour_manual(values = colors) +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 3)))

umap_1 <- DimPlot(pv_cortical_filtered, reduction = "umap",repel = TRUE, pt.size = 0.1) + 
  scale_colour_manual(values = colors) +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 3)))


pdf("umap.pdf", width = 12, height = 10)
umap_1
dev.off()

############### EXTRAS ##############

# DEG Excel list ------
setwd("~/Downloads/Seurat_CellChat/revision/cellchat_new/DE_results_volcano_filtered")

df <- read.csv("SourceData_DEG_Perineuronal_oligodendrocytes.csv")

df <- df |>
  select(-cluster) |>
  mutate(
    diffexpressed = case_when(
      avg_log2FC >= 0.26  & p_val_adj < 0.05 ~ "UP",
      avg_log2FC <= -0.26 & p_val_adj < 0.05 ~ "DOWN",
      TRUE ~ "NO"
    ),
    neg_log10_padj = -log10(p_val_adj)
  ) |>
  select(p_val, avg_log2FC, pct.1, pct.2, p_val_adj, gene, diffexpressed, neg_log10_padj)

write.csv(df, "cluster_Perineuronal_oligodendrocytes_KO_vs_WT_all_genes_filtered.csv", row.names = FALSE)


library(openxlsx)


# Get all CSV files in the current working directory
csv_files <- list.files(pattern = "\\.csv$", full.names = TRUE)

if (length(csv_files) == 0) {
  stop("No CSV files found in the current directory")
}

# Create a workbook
wb <- createWorkbook()

# Loop through each CSV file and add it as a sheet
for (i in seq_along(csv_files)) {
  # Read the CSV
  df <- read.csv(csv_files[i])
  
  # Extract filename without .csv extension
  file_name <- sub("\\.csv$", "", basename(csv_files[i]))
  
  # Generate sheet name as number
  sheet_name <- as.character(i)
  
  # Add worksheet
  addWorksheet(wb, sheet_name)
  
  # Write filename as first row
  writeData(wb, sheet_name, file_name, startRow = 1, startCol = 1)
  
  # Write data starting from row 3 (leaving row 2 blank for spacing)
  writeData(wb, sheet_name, df, startRow = 3)
}

# Save the workbook
saveWorkbook(wb, "combined_data.xlsx", overwrite = TRUE)

cat("Excel file 'combined_data.xlsx' created successfully with", length(csv_files), "sheets\n")

setwd("~/Downloads/Seurat_CellChat/revision/cellchat_new")
