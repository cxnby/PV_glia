# ── 0. Install dependencies ───────────────────────────────────────────────────
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("miloR", "SingleCellExperiment", "scater", "edgeR"))
install.packages(c("dplyr", "ggplot2", "patchwork"))

library(miloR)
library(SingleCellExperiment)
library(scater)
library(edgeR)
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(ggrepel)

# ── 1. Save current Idents and convert Seurat → Milo ─────────────────────────
pv_cortical_filtered$cell_annotation <- Idents(pv_cortical_filtered)

milo_obj <- as.SingleCellExperiment(pv_cortical_filtered)
milo_obj <- Milo(milo_obj)

# ── 2. Build KNN graph ────────────────────────────────────────────────────────
# d = number of PCs used; adjust if you used more/fewer in your Seurat analysis
milo_obj <- buildGraph(milo_obj, k = 30, d = 20, reduced.dim = "PCA")

# ── 3. Define neighbourhoods ──────────────────────────────────────────────────
# prop = 0.1 is recommended for datasets < 30k cells
milo_obj <- makeNhoods(milo_obj, prop = 0.1, k = 30, d = 20,
                       refined = TRUE, reduced_dims = "PCA")

# Check neighbourhood size distribution — mean should be > 5 × n_samples (5×4=20)
# If mean is too low, increase k above and rerun steps 2-3
plotNhoodSizeHist(milo_obj)

# ── 4. Count cells per sample per neighbourhood ───────────────────────────────
milo_obj <- countCells(milo_obj,
                       meta.data = as.data.frame(colData(milo_obj)),
                       sample = "samples")

# ── 5. Build the design matrix (one row per sample) ───────────────────────────
# Extract unique sample-level metadata
design_df <- as.data.frame(colData(milo_obj)) %>%
  dplyr::select(samples, Genotype, Tissue) %>%
  distinct()

rownames(design_df) <- design_df$samples

# Make sure Genotype and Tissue are factors with the right reference level
# WT is the reference (baseline) for Genotype
design_df$Genotype <- factor(design_df$Genotype,
                             levels = c("PV-Cre/tdTom", "PV-Cre/tdTom/Dnmt1 loxP2"))
# Somato-ctx as reference region — change if preferred
design_df$Tissue <- factor(design_df$Tissue)

design_df

# ── 6. Calculate neighbourhood distances (needed for spatial FDR) ─────────────
# This is the slowest step
milo_obj <- calcNhoodDistance(milo_obj, d = 20, reduced.dim = "PCA")

# ── 7. Test for differential abundance ───────────────────────────────────────
# ~ Tissue accounts for region effect; Genotype is the variable of interest
da_results <- testNhoods(milo_obj,
                         design = ~ Tissue + Genotype,
                         design.df = design_df,
                         reduced.dim = "PCA")

# Quick sanity check — p-value histogram should be roughly uniform with a peak near 0
ggplot(da_results, aes(PValue)) + geom_histogram(bins = 50)

# Volcano plot (each point = one neighbourhood, not one cell)
ggplot(da_results, aes(logFC, -log10(SpatialFDR))) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "red") +
  labs(title = "KO vs WT differential abundance",
       subtitle = "Tissue (region) modeled as covariate")

# ── 8. Annotate neighbourhoods with your cell type labels ─────────────────────
da_results <- annotateNhoods(milo_obj, da_results,
                             coldata_col = "cell_annotation")

# Inspect label purity — neighbourhoods with < 70% one cell type are "Mixed"
ggplot(da_results, aes(cell_annotation_fraction)) + geom_histogram(bins = 50)
da_results$cell_annotation <- ifelse(da_results$cell_annotation_fraction < 0.7,
                                     "Mixed", da_results$cell_annotation)

# ── 9. Visualize DA on UMAP ───────────────────────────────────────────────────
milo_obj <- buildNhoodGraph(milo_obj)

umap_pl <- plotReducedDim(milo_obj, dimred = "UMAP",
                          colour_by = "cell_annotation",
                          text_by  = "cell_annotation",
                          text_size = 2.5, point_size = 0.3) +
  guides(fill = "none") +
  ggtitle("Cell type annotations")

nhood_pl <- plotNhoodGraphDA(milo_obj, da_results,
                             layout = "UMAP", alpha = 0.05) +
  ggtitle("DA neighbourhoods (KO vs WT)")

umap_pl + nhood_pl + plot_layout(guides = "collect")

# ── 10. Beeswarm plot: DA fold-changes per cell type ─────────────────────────
# Positive logFC = more abundant in KO; negative = more abundant in WT
plotDAbeeswarm(da_results, group.by = "cell_annotation") +
  labs(title = "DA fold-changes per cell type (KO vs WT)")

# beeswarm_raw <- plotDAbeeswarm(da_results, group.by = "cell_annotation") +
#   labs(title = "DA fold-changes per cell type (KO vs WT)")
# ggsave(
#   filename = "DA_beeswarm_raw.pdf",
#   plot     = beeswarm_raw,
#   width    = 10,
#   height   = 9
# )

beeswarm_raw <- plotDAbeeswarm(da_results %>% filter(cell_annotation != "Mixed"), group.by = "cell_annotation")

# Set dot size — 0.8 is a good starting point, adjust to taste
beeswarm_raw$layers[[1]]$aes_params$size <- 0.8

ggsave("DA_beeswarm_KO_vs_WT.pdf", plot = beeswarm_raw,
       width = 13, height = 9)


### Extract the data
# View what's in it
head(da_results)
colnames(da_results)

# Save as CSV
write.csv(da_results, file = "da_results_milo.csv", row.names = FALSE)

#### Combine with proportion analysis ------
library(dplyr)
library(ggplot2)
library(ggrepel)
library(tidyr)

# ── 1. Load data ──────────────────────────────────────────────────────────────
milo <- read.csv("da_results_milo.csv")
prop <- read.csv("cluster_relative_difference_with_pvalues_cortical.csv")

# ── 2. Summarize Milo per cell type ──────────────────────────────────────────
# Split significant neighbourhoods into enriched (KO) and depleted (KO) separately
milo_filt <- milo %>% filter(cell_annotation != "Mixed")

milo_summary <- milo_filt %>%
  group_by(cell_annotation) %>%
  summarise(
    n_nhoods        = n(),
    median_logFC    = median(logFC),
    n_sig_enriched  = sum(SpatialFDR < 0.05 & logFC > 0),  # more in KO
    n_sig_depleted  = sum(SpatialFDR < 0.05 & logFC < 0),  # less in KO
    .groups = "drop"
  ) %>%
  mutate(
    pct_sig_enriched = n_sig_enriched / n_nhoods * 100,
    pct_sig_depleted = n_sig_depleted / n_nhoods * 100,
    pct_sig_total    = (n_sig_enriched + n_sig_depleted) / n_nhoods * 100,
    # A cell type is "bidirectional" if it has significant nhoods in BOTH directions
    bidirectional    = n_sig_enriched > 0 & n_sig_depleted > 0,
    # Net Milo direction for coloring
    milo_direction   = case_when(
      bidirectional                           ~ "Bidirectional",
      n_sig_enriched > 0                      ~ "Enriched in KO",
      n_sig_depleted > 0                      ~ "Depleted in KO",
      TRUE                                    ~ "Not significant"
    )
  )

# ── 3. Merge with proportion data ─────────────────────────────────────────────
combined <- prop %>%
  left_join(milo_summary, by = c("Cluster" = "cell_annotation")) %>%
  mutate(
    prop_sig = adj_p_value < 0.05,
    # Cap extreme proportion values for display
    RelDiff_display = pmax(pmin(RelativeDifference, 250), -100),
    label = ifelse(
      abs(RelativeDifference) > 250,
      paste0(Cluster, "*"),
      Cluster
    )
  )

write.csv(combined, "milo_proportion_combined.csv", row.names = FALSE)

# ── 4. Color palette ──────────────────────────────────────────────────────────
dir_colors <- c(
  "Enriched in KO"  = "#2166AC",
  "Depleted in KO"  = "#D6604D",
  "Bidirectional"   = "#7B3294",
  "Not significant" = "#AAAAAA"
)

# ── PLOT 1: Main scatter — proportion diff vs Milo median logFC ───────────────
# ── Your custom colors and order ──────────────────────────────────────────────
cluster_ids <- levels(Idents(pv_cortical_filtered))

colors <- setNames(
  c(
    "#6B5B95", "#45B8AC", "#955251", "#4E84C4", "#B565A7",
    "#88B04B", "#7B6888", "#C3447A", "#009B77", "#EFC050",
    "#7FCDCD", "#DD4124", "#5B5EA6", "#E07A5F", "#4BACC6",
    "#E8A0BF", "#9B2335", "#C17BAE", "#DECF3F", "#789262",
    "#BC243C"
  ),
  cluster_ids
)

new_order <- c(
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

# Reorder combined$Cluster as a factor so legend follows new_order
combined$Cluster <- factor(combined$Cluster, levels = new_order)

# ── Updated scatter plot ───────────────────────────────────────────────────────
# Label only where both methods agree AND are significant
combined <- combined %>%
  mutate(
    agreement  = (RelativeDifference > 0 & median_logFC > 0) |
      (RelativeDifference < 0 & median_logFC < 0),
    label_text = ifelse(
      prop_sig & milo_direction != "Not significant" & agreement,
      as.character(Cluster),
      NA
    )
  )

# Then in the plot, add after your geom_point layers:
p1 <- ggplot(combined,
             aes(x = RelativeDifference, y = median_logFC,
                 color = Cluster, size = pct_sig_total)) +
  geom_point(data = filter(combined, milo_direction == "Not significant"),
             alpha = 0.4) +
  geom_point(data = filter(combined, milo_direction != "Not significant"),
             alpha = 0.95) +
  geom_point(data = filter(combined, prop_sig),
             aes(x = RelativeDifference, y = median_logFC),
             shape = 21, fill = NA, color = "black", size = 6,
             stroke = 0.9, inherit.aes = FALSE) +
  geom_text_repel(
    aes(label = label_text),
    size              = 3.2,
    na.rm             = TRUE,
    box.padding       = 2,        # more padding = labels pushed further away
    point.padding     = 0.5,        # minimum distance from the dot
    min.segment.length = 0,         # ALWAYS draw connector line, even if label is close
    segment.color     = "grey40",
    segment.size      = 0.4,
    segment.curvature = 0,        # slight curve looks cleaner than straight lines
    force             = 2,          # stronger repulsion between labels
    max.overlaps      = Inf,        # never drop a label
    show.legend       = FALSE,
    color             = "black"
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey80", linewidth = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey80", linewidth = 0.5) +
  # Symlog scale: linear around 0, log in the tails — handles negatives + large values
  scale_x_continuous(
    trans  = scales::pseudo_log_trans(sigma = 10, base = 10),
    breaks = c(-100, -50, -25, 0, 25, 50, 100, 250, 500, 1000, 1500),
    labels = function(x) paste0(x, "%")
  ) +
  scale_color_manual(values = colors, breaks = new_order, name = "Cell type") +
  scale_size_continuous(range = c(3, 11), name = "% sig nhoods") +
  guides(
    color = guide_legend(override.aes = list(size = 4, alpha = 1), ncol = 1),
    size  = guide_legend(ncol = 1)
  ) +
  labs(
    title    = "Proportion analysis vs Milo DA (KO vs WT)",
    subtitle = "X-axis: pseudo-log scale | Bubble size = % significant Milo neighbourhoods | Black ring = significant in proportion test",
    x        = "Relative proportion difference (%) — pseudo-log scale",
    y        = "Milo median logFC (KO vs WT)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "right",
    legend.text     = element_text(size = 8),
    legend.title    = element_text(size = 9, face = "bold"),
    legend.key.size = unit(0.4, "cm"),
    plot.subtitle   = element_text(size = 9, color = "grey40")
  )

p1
ggsave("plot1_scatter.pdf", plot = p1, width = 10, height = 7)


# ── PLOT 2: Diverging bar for bidirectional cell types ────────────────────────
# For ALL cell types with ANY significant neighbourhoods, show + and - wings

sig_celltypes <- milo_summary %>%
  filter(pct_sig_total > 0) %>%
  arrange(desc(pct_sig_enriched - pct_sig_depleted)) %>%
  mutate(cell_annotation = factor(cell_annotation, levels = cell_annotation))

# Reshape to long format for diverging bar
sig_long <- sig_celltypes %>%
  select(cell_annotation, pct_sig_enriched, pct_sig_depleted) %>%
  mutate(pct_sig_depleted = -pct_sig_depleted) %>%  # flip depleted to negative
  pivot_longer(cols = c(pct_sig_enriched, pct_sig_depleted),
               names_to = "direction", values_to = "pct") %>%
  mutate(direction = recode(direction,
                            "pct_sig_enriched" = "Enriched in KO",
                            "pct_sig_depleted" = "Depleted in KO"))

p2 <- ggplot(sig_long, aes(x = cell_annotation, y = pct, fill = direction)) +
  geom_col(width = 0.7, alpha = 0.9) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
  scale_fill_manual(values = c("Enriched in KO" = "#2166AC",
                               "Depleted in KO" = "#D6604D"),
                    name = "") +
  scale_y_continuous(
    labels = function(x) paste0(abs(x), "%"),
    breaks = seq(-100, 100, 20)
  ) +
  coord_flip() +
  labs(
    title    = "Directionality of significant Milo neighbourhoods",
    subtitle = "Cell types with ≥1 significant neighbourhood | + = more in KO, − = less in KO",
    x        = NULL,
    y        = "% of neighbourhoods (SpatialFDR < 0.05)"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom",
        plot.subtitle = element_text(size = 9, color = "grey40"))

p2
ggsave("plot2_diverging.pdf", plot = p2, width = 9, height = 6)


# ── PLOT 3 (bonus): For bidirectional types — logFC distribution strip ─────────
# Shows the full spread of neighbourhood logFCs to reveal sub-population structure

bidir_types <- milo_summary %>%
  filter(bidirectional) %>%
  pull(cell_annotation)

if (length(bidir_types) > 0) {
  
  milo_bidir <- milo_filt %>%
    filter(cell_annotation %in% bidir_types) %>%
    mutate(sig = SpatialFDR < 0.05,
           direction = case_when(
             sig & logFC > 0  ~ "Enriched in KO",
             sig & logFC < 0  ~ "Depleted in KO",
             TRUE             ~ "Not significant"
           ))
  
  p3 <- ggplot(milo_bidir,
               aes(x = logFC, y = cell_annotation,
                   color = direction, alpha = sig, size = sig)) +
    ggbeeswarm::geom_quasirandom(groupOnX = FALSE, bandwidth = 0.8) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    scale_color_manual(values = c("Enriched in KO"  = "#2166AC",
                                  "Depleted in KO"  = "#D6604D",
                                  "Not significant" = "#AAAAAA"),
                       name = "") +
    scale_alpha_manual(values = c("TRUE" = 0.9, "FALSE" = 0.3), guide = "none") +
    scale_size_manual(values  = c("TRUE" = 1.5, "FALSE" = 0.8), guide = "none") +
    labs(
      title    = "Neighbourhood logFC distribution — bidirectional cell types",
      subtitle = "Each dot = one Milo neighbourhood | Colored = SpatialFDR < 0.05",
      x        = "logFC (KO vs WT)", y = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom",
          plot.subtitle = element_text(size = 9, color = "grey40"))
  
  p3
  ggsave("plot3_bidir_strip.pdf", plot = p3, width = 9, height = 4)
}

#### Check Pvalb and Dnmt1 expression in the corticospinal neurons
library(Seurat)
library(ggplot2)
library(patchwork)

# Subset corticospinal neurons
cs_neurons <- subset(pv_cortical_filtered, 
                     cell_annotation %in% c("Corticospinal neurons (Type I)", 
                                            "Corticospinal neurons (Type II)"))

# ── 1. Violin plots: Pvalb and Dnmt1 by cell type and genotype ────────────────
p1 <- VlnPlot(cs_neurons, 
              features  = c("Pvalb", "Dnmt1"),
              group.by  = "cell_annotation",
              split.by  = "Genotype",
              assay     = "RNA",
              pt.size   = 0.1,
              ncol      = 2)

# ── 2. Feature plots on UMAP of the FULL object — are CS neurons tdTom+? ──────
# This shows where corticospinal neurons sit relative to PV cells on the UMAP
p2 <- FeaturePlot(pv_cortical_filtered,
                  features  = c("Pvalb", "Dnmt1"),
                  order     = TRUE,   # plot highest-expressing cells on top
                  ncol      = 2)

# ── 3. Dot plot: clean overview of expression level + % expressing cells ───────
p3 <- DotPlot(cs_neurons,
              features = c("Pvalb", "Dnmt1", 
                           "Fezf2",   # canonical corticospinal marker
                           "Bcl11b",  # layer 5b / corticofugal marker
                           "Crym"),   # another corticospinal marker
              group.by = "cell_annotation",
              split.by = "Genotype",
              assay    = "RNA") +
  RotatedAxis()

# ── 4. Quick summary stats ─────────────────────────────────────────────────────
DefaultAssay(cs_neurons) <- "RNA"

# % of cells expressing each gene per genotype per cluster
for (gene in c("Pvalb", "Dnmt1")) {
  cat("\n──", gene, "──\n")
  expr_mat <- GetAssayData(cs_neurons, layer = "data")[gene, ]
  print(tapply(expr_mat > 0, 
               paste(cs_neurons$cell_annotation, cs_neurons$Genotype, sep = " | "), 
               function(x) round(mean(x) * 100, 1)))
}

# ── 5. Save ───────────────────────────────────────────────────────────────────
ggsave("cs_neurons_violin.pdf",   plot = p1, width = 10, height = 5)
ggsave("cs_neurons_featureplot.pdf", plot = p2, width = 10, height = 5)
ggsave("cs_neurons_dotplot.pdf",  plot = p3, width = 10, height = 5)

# 1. Confirm with a correlation: within corticospinal neurons, 
#    does Dnmt1 expression correlate with any activity-dependent genes?
FeatureScatter(cs_neurons, feature1 = "Dnmt1", feature2 = "Fos",
               group.by = "Genotype")

# 2. Run DE within corticospinal neurons KO vs WT to get the full picture
Idents(cs_neurons) <- cs_neurons$Genotype
de_cs <- FindMarkers(cs_neurons,
                     ident.1   = "PV-Cre/tdTom/Dnmt1 loxP2",
                     ident.2   = "PV-Cre/tdTom",
                     assay     = "RNA",
                     test.use  = "wilcox",
                     min.pct   = 0.1)
head(de_cs, 20)
write.csv(de_cs, "DE_corticospinal_KO_vs_WT.csv")

#### Subclusters analysis --------
# Define the three cell types to investigate
clusters_of_interest <- c("Layer 5/6 IT neurons", "Astrocytes", "Layer 5a IT neurons")

# Run sub-clustering for each and store results in a named list
subclustered <- list()

for (ct in clusters_of_interest) {
  
  cat("\n── Processing:", ct, "──\n")
  
  sub <- subset(pv_cortical_filtered, cell_annotation == ct)
  sub <- RunPCA(sub, npcs = 20, verbose = FALSE)
  sub <- FindNeighbors(sub, dims = 1:15, verbose = FALSE)
  sub <- FindClusters(sub, resolution = 0.3, verbose = FALSE)
  sub <- RunUMAP(sub, dims = 1:15, verbose = FALSE)
  
  # Add a clean name slot for plot titles
  sub$cell_type <- ct
  
  subclustered[[ct]] <- sub
  
  # Print genotype-by-subcluster breakdown
  cat("Subcluster x Genotype table:\n")
  print(table(sub$seurat_clusters, sub$Genotype))
}

# ── Visualize all three side by side ─────────────────────────────────────────
library(patchwork)

for (ct in clusters_of_interest) {
  
  sub <- subclustered[[ct]]
  
  p1 <- DimPlot(sub, group.by = "seurat_clusters", label = TRUE) +
    ggtitle(paste0(ct, " — subclusters")) +
    NoLegend()
  
  p2 <- DimPlot(sub, group.by = "Genotype") +
    ggtitle("Genotype") +
    scale_color_manual(values = c("PV-Cre/tdTom"             = "#4393C3",
                                  "PV-Cre/tdTom/Dnmt1 loxP2" = "#D6604D"))
  
  print(p1 + p2)
  
  # Save per cell type
  ggsave(
    filename = paste0("subcluster_", gsub("/| ", "_", ct), ".pdf"),
    plot     = p1 + p2,
    width    = 12, height = 5
  )
}

# ── Differential expression between sub-clusters per cell type ────────────────
# Only run this after inspecting the UMAPs above to confirm sub-clusters
# separate by genotype

de_results <- list()

for (ct in clusters_of_interest) {
  
  sub <- subclustered[[ct]]
  Idents(sub) <- sub$seurat_clusters
  
  # Find markers distinguishing each sub-cluster
  markers <- FindAllMarkers(
    sub,
    assay      = "RNA",    # change to "SCT" if needed
    only.pos   = FALSE,
    min.pct    = 0.25,
    logfc.threshold = 0.25,
    test.use   = "wilcox",
    verbose    = FALSE
  )
  
  de_results[[ct]] <- markers
  
  # Save DE results as CSV
  write.csv(markers,
            file = paste0("DE_subclusters_", gsub("/| ", "_", ct), ".csv"),
            row.names = TRUE)
  
  cat("\nTop markers for", ct, ":\n")
  print(markers %>% group_by(cluster) %>% slice_max(avg_log2FC, n = 3))
}
