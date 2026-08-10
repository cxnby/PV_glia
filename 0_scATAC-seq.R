options(timeout = 1200)

install.packages("TFMPvalue")
BiocManager::install("EnsDb.Mmusculus.v79")
BiocManager::install("JASPAR2020")
BiocManager::install("DirichletMultinomial")
BiocManager::install("BSgenome")
BiocManager::install("rtracklayer")
BiocManager::install("motifmatchr")
BiocManager::install("TFBSTools")
BiocManager::install("hdf5r")
BiocManager::install("BSgenome.Mmusculus.UCSC.mm10", ask = FALSE)
BiocManager::install("biovizBase")

library(ggplot2)
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(GenomicRanges)
library(EnsDb.Mmusculus.v79)
library(patchwork)
library(dplyr)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Mmusculus.UCSC.mm10)
library(biovizBase)

setwd("~/Downloads/Seurat_scATAC-seq")

###########################

set.seed(42)

# ── 0. Sample mapping from aggregation CSV ───────────────────
agg <- read.csv("aggregation_csv.csv", header = TRUE)
# Row order in agg corresponds to barcode suffixes -1, -2, -3 ...
# The column with sample names is usually called "library_id" or "sample_id"
# Adjust the column name below if yours differs
sample_col <- "library_id"   # <- change if needed
sample_map <- setNames(agg[[sample_col]], seq_len(nrow(agg)))
message("Samples found: ", paste(sample_map, collapse = ", "))

# ── 1. Load count matrix & per-cell metadata ────────────────
counts   <- Read10X_h5("filtered_peak_bc_matrix.h5")
metadata <- read.csv("singlecell.csv", header = TRUE, row.names = 1)

# ── 2. Create ChromatinAssay ─────────────────────────────────
chrom_assay <- CreateChromatinAssay(
  counts    = counts,
  sep       = c(":", "-"),          # peak name format  chr:start-end
  fragments = "fragments.tsv.gz",   # .tbi index must be in same folder
  min.cells = 10,
  min.features = 200
)

obj <- CreateSeuratObject(
  counts    = chrom_assay,
  assay     = "peaks",
  meta.data = metadata
)

# ── 3. Add sample-of-origin from barcode suffix ──────────────
# Barcodes look like ACGTACGT-1 / ACGTACGT-2 etc.
barcode_suffix      <- as.integer(sub(".*-", "", colnames(obj)))
obj$sample_number   <- barcode_suffix
obj$sample        <- unname(sample_map[obj$sample_number])   # human-readable name
message("Cell counts per sample:")
print(table(obj$sample))

# ── 4. Gene annotations ──────────────────────────────────────
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevels(annotations) <- paste0("chr", seqlevels(annotations))   # add "chr" prefix
genome(annotations)     <- "mm10"                                  # adjust if human: hg38
Annotation(obj)         <- annotations

# ── 5. QC metrics ────────────────────────────────────────────
obj <- NucleosomeSignal(obj)
obj <- TSSEnrichment(obj, fast = FALSE)

# Visualise QC (optional — saves to file)
qc_plot <- VlnPlot(
  obj,
  features  = c("nCount_peaks", "TSS.enrichment", "nucleosome_signal"),
  group.by  = "sample",
  pt.size   = 0,
  ncol      = 3
)
ggsave("qc_violin.pdf", qc_plot, width = 12, height = 4)

obj_raw <- obj
saveRDS(obj_raw, "scATAC_00_preQC.rds")

# ── 6. Filter low-quality cells ──────────────────────────────
# Inspect your actual distributions
summary(obj@meta.data[, c("nCount_peaks", "TSS.enrichment", "nucleosome_signal")])

VlnPlot(
  obj,
  features = c("nCount_peaks", "TSS.enrichment", "nucleosome_signal"),
  group.by = "sample",
  pt.size  = 0,
  ncol     = 3,
  log      = TRUE   # log scale helps for nCount_peaks
)

obj <- subset(
  obj,
  subset =
    nCount_peaks       > 1000  &
    nCount_peaks       < 50000 &
    TSS.enrichment     > 1     &
    nucleosome_signal  < 4
)
message("Cells after QC: ", ncol(obj))

# ── 7. Normalisation → Feature selection → LSA ───────────────
obj <- RunTFIDF(obj)
obj <- FindTopFeatures(obj, min.cutoff = "q5")   # top 95% most variable peaks
obj <- RunSVD(obj)                               # LSA / dimensionality reduction

# Check which LSI components correlate with sequencing depth
# (component 1 almost always does — exclude it)
depth_cor <- DepthCor(obj)
print(depth_cor)
# Rule of thumb: drop dims that show |correlation| > 0.75 with depth
# Default: exclude dim 1 (almost -1) and dim 2 (-0.58 which is unusually high)

# ── 8. Clustering + UMAP ─────────────────────────────────────
obj <- RunUMAP(obj, reduction = "lsi", dims = 3:30)
obj <- FindNeighbors(obj, reduction = "lsi", dims = 3:30)
obj <- FindClusters(obj, algorithm = 3, resolution = 0.5)   # SLM algorithm

# ── 9. Quick plots ───────────────────────────────────────────
p1 <- DimPlot(obj, group.by = "seurat_clusters", label = TRUE, pt.size = 0.1) +
  ggtitle("Clusters")

p2 <- DimPlot(obj, group.by = "sample", pt.size = 0.1) +
  ggtitle("Sample of origin")

ggsave("umap_clusters.pdf", p1, width = 12, height = 10)
ggsave("umap_samples.pdf",  p2, width = 12, height = 10)

# ── 10. Save object ──────────────────────────────────────────
saveRDS(obj, "scATAC_01_processed.rds")

message("Done! Object saved to scATAC_processed.rds")
message("Key metadata columns: seurat_clusters | sample | sample_number")
message("  obj$sample            -> sample name")
message("  obj$seurat_clusters   -> cluster identity")
message("  Idents(obj)           -> active identity (clusters by default)")

# 11. Read the RNA data ---------------
rna <- readRDS("rna_for_atac.rds")

#################### ATAC-only analysis ##################
# 1a. Gene activity matrix ────────────────────────────────────
# Summarises chromatin accessibility around each gene body as a
# proxy for gene expression — needed for label transfer
gene_activities <- GeneActivity(obj)

obj[["RNA"]] <- CreateAssayObject(counts = gene_activities)
obj <- NormalizeData(
  obj,
  assay          = "RNA",
  normalization.method = "LogNormalize",
  scale.factor   = median(obj$nCount_peaks)
)

# Quick sanity check: known cortical marker genes
DefaultAssay(obj) <- "RNA"
marker_plot <- FeaturePlot(
  obj,
  features  = c("Slc17a7",  # excitatory neurons
                "Gad2",     # inhibitory neurons
                "Mbp",      # oligodendrocytes
                "S100b",     # astrocytes
                "P2ry12",   # microglia
                "Pdgfra",   # OPCs
                "Vwf"),  # endothelial
  pt.size   = 0.1,
  ncol      = 3,
  max.cutoff = "q95"
)
ggsave("gene_activity_markers.pdf", marker_plot, width = 14, height = 9)

# ── PART 2: LABEL TRANSFER FROM scRNA-seq ────────────────────

# ══════════════════════════════════════════════════════════════
# scATAC-seq ANALYSIS — PART 2: LABEL TRANSFER FROM scRNA-seq
#
# HOW TO USE THIS SCRIPT:
#   Run section by section. Each CHECKPOINT block will print
#   diagnostic output — read it before continuing. Instructions
#   for what to look for and how to adjust are in the comments
#   immediately below each checkpoint.
# ══════════════════════════════════════════════════════════════

DefaultAssay(obj) <- "peaks"

# ──────────────────────────────────────────────────────────────
# 2a. Prepare RNA reference
# ──────────────────────────────────────────────────────────────
# Change "cell_type" if your RNA annotation column is named differently
# e.g. rna$seurat_clusters, rna$celltype_annotation, rna$CellType
rna$cell_type      <- as.character(Idents(rna))
rna_annotation_col <- "cell_type"

rna <- FindVariableFeatures(rna, nfeatures = 5000)

# Find common genes between gene activity matrix and RNA variable features
DefaultAssay(obj) <- "RNA"
common_features   <- intersect(
  rownames(obj[["RNA"]]),
  VariableFeatures(rna)
)

# ── CHECKPOINT 1: Common features ─────────────────────────────
message("\n══ CHECKPOINT 1: Common features for label transfer ══")
message("Common features found: ", length(common_features))
message("RNA variable features: ", length(VariableFeatures(rna)))
message("Gene activity genes:   ", nrow(obj[["RNA"]]))

# WHAT TO LOOK FOR:
#   Good:    500–2000 common features
#   Warning: 100–499  — transfer will work but may be noisy
#   Problem: <100     — something is wrong with gene naming (see below)
#
# IF TOO FEW:
#   Check whether gene names match (case, prefix, version number):
#     head(rownames(obj[["RNA"]]))   # e.g. "Gfap", "Snap25"
#     head(VariableFeatures(rna))    # should match format exactly
#   If RNA uses Ensembl IDs and ATAC uses symbols (or vice versa),
#   you need to convert one to match the other before proceeding.
#
# STOP HERE IF < 100 common features. Fix gene name matching first.

if (length(common_features) < 100) {
  stop("Too few common features (", length(common_features), "). ",
       "Check gene naming format between obj[['RNA']] and rna. ",
       "Run: head(rownames(obj[['RNA']])) and head(VariableFeatures(rna))")
}

obj <- ScaleData(obj, features = common_features)

# ──────────────────────────────────────────────────────────────
# 2b. Find anchors
# ──────────────────────────────────────────────────────────────
# k.anchor = 10 (default 5): higher = more anchors found for rare cell types
# k.filter = 100 (default 200): lower = less strict, helps rare types form anchors
# Adjust these if CHECKPOINT 2 shows too few anchors for specific cell types

transfer_anchors <- FindTransferAnchors(
  reference       = rna,
  query           = obj,
  features        = common_features,
  reference.assay = "RNA",
  query.assay     = "RNA",
  reduction       = "cca",
  dims            = 1:30,
  k.anchor        = 10,    # increase to 15–20 if rare types are missed
  k.filter        = 100,   # decrease to 50 if rare types still not anchored
  k.score         = 30
)

# ── CHECKPOINT 2: Anchor quality ──────────────────────────────
message("\n══ CHECKPOINT 2: Transfer anchors ══")
message("Total anchors found: ", nrow(transfer_anchors@anchors))

# WHAT TO LOOK FOR:
#   Good:    > 500 anchors
#   Warning: 100–500 — will work, monitor prediction scores
#   Problem: < 100  — increase k.anchor, decrease k.filter, or check
#                     that CCA is appropriate for your data
#
# IF TOO FEW:
#   Try: k.anchor = 20, k.filter = 50
#   Or switch reduction = "rpca" if CCA is too slow / not converging

# ──────────────────────────────────────────────────────────────
# 2c. Transfer labels
# ──────────────────────────────────────────────────────────────
predicted_labels <- TransferData(
  anchorset        = transfer_anchors,
  refdata          = rna$cell_type,
  weight.reduction = obj[["lsi"]],
  dims             = 3:30
)

DefaultAssay(obj) <- "peaks"
obj <- AddMetaData(obj, metadata = predicted_labels)

# ── CHECKPOINT 3: Per-cell prediction quality ─────────────────
message("\n══ CHECKPOINT 3: Per-cell prediction scores ══")

score_summary <- summary(obj$prediction.score.max)
print(score_summary)
message("\nConfidence breakdown:")
message("  High confidence (>= 0.75): ",
        sum(obj$prediction.score.max >= 0.75), " cells (",
        round(mean(obj$prediction.score.max >= 0.75) * 100, 1), "%)")
message("  Medium (0.50–0.75):        ",
        sum(obj$prediction.score.max >= 0.50 & obj$prediction.score.max < 0.75),
        " cells")
message("  Low confidence (< 0.50):   ",
        sum(obj$prediction.score.max < 0.50), " cells (",
        round(mean(obj$prediction.score.max < 0.50) * 100, 1), "%)")

# WHAT TO LOOK FOR:
#   Good:    Median score > 0.7, <20% low-confidence cells
#   Warning: Median 0.5–0.7 — transfer is uncertain, inspect per-cluster
#   Problem: Median < 0.5  — label transfer is not reliable overall;
#            consider using more variable features (nfeatures = 7000) or
#            a different reduction method
#
# FLAG low-confidence cells so you can track them downstream
obj$prediction_quality <- ifelse(
  obj$prediction.score.max >= 0.5, "high", "low"
)
print(table(obj$prediction_quality))

# DIAGNOSTIC PLOTS — run all of these before continuing
# Plot 1: Prediction score on UMAP — red = unreliable
# Look for: any cluster that is mostly red → its label will be wrong
p1 <- FeaturePlot(obj, features = "prediction.score.max", pt.size = 0.1) +
  scale_color_gradient(low = "lightblue", high = "darkblue") +
  ggtitle("Label transfer confidence")
print(p1)

ggsave("label_transfer_confidence.pdf", p1, width = 10, height = 8)

# Plot 2: Per-cell predicted labels before majority vote
# Look for: clusters that show a MIX of colors → majority vote will
# pick the dominant one even if 60% of cells belong to a minority type
p2 <- DimPlot(obj, group.by = "predicted.id",
              label = TRUE, repel = TRUE, pt.size = 0.1) +
  NoLegend() +
  ggtitle("Per-cell predictions (pre-majority-vote) — inspect carefully")
print(p2)

ggsave("per-cell_predictions.pdf", p2, width = 10, height = 8)


# STOP HERE — examine p1 and p2 before continuing.
# Key question: Does any clearly separated ATAC cluster contain mostly
# one color in p2 but surrounded by another cell type?
# If yes: that cluster is at risk of being absorbed — note its cluster
# number so you can manually override it in section 2e below.

saveRDS(obj, "scATAC_02_labeled.rds")
# obj <- readRDS("scATAC_02_labeled.rds")

# ──────────────────────────────────────────────────────────────
# 2d. Majority-vote annotation (score-weighted)
# ──────────────────────────────────────────────────────────────
# Using SCORE-WEIGHTED vote: cells with high confidence count more.
# This protects small clusters from being outvoted by a dominant type.

cluster_vote_summary <- obj@meta.data %>%
  group_by(seurat_clusters, predicted.id) %>%
  summarise(
    n            = n(),
    mean_score   = round(mean(prediction.score.max), 3),
    weighted_n   = sum(prediction.score.max),   # key: weight by confidence
    .groups      = "drop"
  ) %>%
  group_by(seurat_clusters) %>%
  mutate(
    pct_of_cluster = round(n / sum(n) * 100, 1)
  ) %>%
  arrange(seurat_clusters, desc(weighted_n))

# ── CHECKPOINT 4: Cluster composition ─────────────────────────
message("\n══ CHECKPOINT 4: Cluster composition (top 2 candidates per cluster) ══")

# Shows top 2 competing labels per cluster with their weighted scores
top2_per_cluster <- cluster_vote_summary %>%
  group_by(seurat_clusters) %>%
  slice_max(weighted_n, n = 2, with_ties = FALSE)

print(as.data.frame(top2_per_cluster), row.names = FALSE)

# WHAT TO LOOK FOR:
#   Healthy cluster: top label has >> 50% of cells, mean_score > 0.6
#   PROBLEM cluster: top two labels are close in weighted_n AND/OR
#                    pct_of_cluster is split (e.g. 55% vs 40%)
#                    → this cluster is a candidate for wrong annotation
#
# Also check median prediction score per cluster:
message("\n── Median prediction score per cluster ──")
cluster_score_check <- obj@meta.data %>%
  group_by(seurat_clusters) %>%
  summarise(
    median_score = round(median(prediction.score.max), 3),
    n_cells      = n(),
    top_label    = names(sort(table(predicted.id), decreasing = TRUE))[1],
    top_label_pct = round(max(table(predicted.id)) / n() * 100, 1)
  ) %>%
  arrange(median_score)

print(as.data.frame(cluster_score_check), row.names = FALSE)

# WHAT TO LOOK FOR:
#   median_score < 0.4 → label transfer failed for this cluster;
#                        use manual override in 2e
#   top_label_pct < 50 → no clear winner; inspect UMAP for this cluster

# Assign best label using weighted vote
best_labels <- cluster_vote_summary %>%
  group_by(seurat_clusters) %>%
  slice_max(weighted_n, n = 1, with_ties = FALSE) %>%
  select(seurat_clusters, predicted.id, pct_of_cluster, mean_score)

label_map <- setNames(best_labels$predicted.id,
                      best_labels$seurat_clusters)

obj$cell_type <- unname(label_map[as.character(obj$seurat_clusters)])

# ──────────────────────────────────────────────────────────────
# 2e. Manual overrides — EDIT THIS SECTION based on Checkpoint 4
# ──────────────────────────────────────────────────────────────
# If Checkpoint 4 shows a cluster with wrong or uncertain annotation,
# override it here. Comment out any lines that don't apply to your data.
#
# HOW TO DECIDE WHAT TO PUT HERE:
#   1. Note which cluster numbers had low median_score or split vote
#   2. Look at that cluster in p1 (confidence) and p2 (per-cell labels)
#   3. Look at gene activity (FeaturePlot) of known marker genes for
#      the cell type you expect it to be:
#        FeaturePlot(obj, features = c("Rbfox3", "Gfap", "Olig2"), pt.size = 0.1)
#   4. Assign the biologically correct label below

# Run this to visualize cluster 18 specifically
DimPlot(obj, cells.highlight = WhichCells(obj, idents = NULL,
                                          expression = seurat_clusters == "18"), pt.size = 0.1) +
  ggtitle("Cluster 18 location")

# ── Step 1: Remove cluster 18 (contamination) ─────────────────
# Cluster 18 is scattered across clusters (0, 1, 3, 5, 6, 9, 11, 14, 20) with
# mean_score = 0.457 (below threshold). These are likely doublets
# or ambient RNA contamination. Remove them entirely.
cells_to_remove <- WhichCells(obj, expression = seurat_clusters == 18)
message("Removing ", length(cells_to_remove), " cells from cluster 18 (contamination)")
obj <- obj[, !colnames(obj) %in% cells_to_remove]

# ── Step 2: Recover Layer 2/3 IT neurons ──────────────────────
# These cells are split across clusters 0 (Layer 5a IT) and 5
# (Layer 4 sensory). Rather than losing them to majority vote,
# we recover individual cells with high prediction confidence.
#
# TUNING: Adjust the score threshold (currently 0.75) if you
# recover too few or too many cells:
#   - Too few recovered: lower to 0.65
#   - Too many (biologically implausible): raise to 0.85
#

# Diagnostic — see the score distribution for Layer 2/3 IT cells in cluster 0
scores_c0_l23 <- obj@meta.data[
  obj$seurat_clusters == "0" & obj$predicted.id == "Layer 2/3 IT neurons",
  "prediction.score.max"
]
summary(scores_c0_l23)
hist(scores_c0_l23, breaks = 20,
     main = "Layer 2/3 IT scores in cluster 0",
     xlab = "prediction.score.max")

# Unimodal distribution -> do not subcluster here

# Check cluster 5 score distribution for Layer 2/3 IT cells
scores_c5_l23 <- obj@meta.data[
  obj$seurat_clusters == "5" & obj$predicted.id == "Layer 2/3 IT neurons",
  "prediction.score.max"
]
summary(scores_c5_l23)
hist(scores_c5_l23, breaks = 20,
     main = "Layer 2/3 IT scores in cluster 5",
     xlab = "prediction.score.max")

### Valley sits around 0.55, only subcluster for cluster 5

layer23_recovery_threshold <- 0.55

layer23_cells <- rownames(obj@meta.data[
  obj$seurat_clusters %in% c("5") &
    obj$predicted.id == "Layer 2/3 IT neurons" &
    obj$prediction.score.max >= layer23_recovery_threshold,
])

obj$cell_type[layer23_cells] <- "Layer 2/3 IT neurons"

message("Recovered Layer 2/3 IT neurons: ", length(layer23_cells), " cells")
message("  From cluster 5: ",
        sum(obj$seurat_clusters[layer23_cells] == "5"), " cells")

# ── Recover Corticospinal neurons (Type I) ────────────────────
# Type I was split across clusters 6 and 20 and outvoted in both.
# User confirmed they are spatially distinct on UMAP from Type II.
# Using a lower threshold (0.65) because this is a rare population
# with fewer cells to draw anchors from during transfer.
#
# TUNING: If recovered cells look spatially wrong on UMAP, raise to 0.70
#         If you recover 0 cells, lower to 0.60

typeI_recovery_threshold <- 0.65

typeI_cells <- rownames(obj@meta.data[
  obj$predicted.id    == "Corticospinal neurons (Type I)" &
    obj$prediction.score.max >= typeI_recovery_threshold,
])

obj$cell_type[typeI_cells] <- "Corticospinal neurons (Type I)"

message("Recovered Corticospinal neurons (Type I): ", length(typeI_cells), " cells")
message("  From cluster 6:  ", sum(obj$seurat_clusters[typeI_cells] == "6"),  " cells")
message("  From cluster 20: ", sum(obj$seurat_clusters[typeI_cells] == "20"), " cells")
message("  From other clusters: ",
        sum(!obj$seurat_clusters[typeI_cells] %in% c("6", "20")), " cells")

# ── Step 3: Explicit cluster-level overrides ──────────────────

# Cluster 2: Oligodendrocyte states — ATAC cannot distinguish maturation
# stages reliably (mean_score only 0.563). Collapse to one label.
obj$cell_type[obj$seurat_clusters == "2"] <- "Oligodendrocytes"

# Cluster 15: Leptomeningeal + meningeal fibroblasts — biologically
# appropriate to collapse since ATAC cannot resolve transcriptional
# fibroblast subtypes. Use the broader term.
obj$cell_type[obj$seurat_clusters == "15"] <- "Meningeal fibroblasts"

# Cluster 19: Layer 6b neurons — spatially distinct from cluster 3
# (Layer 6a CT neurons) on UMAP. Keep as separate population.
obj$cell_type[obj$seurat_clusters == "19"] <- "Layer 6b neurons"

# Cluster 20: Corticospinal Type I and II are spatially distinct on UMAP.
# Do NOT collapse — keep them separate. The majority vote already assigns
# Type II here (77.3%), which is correct per the UMAP position.
# No action needed.

# ── Step 4: Verification ──────────────────────────────────────
message("\n══ CHECKPOINT 5: Final cell type composition ══")
ct_table <- sort(table(obj$cell_type), decreasing = TRUE)
print(ct_table)

# WHAT TO LOOK FOR:
#   - "Layer 2/3 IT neurons" should now appear with recovered cell count
#   - "Meningeal fibroblasts" should contain both former cluster 15 types
#   - "Oligodendrocytes" should be the merged cluster 2
#   - Both "Corticospinal neurons (Type I)" and "(Type II)" appear separately
#   - Cluster 18 cells are gone (total cell count should be ~23k not ~23.4k)
#   - No NA values

if (any(is.na(obj$cell_type))) {
  warning("NA cell_type values found for clusters: ",
          paste(unique(obj$seurat_clusters[is.na(obj$cell_type)]),
                collapse = ", "))
}

# ── Step 5: Set Idents and final plots ────────────────────────
Idents(obj) <- "cell_type"

p_final <- DimPlot(obj, label = TRUE, repel = TRUE, pt.size = 0.1) +
  NoLegend() +
  ggtitle("Final annotation — post-override")
print(p_final)

# FINAL SANITY CHECKS:
#   1. Layer 2/3 IT neurons form a visible group (may straddle clusters 0/5)
#   2. Layer 6b neurons are clearly separate from cluster 3 (Layer 6a CT)
#   3. Both Corticospinal subtypes appear as distinct labelled groups
#   4. Cluster 18 area (was scattered across 0/1/5/20) shows no stray label
#   5. No isolated "island" of cells with a mismatched label
#
# IF Layer 2/3 IT recovery looks too sparse or fragmented:
#   Lower threshold: layer23_recovery_threshold <- 0.65
#   Then re-run Steps 3 and 6 only (no need to re-run transfer)
#
# IF you see biologically wrong cells labelled as Layer 2/3 IT:
#   Raise threshold: layer23_recovery_threshold <- 0.85

saveRDS(obj, "scATAC_03_labeled_clean.rds")
message("\n=== Part 2 complete. Saved: scATAC_03_labeled_clean.rds ===")


# ── PART 3: VISUALIZATION ────────────────────────────────────

# Original RNA label order (for predicted.id and RNA UMAP)
original_new_order <- c(
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

# ATAC label order (with renames)
atac_new_order <- c(
  "Layer 2/3 IT neurons",
  "Layer 4 sensory neurons",
  "Layer 5a IT neurons",
  "Layer 5b PT neurons",
  "Layer 5/6 IT neurons",
  "Layer 6a corticothalamic neurons",
  "Layer 6b neurons",
  "Deep-layer extratelencephalic neurons",
  "Corticospinal neurons (Type I)",
  "Corticospinal neurons (Type II)",               # renamed from Type II
  "PV+ interneurons",
  "SST+ interneurons",
  "VIP+ interneurons",
  "Astrocytes",
  "Oligodendrocyte precursor cells",
  "Oligodendrocytes",             # renamed, gets myelinating color
  "Perineuronal oligodendrocytes",
  "Myelinating oligodendrocytes",
  "Microglia",
  "Endothelial cells",
  "Leptomeningeal cells",
  "Meningeal fibroblasts"
)

# Shared hex colors — same position = same color in both palettes
hex_colors <- c(
  "#6B5B95", "#45B8AC", "#955251", "#4E84C4", "#B565A7",
  "#88B04B", "#7B6888", "#C3447A", "#009B77", "#EFC050",
  "#7FCDCD", "#DD4124", "#5B5EA6", "#E07A5F", "#4BACC6",
  "#E8A0BF", "#ea8a33", "#9B2335", "#C17BAE", "#DECF3F",
  "#789262", "#BC243C"
)

rna_colors_original <- setNames(hex_colors, original_new_order)
rna_colors_atac     <- setNames(hex_colors, atac_new_order)

# predicted.id colors (original RNA names)
predicted_types  <- original_new_order[original_new_order %in% unique(as.character(obj$predicted.id))]
obj$predicted.id <- factor(obj$predicted.id, levels = predicted_types)
predicted_colors <- rna_colors_original[predicted_types]

# cell_type colors (ATAC renamed names)
atac_cell_types <- atac_new_order[atac_new_order %in% unique(as.character(obj$cell_type))]
obj$cell_type   <- factor(obj$cell_type, levels = atac_cell_types)
Idents(obj)     <- "cell_type"
atac_colors     <- rna_colors_atac[atac_cell_types]

# Sanity check — will error immediately if any NAs sneak in
stopifnot(!any(is.na(atac_colors)))
stopifnot(!any(is.na(predicted_colors)))

# UMAP plots
p1 <- DimPlot(obj, group.by = "seurat_clusters",
              label = TRUE, pt.size = 0.1) +
  ggtitle("ATAC clusters") + NoLegend()

p2 <- DimPlot(obj,
              group.by = "predicted.id",
              cols     = predicted_colors,
              label    = TRUE, pt.size = 0.1, repel = TRUE) +
  ggtitle("Predicted cell type (label transfer)") + NoLegend()

p3 <- DimPlot(obj,
              group.by = "cell_type",
              cols     = atac_colors,
              label    = TRUE, pt.size = 0.1, repel = TRUE) +
  ggtitle("Cluster annotation") + NoLegend()

p4 <- DimPlot(obj, group.by = "sample",
              pt.size = 0.1) +
  ggtitle("Sample of origin")

p5 <- FeaturePlot(obj, features = "prediction.score.max",
                  pt.size = 0.1) +
  ggtitle("Label transfer confidence")

layout <- (p1 | p2) / (p3 | p4) / p5 +
  plot_layout(heights = c(1, 1, 1.5))
ggsave("integration_overview.pdf", layout, width = 12, height = 20)

# Side-by-side RNA vs ATAC UMAP
# RNA UMAP — no labels, legend kept for collection

rna_types  <- original_new_order[original_new_order %in% unique(as.character(rna$cell_type))]
rna$cell_type <- factor(rna$cell_type, levels = rna_types)
rna_colors <- rna_colors_original[rna_types]

p_rna <- DimPlot(rna,
                 group.by = "cell_type",
                 cols     = rna_colors,
                 pt.size  = 0.1,
                 label    = FALSE) +
  ggtitle("scRNA-seq") +
  theme(
    legend.text     = element_text(size = 9),
    legend.key.size = unit(0.4, "cm")
  ) +
  guides(color = guide_legend(
    ncol     = 1,
    override.aes = list(size = 3)
  ))

# ATAC UMAP — no labels, legend kept for collection
p_atac <- DimPlot(obj,
                  group.by = "cell_type",
                  cols     = atac_colors,
                  pt.size  = 0.1,
                  label    = FALSE) +
  ggtitle("scATAC-seq (projected)") +
  theme(
    legend.text     = element_text(size = 9),
    legend.key.size = unit(0.4, "cm")
  ) +
  guides(color = guide_legend(
    ncol     = 1,
    override.aes = list(size = 3)
  ))

# Combine with patchwork — collect_guides merges the two legends into one
# on the far right. Because both use the same color→cell type mapping,
# patchwork will deduplicate them into a single unified legend.
combined <- (p_rna | p_atac)

ggsave("rna_vs_atac_umap.pdf", combined, width = 20, height = 7)

#### DIFFERENTIAL ACCESSIBILITY #####
# ── PART 4: DIFFERENTIAL ACCESSIBILITY ───────────────────────
# Five contrasts per cell type:
#   A. KO vs Ctrl — all sexes combined (sex regressed out)
#   B. KO vs Ctrl — females only
#   C. KO vs Ctrl — males only
#   D. Ctrl female vs Ctrl male
#   E. KO female vs KO male

# ADD SEX & GENOTYPE METADATA ──────────────────────
# Must run before anything else.
# Sample names are expected to contain "Ctrl" or "KO" and end in "_f" or "_m"
# e.g. "Ctrl_1_f", "KO_2_m", "Ctrl_3_m"

obj$genotype <- ifelse(grepl("Ctrl", obj$sample), "Ctrl", "KO")
obj$sex      <- ifelse(grepl("_f",   obj$sample), "female", "male")

# Sanity check — should show 3 Ctrl + 3 KO per sex
message("=== Sex x Genotype cell counts ===")
print(table(obj$genotype, obj$sex))
stopifnot("genotype" %in% colnames(obj@meta.data))
stopifnot("sex"      %in% colnames(obj@meta.data))

# RUN THE DIFFERENTIAL ANALYSIS -------------

DefaultAssay(obj) <- "peaks"
cell_types        <- levels(obj$cell_type)
dir.create("DA_results", showWarnings = FALSE)

# Helper — LR test with field-standard thresholds
run_da_fast <- function(obj, cells_1, cells_2, label, regress_sex = FALSE) {
  n1 <- length(cells_1); n2 <- length(cells_2)
  if (n1 < 20 | n2 < 20) {
    message("  -> SKIP [", label, "] — too few cells (", n1, " vs ", n2, ")")
    return(NULL)
  }
  lv <- if (regress_sex) c("nCount_peaks", "sex") else "nCount_peaks"
  tryCatch({
    FindMarkers(
      obj,
      ident.1         = cells_1,
      ident.2         = cells_2,
      min.pct         = 0.05,   # Field standard for scATAC-seq
      logfc.threshold = 0.05,   # Field standard for scATAC-seq
      test.use        = "LR",
      latent.vars     = lv
    )
  }, error = function(e) {
    message("  -> ERROR [", label, "]: ", e$message)
    NULL
  })
}

md <- obj@meta.data

# ──────────────────────────────────────────────────────────────
# CONTRAST 1: KO vs Ctrl — All Sexes (sex regressed out)
# ──────────────────────────────────────────────────────────────
message("\n=== CONTRAST 1: KO vs Ctrl (All) ===")
for (ct in cell_types) {
  ct_safe <- gsub("[^A-Za-z0-9]", "_", ct)
  outfile <- file.path("DA_results", paste0(ct_safe, "__KOvCtrl_all.csv"))
  if (file.exists(outfile)) { message("Skipping ", ct, " (done)"); next }
  message("Processing: ", ct)
  res <- run_da_fast(obj,
                     rownames(md[md$cell_type == ct & md$genotype == "KO",   ]),
                     rownames(md[md$cell_type == ct & md$genotype == "Ctrl", ]),
                     paste(ct, "KO vs Ctrl all"), regress_sex = TRUE)
  if (!is.null(res)) write.csv(res, outfile)
}

# ──────────────────────────────────────────────────────────────
# CONTRAST 3: KO vs Ctrl — Males Only
# ──────────────────────────────────────────────────────────────
message("\n=== CONTRAST 3: KO vs Ctrl (Males) ===")
for (ct in cell_types) {
  ct_safe <- gsub("[^A-Za-z0-9]", "_", ct)
  outfile <- file.path("DA_results", paste0(ct_safe, "__KOvCtrl_male.csv"))
  if (file.exists(outfile)) { message("Skipping ", ct, " (done)"); next }
  message("Processing: ", ct)
  res <- run_da_fast(obj,
                     rownames(md[md$cell_type == ct & md$genotype == "KO"   & md$sex == "male", ]),
                     rownames(md[md$cell_type == ct & md$genotype == "Ctrl" & md$sex == "male", ]),
                     paste(ct, "KO vs Ctrl male"), regress_sex = FALSE)
  if (!is.null(res)) write.csv(res, outfile)
}

# ──────────────────────────────────────────────────────────────
# CONTRAST 2: KO vs Ctrl — Females Only
# ──────────────────────────────────────────────────────────────
message("\n=== CONTRAST 2: KO vs Ctrl (Females) ===")
for (ct in cell_types) {
  ct_safe <- gsub("[^A-Za-z0-9]", "_", ct)
  outfile <- file.path("DA_results", paste0(ct_safe, "__KOvCtrl_female.csv"))
  if (file.exists(outfile)) { message("Skipping ", ct, " (done)"); next }
  message("Processing: ", ct)
  res <- run_da_fast(obj,
                     rownames(md[md$cell_type == ct & md$genotype == "KO"   & md$sex == "female", ]),
                     rownames(md[md$cell_type == ct & md$genotype == "Ctrl" & md$sex == "female", ]),
                     paste(ct, "KO vs Ctrl female"), regress_sex = FALSE)
  if (!is.null(res)) write.csv(res, outfile)
}



# ──────────────────────────────────────────────────────────────
# CONTRAST 4: Ctrl Female vs Ctrl Male
# ──────────────────────────────────────────────────────────────
message("\n=== CONTRAST 4: Ctrl Female vs Male ===")
for (ct in cell_types) {
  ct_safe <- gsub("[^A-Za-z0-9]", "_", ct)
  outfile <- file.path("DA_results", paste0(ct_safe, "__Ctrl_FvM.csv"))
  if (file.exists(outfile)) { message("Skipping ", ct, " (done)"); next }
  message("Processing: ", ct)
  res <- run_da_fast(obj,
                     rownames(md[md$cell_type == ct & md$genotype == "Ctrl" & md$sex == "female", ]),
                     rownames(md[md$cell_type == ct & md$genotype == "Ctrl" & md$sex == "male",   ]),
                     paste(ct, "Ctrl F vs M"), regress_sex = FALSE)
  if (!is.null(res)) write.csv(res, outfile)
}

# ──────────────────────────────────────────────────────────────
# CONTRAST 5: KO Female vs KO Male
# ──────────────────────────────────────────────────────────────
message("\n=== CONTRAST 5: KO Female vs Male ===")
for (ct in cell_types) {
  ct_safe <- gsub("[^A-Za-z0-9]", "_", ct)
  outfile <- file.path("DA_results", paste0(ct_safe, "__KO_FvM.csv"))
  if (file.exists(outfile)) { message("Skipping ", ct, " (done)"); next }
  message("Processing: ", ct)
  res <- run_da_fast(obj,
                     rownames(md[md$cell_type == ct & md$genotype == "KO" & md$sex == "female", ]),
                     rownames(md[md$cell_type == ct & md$genotype == "KO" & md$sex == "male",   ]),
                     paste(ct, "KO F vs M"), regress_sex = FALSE)
  if (!is.null(res)) write.csv(res, outfile)
}

message("\n=== ALL CONTRASTS COMPLETE ===")

# ══════════════════════════════════════════════════════════════
# POST-PROCESSING: Replace Bonferroni with BH FDR correction
# Bonferroni in FindMarkers uses ALL peaks as denominator (~177k)
# but only ~28k were actually tested — this overcorrects ~6x.
# BH correction on the tested peaks is the appropriate fix.
# ══════════════════════════════════════════════════════════════

library(dplyr)

da_files <- list.files("DA_results", pattern = "\\.csv$", full.names = TRUE)

for (f in da_files) {
  df <- read.csv(f, row.names = 1)
  
  # Add BH-corrected FDR column based only on peaks actually tested
  df$p_val_fdr <- p.adjust(df$p_val, method = "BH")
  
  # Write back to same file (overwrites with extra column)
  write.csv(df, f)
}

message("BH FDR correction applied to ", length(da_files), " files")

# Quick summary of how many significant peaks per contrast
summary_df <- lapply(da_files, function(f) {
  df  <- read.csv(f, row.names = 1)
  name <- gsub("DA_results/|.csv", "", f)
  data.frame(
    contrast    = name,
    n_tested    = nrow(df),
    sig_bonf    = sum(df$p_val_adj  < 0.05, na.rm = TRUE),
    sig_fdr05   = sum(df$p_val_fdr  < 0.05, na.rm = TRUE),
    sig_fdr10   = sum(df$p_val_fdr  < 0.10, na.rm = TRUE)
  )
}) |> bind_rows()

print(summary_df)
write.csv(summary_df, "DA_results/significance_summary.csv", row.names = FALSE)

## Plot interesting findings
# ══════════════════════════════════════════════════════════════
# SHARED SETUP — colors & order (used by both plots)
# ══════════════════════════════════════════════════════════════

library(ggplot2)
library(ggrepel)
library(dplyr)
library(tidyr)
library(stringr)
library(forcats)
library(scales)
library(ragg)   # for PNG — install.packages("ragg") if missing

atac_new_order <- c(
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
  "Oligodendrocytes",
  "Perineuronal oligodendrocytes",
  "Myelinating oligodendrocytes",
  "Microglia",
  "Endothelial cells",
  "Leptomeningeal cells",
  "Meningeal fibroblasts"
)
hex_colors <- c(
  "#6B5B95","#45B8AC","#955251","#4E84C4","#B565A7",
  "#88B04B","#7B6888","#C3447A","#009B77","#EFC050",
  "#7FCDCD","#DD4124","#5B5EA6","#E07A5F","#4BACC6",
  "#E8A0BF","#ea8a33","#9B2335","#C17BAE","#DECF3F",
  "#789262","#BC243C"
)
atac_colors <- setNames(hex_colors, atac_new_order)

known_suffixes <- "KOvCtrl_all|KOvCtrl_female|KOvCtrl_male|KO_FvM|Ctrl_FvM"

name_map <- c(
  "Astrocytes"                            = "Astrocytes",
  "Deep_layer_extratelencephalic_neurons" = "Deep-layer extratelencephalic neurons",
  "Endothelial_cells"                     = "Endothelial cells",
  "Layer_2_3_IT_neurons"                  = "Layer 2/3 IT neurons",
  "Layer_4_sensory_neurons"               = "Layer 4 sensory neurons",
  "Layer_5_6_IT_neurons"                  = "Layer 5/6 IT neurons",
  "Layer_5a_IT_neurons"                   = "Layer 5a IT neurons",
  "Layer_5b_PT_neurons"                   = "Layer 5b PT neurons",
  "Layer_6a_corticothalamic_neurons"      = "Layer 6a corticothalamic neurons",
  "Layer_6b_neurons"                      = "Layer 6b neurons",
  "Meningeal_fibroblasts"                 = "Meningeal fibroblasts",
  "Microglia"                             = "Microglia",
  "Oligodendrocyte_precursor_cells"       = "Oligodendrocyte precursor cells",
  "Oligodendrocytes"                      = "Oligodendrocytes",
  "PV__interneurons"                      = "PV+ interneurons",
  "SST__interneurons"                     = "SST+ interneurons",
  "VIP__interneurons"                     = "VIP+ interneurons"
)

# Helper: save plot as both PDF (quartz, macOS native) and PNG (ragg)
save_plot <- function(plot, stem, width = 10, height = 7) {
  # PDF — quartz uses CoreText; all system fonts including Unicode work natively
  quartz(type = "pdf", file = paste0(stem, ".pdf"), width = width, height = height)
  print(plot)
  dev.off()
  message("Saved: ", stem, ".pdf")
}

# ══════════════════════════════════════════════════════════════
# PLOT 1: Grouped bar chart
# ══════════════════════════════════════════════════════════════

bar_df <- summary_df |>
  filter(!grepl("significance_summary", contrast)) |>
  mutate(
    contrast_type = str_extract(contrast, known_suffixes),
    cell_prefix   = str_remove(contrast, paste0("__(", known_suffixes, ")$"))
  ) |>
  filter(contrast_type %in% c("KOvCtrl_all", "KOvCtrl_female",
                              "KOvCtrl_male", "KO_FvM")) |>
  mutate(cell_type = name_map[cell_prefix]) |>
  filter(!is.na(cell_type)) |>
  select(cell_type, contrast_type, sig_fdr05)

active_types <- bar_df |>
  group_by(cell_type) |>
  summarise(total = sum(sig_fdr05)) |>
  filter(total > 0) |>
  pull(cell_type)

bar_df <- bar_df |>
  filter(cell_type %in% active_types) |>
  mutate(
    contrast_type = factor(contrast_type,
                           levels = c("KOvCtrl_all", "KOvCtrl_female", "KOvCtrl_male", "KO_FvM"),
                           labels = c("KO vs Ctrl (all)", "KO vs Ctrl (\u2640)",
                                      "KO vs Ctrl (\u2642)", "KO \u2640 vs KO \u2642")
    ),
    # y-axis order = ATAC UMAP order (rev so top of legend = top of y-axis)
    cell_type = factor(cell_type,
                       levels = rev(intersect(atac_new_order, active_types)))
  )

bar_colors <- c(
  "KO vs Ctrl (all)"      = "#4E84C4",
  "KO vs Ctrl (\u2640)"   = "#C3447A",
  "KO vs Ctrl (\u2642)"   = "#009B77",
  "KO \u2640 vs KO \u2642" = "#EFC050"
)

p_bar <- ggplot(bar_df, aes(x = sig_fdr05, y = cell_type, fill = contrast_type)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7, color = NA) +
  scale_fill_manual(values = bar_colors, name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05)), labels = comma) +
  labs(
    title    = "Sex-dependent chromatin differences in KO exceed genotype effect",
    subtitle = "Significant DA peaks (FDR < 0.05) per contrast per cell type",
    x        = "Significant peaks (FDR < 0.05)",
    y        = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title         = element_text(face = "bold", size = 13),
    plot.subtitle      = element_text(size = 11, color = "grey40",
                                      margin = margin(b = 10)),
    axis.text.y        = element_text(size = 11),
    axis.text.x        = element_text(size = 10),
    legend.position    = "bottom",
    legend.text        = element_text(size = 11,
                                      family = "Arial Unicode MS"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    plot.margin        = margin(10, 15, 10, 10)
  ) +
  guides(fill = guide_legend(nrow = 1))

save_plot(p_bar, "DA_sex_vs_genotype_effect", width = 10, height = 7)

# ══════════════════════════════════════════════════════════════
# PLOT 2: Scatter — sex-driven vs KO-driven
# ══════════════════════════════════════════════════════════════

library(ggtext)   # install.packages("ggtext") if missing

scatter_df <- summary_df |>
  filter(!grepl("significance_summary", contrast)) |>
  mutate(
    contrast_type = str_extract(contrast, known_suffixes),
    cell_prefix   = str_remove(contrast, paste0("__(", known_suffixes, ")$"))
  ) |>
  filter(contrast_type %in% c("KOvCtrl_all", "KO_FvM")) |>
  mutate(cell_type = name_map[cell_prefix]) |>
  filter(!is.na(cell_type)) |>
  select(cell_type, contrast_type, sig_fdr05) |>
  pivot_wider(names_from = contrast_type, values_from = sig_fdr05,
              values_fill = 0) |>
  filter((KOvCtrl_all + KO_FvM) > 0) |>
  mutate(cell_type = factor(cell_type, levels = atac_new_order))

max_val <- max(scatter_df$KOvCtrl_all, scatter_df$KO_FvM)

p_scatter <- ggplot(scatter_df,
                    aes(x = KOvCtrl_all + 1, y = KO_FvM + 1,
                        color = cell_type)) +
  
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", color = "grey55", linewidth = 0.6) +
  
  annotate("richtext",
           x = 2, y = 2400,
           label = "*Sex-driven*<br>(KO &#9792; &ne; KO &#9794;)",
           color = "grey40", size = 3.8, hjust = 0.5,
           fill = NA, label.color = NA) +
  annotate("richtext",
           x = 2400, y = 2,
           label = "*KO-driven*<br>(KO &ne; Ctrl)",
           color = "grey40", size = 3.8, hjust = 0.5,
           fill = NA, label.color = NA) +
  
  geom_point(size = 4.5, alpha = 0.92) +
  
  scale_color_manual(
    values = atac_colors,
    breaks = atac_new_order,
    drop   = TRUE,
    name   = "Cell type"
  ) +
  
  scale_x_continuous(
    trans  = "log1p",
    limits = c(0, max_val + 1),
    breaks = c(0, 1, 10, 100, 1000, max_val),
    labels = c("0", "1", "10", "100", "1000", as.character(max_val)),
    expand = expansion(mult = c(0.05, 0.08))
  ) +
  scale_y_continuous(
    trans  = "log1p",
    limits = c(0, max_val + 1),
    breaks = c(0, 1, 10, 100, 1000, max_val),
    labels = c("0", "1", "10", "100", "1000", as.character(max_val)),
    expand = expansion(mult = c(0.05, 0.08))
  ) +
  
  labs(
    title    = "Sex-dependent chromatin remodeling exceeds KO genotype effect",
    subtitle = "Points above the diagonal: more DA peaks between KO &#9792; and KO &#9794; than between KO and Ctrl",
    x        = "KO vs Ctrl (all) \u2014 significant peaks (FDR < 0.05)",
    y        = "KO &#9792; vs KO &#9794; \u2014 significant peaks (FDR < 0.05)"
  ) +
  
  theme_minimal(base_size = 13) +
  theme(
    aspect.ratio     = 1,
    plot.title       = element_text(face = "bold", size = 13),
    plot.subtitle    = element_markdown(size = 10.5, color = "grey40",
                                        margin = margin(b = 8)),
    axis.title.x     = element_text(color = "black", family = "Arial Unicode MS"),  # <- updated
    axis.title.y     = element_markdown(color = "black"),                            # <- updated
    axis.text        = element_text(color = "black"),                                # <- add
    axis.line        = element_line(color = "black", linewidth = 0.5),
    axis.ticks       = element_line(color = "black", linewidth = 0.4),
    axis.ticks.length = unit(0.2, "cm"),
    legend.position  = "right",
    legend.text      = element_text(size = 9.5),
    legend.key.size  = unit(0.4, "cm"),
    panel.grid.major = element_line(color = "grey88", linetype = "dashed", linewidth = 0.2),
    panel.grid.minor = element_blank(),
    panel.border     = element_blank(),
    plot.margin      = margin(10, 10, 10, 10)
  )+
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 3.5)))

save_plot(p_scatter, "DA_scatter_sex_vs_ko", width = 10, height = 8)

###### MOTIF ######
# ══════════════════════════════════════════════════════════════════════════════
# STEP 5: MOTIF ENRICHMENT — full script (v3)
# Changes:
#   • FindMotifs run separately for gained (log2FC > 0.1) AND
#     lost (log2FC < -0.1) peaks per contrast
#   • Mirrored lollipop combining both directions:
#       right (teal)   = gained accessibility
#       left  (orange) = lost accessibility
#       dot alpha      = significance tier (FDR<0.05 / FDR<0.2 / nominal)
#       dot size       = fold enrichment
#   • No more shape encoding
# ══════════════════════════════════════════════════════════════════════════════

library(Signac)
library(Seurat)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Mmusculus.UCSC.mm10)
library(GenomicRanges)
library(ggplot2)
library(ggtext)
library(dplyr)
library(tidyr)
library(stringr)
library(patchwork)

# ── 0. Genome setup & peak filtering ─────────────────────────────────────────
main_chroms <- c(paste0("chr", 1:19), "chrX", "chrY", "chrM")
chrom_sizes <- seqlengths(BSgenome.Mmusculus.UCSC.mm10)[main_chroms]

peaks_gr <- StringToGRanges(rownames(obj))
peaks_gr <- peaks_gr[as.character(seqnames(peaks_gr)) %in% main_chroms]
seqlevels(peaks_gr)  <- main_chroms
seqlengths(peaks_gr) <- chrom_sizes

in_bounds   <- end(peaks_gr) <= seqlengths(peaks_gr)[as.character(seqnames(peaks_gr))]
peaks_gr    <- peaks_gr[in_bounds]
valid_peaks <- GRangesToString(peaks_gr)

message("Peaks removed (out-of-bounds or alt chroms): ", nrow(obj) - length(valid_peaks))
message("Peaks kept: ", length(valid_peaks))

DefaultAssay(obj) <- "peaks"
obj <- obj[rownames(obj) %in% valid_peaks, ]

# ── 1. Sex-effect baseline from Ctrl_FvM (FDR < 0.1) ─────────────────────────
ctrl_fvm_files     <- list.files("DA_results", pattern = "Ctrl_FvM\\.csv$", full.names = TRUE)
sex_baseline_peaks <- character(0)
for (f in ctrl_fvm_files) {
  da  <- read.csv(f, row.names = 1)
  sig <- rownames(da)[abs(da$avg_log2FC) > 0.1 & da$p_val_fdr < 0.1]
  sex_baseline_peaks <- union(sex_baseline_peaks, sig)
}
sex_baseline_peaks <- intersect(sex_baseline_peaks, valid_peaks)
message("Sex baseline peaks (union Ctrl_FvM FDR<0.1): ", length(sex_baseline_peaks))

# ── 2. Expressed TF filter ────────────────────────────────────────────────────
if ("RNA" %in% Assays(obj)) {
  rna_counts      <- GetAssayData(obj, assay = "RNA", layer = "counts")
  expressed_genes <- rownames(rna_counts)[rowSums(rna_counts > 0) >= 10]
  message("Expressed genes: ", length(expressed_genes))
} else {
  expressed_genes <- NULL
  message("No RNA assay — expression filter skipped")
}

# ── 3. JASPAR → filter → AddMotifs ───────────────────────────────────────────
pwm_set <- getMatrixSet(JASPAR2020,
                        opts = list(species = 10090, all_versions = FALSE))

if (!is.null(expressed_genes)) {
  tf_names_all <- sapply(pwm_set@listData, function(x) x@name)
  keep_motif   <- sapply(tf_names_all, function(tf) {
    parts <- sub("\\(.*\\)", "", unlist(strsplit(tf, "::"))) |> trimws()
    any(tolower(parts) %in% tolower(expressed_genes))
  })
  pwm_filtered <- pwm_set[keep_motif]
  message("Motifs retained after expression filter: ",
          length(pwm_filtered), " of ", length(pwm_set))
} else {
  pwm_filtered <- pwm_set
}

obj <- AddMotifs(obj, genome = BSgenome.Mmusculus.UCSC.mm10, pfm = pwm_filtered)

# ── 4. GC/length stats for matched background ─────────────────────────────────
message("Computing region stats (GC content)...")
obj <- RegionStats(obj, genome = BSgenome.Mmusculus.UCSC.mm10, assay = "peaks")

# ── 5. Helper: run FindMotifs for one direction ───────────────────────────────
run_motifs <- function(obj, da_df, direction = c("gained", "lost"),
                       valid_peaks, meta_feat, all_open,
                       sex_baseline = NULL) {
  direction <- match.arg(direction)
  fc_filter <- if (direction == "gained") da_df$avg_log2FC > 0.1 else da_df$avg_log2FC < -0.1
  peaks     <- rownames(da_df)[fc_filter & da_df$p_val_fdr < 0.1]
  peaks     <- intersect(peaks, valid_peaks)
  
  if (!is.null(sex_baseline)) {
    n_before <- length(peaks)
    peaks    <- setdiff(peaks, sex_baseline)
    message("    Sex-baseline excluded: ", n_before - length(peaks),
            " → ", length(peaks), " ", direction, " peaks remain")
  }
  
  if (length(peaks) < 10) {
    message("    Too few peaks (", length(peaks), ") — skipping ", direction)
    return(NULL)
  }
  
  peaks_in_meta <- intersect(peaks, rownames(meta_feat))
  bg <- tryCatch(
    MatchRegionStats(
      meta.feature  = meta_feat[all_open, , drop = FALSE],
      query.feature = meta_feat[peaks_in_meta, , drop = FALSE],
      n             = min(50000, length(all_open))
    ),
    error = function(e) { message("    MatchRegionStats failed: ", e$message); NULL }
  )
  
  res <- tryCatch(
    if (is.null(bg)) FindMotifs(obj, features = peaks)
    else             FindMotifs(obj, features = peaks, background = bg),
    error = function(e) { message("    FindMotifs FAILED: ", e$message); NULL }
  )
  
  if (!is.null(res)) {
    res$direction  <- direction
    res$n_peaks    <- length(peaks)
  }
  res
}

# ── 6. Enrichment loop ────────────────────────────────────────────────────────
dir.create("motif_results", showWarnings = FALSE)
dir.create("motif_plots",   showWarnings = FALSE)

da_files <- list.files("DA_results", pattern = "\\.csv$", full.names = TRUE)
if (length(da_files) == 0) stop("No DA result files found. Run Part 4 first.")
message("Found ", length(da_files), " DA result files...")

all_open  <- rownames(obj)
meta_feat <- obj[["peaks"]]@meta.features

for (f in da_files) {
  nm <- gsub("\\.csv$", "", basename(f))
  message("\nProcessing: ", nm)
  
  is_ko_ctrl <- grepl("KOvCtrl", nm)
  baseline   <- if (is_ko_ctrl) sex_baseline_peaks else NULL
  
  da_df <- read.csv(f, row.names = 1)
  
  for (dir in c("gained", "lost")) {
    outfile <- file.path("motif_results", paste0(nm, "_", dir, "_motifs.csv"))
    if (file.exists(outfile)) { message("  Skipping (done): ", dir); next }
    
    res <- run_motifs(obj, da_df, direction = dir,
                      valid_peaks = valid_peaks,
                      meta_feat   = meta_feat,
                      all_open    = all_open,
                      sex_baseline = baseline)
    if (!is.null(res)) {
      write.csv(res, outfile)
      message("  Saved [", dir, "]: ",
              sum(res$pvalue < 0.05), " nominal-sig | ",
              sum(res$p.adjust < 0.05), " FDR-sig motifs")
    }
  }
}
message("\n>>> Motif enrichment complete.")

# ══════════════════════════════════════════════════════════════════════════════
# PLOTTING
# ══════════════════════════════════════════════════════════════════════════════

clean_tf <- function(x) sub("\\(.*\\)", "", x) |> trimws()

# Palette
COL_GAINED <- "#01696f"   # teal  — gained accessibility
COL_LOST   <- "#bb653b"   # orange-brown — lost accessibility

motif_theme <- theme_minimal(base_size = 11) +
  theme(
    plot.title         = element_text(face = "bold", size = 11),
    plot.subtitle      = element_text(size = 9, color = "grey40",
                                      margin = margin(b = 6)),
    axis.text          = element_text(color = "black"),
    axis.title         = element_text(color = "black"),
    axis.line.x        = element_line(color = "black", linewidth = 0.4),
    axis.ticks.x       = element_line(color = "black", linewidth = 0.3),
    axis.ticks.length  = unit(0.15, "cm"),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(color = "grey88", linetype = "dashed",
                                      linewidth = 0.3),
    panel.grid.minor   = element_blank(),
    legend.text        = element_text(size = 9),
    legend.title       = element_text(size = 9, face = "bold")
  )

# ── Mirrored lollipop (gained right / lost left) ──────────────────────────────
# sig_col encodes p-value type shown: "rawp" or "padj"
plot_mirrored_lollipop <- function(gained_df, lost_df,
                                   title, top_n = 15,
                                   sig_col = "rawp") {
  
  p_col   <- if (sig_col == "padj") "p.adjust" else "pvalue"
  p_label <- if (sig_col == "padj") "-log10(p.adjust)" else "-log10(p-value)"
  p_sub   <- if (sig_col == "padj")
    "BH-adjusted — right = gained accessibility / left = lost accessibility"
  else
    "Nominal — right = gained accessibility / left = lost accessibility"
  
  # Color palette: 3 shades per direction (full → medium → light)
  color_map <- c(
    "Gained — FDR < 0.05"  = "#01696f",
    "Gained — FDR < 0.20"  = "#5fa8ad",
    "Gained — Nominal"     = "#b0d8da",
    "Lost — FDR < 0.05"    = "#bb653b",
    "Lost — FDR < 0.20"    = "#d4956d",
    "Lost — Nominal"       = "#ecc9b2"
  )
  
  prep <- function(df, dir) {
    if (is.null(df)) return(NULL)
    label_prefix <- if (dir == "gained") "Gained" else "Lost"
    df |>
      filter(pvalue < 0.05) |>
      arrange(.data[[p_col]]) |>
      mutate(
        tf_label    = clean_tf(motif.name),
        log_p       = -log10(.data[[p_col]]),
        fold_enrich = percent.observed / percent.background,
        sig_color   = factor(
          case_when(
            p.adjust < 0.05 ~ paste0(label_prefix, " — FDR < 0.05"),
            p.adjust < 0.20 ~ paste0(label_prefix, " — FDR < 0.20"),
            TRUE            ~ paste0(label_prefix, " — Nominal")
          ),
          levels = names(color_map)
        ),
        direction = dir
      ) |>
      # keep only the most significant motif per TF name
      group_by(tf_label) |>
      slice_min(order_by = .data[[p_col]], n = 1, with_ties = FALSE) |>
      ungroup() |>
      head(top_n)
  }
  
  g_df <- prep(gained_df, "gained")
  l_df <- prep(lost_df,   "lost")
  if (is.null(g_df) && is.null(l_df)) return(NULL)
  
  # Order TFs by max log_p across both directions
  all_tf <- bind_rows(g_df, l_df) |>
    group_by(tf_label) |>
    summarise(max_logp = max(log_p), .groups = "drop") |>
    arrange(max_logp) |>
    pull(tf_label)
  
  plot_df <- bind_rows(g_df, l_df) |>
    mutate(
      tf_label = factor(tf_label, levels = all_tf),
      x_val    = ifelse(direction == "gained", log_p, -log_p)
    )
  
  sig_line <- -log10(0.05)
  
  # x-axis upper limit with a bit of padding
  x_max <- max(abs(plot_df$x_val), na.rm = TRUE) * 1.08
  
  ggplot(plot_df, aes(x = x_val, y = tf_label)) +
    
    # zero spine
    geom_vline(xintercept = 0,
               color = "grey30", linewidth = 0.5) +
    
    # Stems — colored to match dot, fully opaque
    geom_segment(aes(xend = 0, yend = tf_label, color = sig_color),
                 linewidth = 0.7) +
    
    # Dots
    geom_point(aes(size = fold_enrich, color = sig_color)) +
    
    # p = 0.05 reference lines — appear in legend via linetype
    geom_vline(aes(xintercept = sig_line,  linetype = "p = 0.05"),
               color = "darkred", linewidth = 0.3,
               key_glyph = "path") +
    geom_vline(aes(xintercept = -sig_line, linetype = "p = 0.05"),
               color = "darkred", linewidth = 0.3,
               key_glyph = "path") +
    # zero spine (no legend entry needed)
    geom_vline(xintercept = 0, color = "grey30", linewidth = 0.5) +
    
    scale_color_manual(
      values = color_map,
      drop   = TRUE,
      name   = NULL,
      guide  = guide_legend(
        override.aes = list(size = 3.5, shape = 16),
        ncol = 1
      )
    ) +
    scale_linetype_manual(
      values = c("p = 0.05" = "dashed"),
      name   = NULL
    ) +
    scale_size_continuous(range = c(2, 7), name = "Fold\nenrichment") +
    scale_x_continuous(
      labels = function(x) round(abs(x), 1),
      limits = c(-x_max, x_max),
      name   = p_label,
      expand = expansion(mult = c(0.02, 0.02))
    ) +
    # clip = "off" lets annotations render outside the panel boundary
    labs(title = title, subtitle = p_sub, y = NULL) +
    motif_theme +
    theme(
      legend.position  = "right",
      legend.key.size  = unit(0.4, "cm"),
      legend.spacing.y = unit(0.15, "cm"),
      plot.margin      = margin(t = 10, r = 10, b = 10, l = 10)
    )
}

# ── Per-contrast heatmap (unchanged from before, now uses _gained_ files) ─────
plot_motif_heatmap <- function(contrast_pattern, direction = "gained",
                               top_n_per_ct = 5, title = contrast_pattern,
                               p_cap = 1e-10, use_padj = FALSE) {
  
  pattern    <- paste0(contrast_pattern, "_", direction, "_motifs\\.csv$")
  motif_files <- list.files("motif_results", pattern = pattern, full.names = TRUE)
  if (length(motif_files) == 0) {
    message("  No files for: ", contrast_pattern, " [", direction, "]")
    return(NULL)
  }
  
  p_col <- if (use_padj) "p.adjust" else "pvalue"
  p_cut <- if (use_padj) 0.20       else 0.05
  p_lbl <- if (use_padj) "-log10(p.adjust)" else "-log10(p-value)"
  p_sub <- paste0(if (use_padj) "BH-adjusted" else "Nominal",
                  " | ", direction, " accessibility")
  fill_col <- if (direction == "gained") "#01696f" else "#bb653b"
  
  all_data <- lapply(motif_files, function(f) {
    ct  <- gsub(paste0("__", contrast_pattern, "_", direction,
                       "_motifs\\.csv$"), "", basename(f))
    df  <- read.csv(f, row.names = 1)
    df$cell_type <- gsub("_", " ", ct)
    df
  }) |> bind_rows()
  
  top_motifs <- all_data |>
    filter(.data[[p_col]] < p_cut) |>
    group_by(cell_type) |>
    slice_min(.data[[p_col]], n = top_n_per_ct, with_ties = FALSE) |>
    ungroup() |>
    pull(motif.name) |>
    unique()
  
  if (length(top_motifs) == 0) {
    message("  No sig motifs for: ", contrast_pattern, " [", direction, "]")
    return(NULL)
  }
  
  heat_long <- all_data |>
    filter(motif.name %in% top_motifs) |>
    mutate(
      log_p    = -log10(pmax(.data[[p_col]], p_cap)),
      tf_label = clean_tf(motif.name)
    ) |>
    select(cell_type, tf_label, log_p) |>
    complete(cell_type, tf_label, fill = list(log_p = 0))
  
  ggplot(heat_long, aes(x = cell_type, y = tf_label, fill = log_p)) +
    geom_tile(color = "white", linewidth = 0.35) +
    scale_fill_gradient(low = "white", high = fill_col, name = p_lbl) +
    labs(title = title, subtitle = p_sub, x = NULL, y = NULL) +
    theme_minimal(base_size = 10) +
    theme(
      plot.title    = element_text(face = "bold", size = 11),
      plot.subtitle = element_text(size = 9, color = "grey40"),
      axis.text.x   = element_text(color = "black", angle = 35,
                                   hjust = 1, size = 9),
      axis.text.y   = element_text(color = "black", size = 8),
      panel.grid    = element_blank()
    )
}

# ── 7. Generate plots — 2 per result file (rawp + padj) ───────────────────────
# Get unique contrast names (without _gained / _lost suffix)
all_motif_files <- list.files("motif_results", pattern = "_motifs\\.csv$",
                              full.names = TRUE)
contrast_nms <- all_motif_files |>
  basename() |>
  gsub("_(gained|lost)_motifs\\.csv$", "", x = _) |>
  unique()

message("Generating mirrored lollipop plots for ",
        length(contrast_nms), " contrasts...")

for (nm in contrast_nms) {
  gained_file <- file.path("motif_results", paste0(nm, "_gained_motifs.csv"))
  lost_file   <- file.path("motif_results", paste0(nm, "_lost_motifs.csv"))
  
  g_df <- if (file.exists(gained_file)) read.csv(gained_file, row.names = 1) else NULL
  l_df <- if (file.exists(lost_file))   read.csv(lost_file,   row.names = 1) else NULL
  
  if (is.null(g_df) && is.null(l_df)) {
    message("  No data for: ", nm); next
  }
  
  label <- gsub("_", " ", gsub("__", "  |  ", nm))
  
  for (sc in c("rawp", "padj")) {
    p <- plot_mirrored_lollipop(g_df, l_df, title = label,
                                top_n = 15, sig_col = sc)
    if (!is.null(p)) {
      ggsave(file.path("motif_plots",
                       paste0(nm, "_lollipop_", sc, ".pdf")),
             p, width = 9, height = 6)
      message("  Saved [", sc, "]: ", nm)
    } else {
      message("  No sig motifs for plot [", sc, "]: ", nm)
    }
  }
}

# ── 8. Cross-cell-type heatmaps ───────────────────────────────────────────────
contrast_labels <- list(
  KOvCtrl_all    = "KO vs Ctrl (all)",
  KO_FvM         = "KO female vs male",
  Ctrl_FvM       = "Ctrl female vs male (sex baseline)",
  KOvCtrl_female = "KO vs Ctrl (female)",
  KOvCtrl_male   = "KO vs Ctrl (male)"
)

for (contrast in names(contrast_labels)) {
  lbl <- contrast_labels[[contrast]]
  for (dir in c("gained", "lost")) {
    for (adj in c(FALSE, TRUE)) {
      suffix <- paste0(dir, "_", if (adj) "padj" else "rawp")
      p_heat <- plot_motif_heatmap(
        contrast_pattern = contrast,
        direction        = dir,
        top_n_per_ct     = 5,
        title            = paste0("Top enriched motifs — ", lbl),
        use_padj         = adj
      )
      if (!is.null(p_heat)) {
        ggsave(file.path("motif_plots",
                         paste0("heatmap_", suffix, "_", contrast, ".pdf")),
               p_heat, width = 11, height = 8)
        message("Heatmap [", suffix, "]: ", contrast)
      }
    }
  }
}

message("\n>>> All done.")
message("CSVs  → motif_results/  (*_gained_motifs.csv + *_lost_motifs.csv)")
message("Plots → motif_plots/")
message("  Per contrast: <name>_lollipop_rawp.pdf + <name>_lollipop_padj.pdf")
message("  Heatmaps:     heatmap_{gained|lost}_{rawp|padj}_<contrast>.pdf")



####### SAVE ######
# ── STEP 6: SAVE OBJECT ───────────────────────────────────────

saveRDS(obj, "scATAC_04_final.rds")

message("\n=== DONE ===")
message("Output folders:")
message("  DA_results/    ->  5 contrast CSVs per cell type")
message("  motif_results/ ->  motif enrichment CSVs per contrast")