#install.packages('devtools')
#devtools::install_github('immunogenomics/presto')

# ══════════════════════════════════════════════════════════════════════════════
# PART 5: INTEGRATIVE snATAC-seq + snRNA-seq ANALYSIS
# Males only — KO vs Ctrl
#
# Requires in session (or reload):
#   obj          — full ATAC object (scATAC_04_final.rds)
#   rna          — RNA object (rna_for_atac.rds) [male-only by design]
#   atac_colors  — named color vector from Part 3
#   name_map     — cell type key map from Part 3
#   pwm_filtered — expressed JASPAR motifs from Part 5 motif section
#
# Flow:
#   Step 1: Subset male ATAC cells
#   Step 2: Co-embedding (imputed RNA from ATAC + RNA on shared UMAP)
#   Step 3: Peak–gene linking (LinkPeaks)
#   Step 4: DA peaks × DEGs overlap + concordance plot
#   Step 5: chromVAR TF activity scores (male, per cell type)
#   Step 6: Double-hit TF analysis (chromVAR activity + RNA expression)
#   Step 7: Coverage plots for top concordant loci
#   Step 8: Per-cell-type integration summary table
# ══════════════════════════════════════════════════════════════════════════════

# ── Additional packages ───────────────────────────────────────────────────────
#BiocManager::install("betterChromVAR")
library(betterChromVAR)
library(Signac)
library(Seurat)
library(GenomicRanges)
library(EnsDb.Mmusculus.v79)
library(BSgenome.Mmusculus.UCSC.mm10)
library(JASPAR2020)
library(TFBSTools)
library(SummarizedExperiment)
library(ggplot2)
library(ggrepel)
library(ggtext)
library(dplyr)
library(tidyr)
library(patchwork)
library(scales)
library(motifmatchr)

setwd("~/Downloads/Seurat_scATAC-seq")
set.seed(42)

# ── Reload objects if starting a new session ─────────────────────────────────
obj <- readRDS("scATAC_04_final.rds")
rna <- readRDS("rna_for_atac.rds")
rna$cell_type <- as.character(Idents(rna))
rna$sex <- "male"
rna$genotype <- ifelse(
  rna$Genotype == "PV-Cre/tdTom",
  "Ctrl",
  "KO"
)

# ── Redefine name_map if not in session (copy from Part 3) ───────────────────
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

# ── Helper ────────────────────────────────────────────────────────────────────
clean_tf <- function(x) sub("\\(.*\\)", "", x) |> trimws()

save_plot <- function(plot, stem, width = 12, height = 8) {
  pdf(paste0(stem, ".pdf"), width = width, height = height)
  print(plot)
  dev.off()
  message("Saved: ", stem, ".pdf")
}

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

known_suffixes <- "KOvCtrl_all|KOvCtrl_female|KOvCtrl_male|KO_FvM|Ctrl_FvM"

# atac_colors is rna_colors_atac — just an alias used in Part 3
atac_colors <- rna_colors_atac

# ══════════════════════════════════════════════════════════════════════════════
# PATCH: correct clean_tf regex — header has escaped backslashes from copy-paste
# ══════════════════════════════════════════════════════════════════════════════
clean_tf <- function(x) sub("\\(.*\\)", "", x) |> trimws()

# ══════════════════════════════════════════════════════════════════════════════
# STEP 1: SUBSET MALE ATAC DATA
# ══════════════════════════════════════════════════════════════════════════════

obj_male <- subset(obj, subset = sex == "male")

message("Male ATAC cells: ", ncol(obj_male))
message("  KO:   ", sum(obj_male$genotype == "KO"))
message("  Ctrl: ", sum(obj_male$genotype == "Ctrl"))

stopifnot(
  "cell_type" %in% colnames(obj_male@meta.data),
  "genotype"  %in% colnames(obj_male@meta.data),
  ncol(obj_male) > 100
)

message("\nATAC cells per cell type × genotype (male):")
print(table(obj_male$cell_type, obj_male$genotype))

# ══════════════════════════════════════════════════════════════════════════════
# STEP 2: SIDE-BY-SIDE UMAP — label transfer already done, no re-anchoring
# ══════════════════════════════════════════════════════════════════════════════

# atac_colors and rna_colors_original must be in session from Part 3
# Factor levels — drop any cell types not present in subset
rna_types  <- intersect(names(rna_colors_original), unique(as.character(rna$cell_type)))
atac_types <- intersect(names(atac_colors), as.character(levels(obj_male$cell_type)))

rna$cell_type       <- factor(rna$cell_type,       levels = rna_types)
obj_male$cell_type  <- factor(obj_male$cell_type,  levels = atac_types)
Idents(obj_male)    <- "cell_type"

p_rna_umap <- DimPlot(rna, group.by = "cell_type",
                      cols = rna_colors_original[rna_types],
                      pt.size = 0.1, label = FALSE) +
  ggtitle("snRNA-seq (male)") +
  theme(legend.text = element_text(size = 8),
        legend.key.size = unit(0.35, "cm")) +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 3)))

p_atac_umap <- DimPlot(obj_male, group.by = "cell_type",
                       cols = atac_colors[atac_types],
                       pt.size = 0.1, label = FALSE) +
  ggtitle("snATAC-seq (male)") +
  theme(legend.text = element_text(size = 8),
        legend.key.size = unit(0.35, "cm")) +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 3)))

p_geno <- DimPlot(obj_male, group.by = "genotype", pt.size = 0.1,
                  cols = c("Ctrl" = "#009B77", "KO" = "#C3447A")) +
  ggtitle("snATAC-seq — genotype (male)")

save_plot(
  (p_rna_umap | p_atac_umap) / p_geno + plot_layout(heights = c(1.5, 1)),
  "integration_umap_overview", width = 18, height = 14
)

# ══════════════════════════════════════════════════════════════════════════════
# STEP 3: SLICE GENE ACTIVITY ASSAY + LINKPEAKS
# Slicing from full obj is equivalent to recomputing but takes seconds
# ══════════════════════════════════════════════════════════════════════════════

# 3a. Slice gene activity counts for male cells only
# Run GeneActivity on the male subset directly
DefaultAssay(obj_male) <- "peaks"

gene_activities_male <- GeneActivity(obj_male)

obj_male[["RNA"]] <- CreateAssayObject(counts = gene_activities_male)
obj_male <- NormalizeData(
  obj_male,
  assay                = "RNA",
  normalization.method = "LogNormalize",
  scale.factor         = median(obj_male$nCount_peaks)
)
message("Gene activity added: ", nrow(obj_male[["RNA"]]),
        " genes \u00d7 ", ncol(obj_male), " cells")

saveRDS(obj_male, "scATAC_05_male_geneactivity.rds")

Assays(obj_male)                          # should show: peaks  RNA
dim(obj_male[["RNA"]])                    # should show: genes × 13036 cells

# 3b. Transfer region GC/length stats from full obj (avoid rerunning RegionStats)
meta_full <- obj[["peaks"]]@meta.features
if (!is.null(meta_full) && "GC.percent" %in% colnames(meta_full)) {
  shared_peaks <- intersect(rownames(obj_male), rownames(meta_full))
  if (length(shared_peaks) == nrow(obj_male)) {
    obj_male[["peaks"]]@meta.features <- meta_full[rownames(obj_male), , drop = FALSE]
    message("Region stats transferred from full obj")
  } else {
    message("Peak mismatch detected — recomputing RegionStats on obj_male")
    obj_male <- RegionStats(obj_male, genome = BSgenome.Mmusculus.UCSC.mm10)
  }
} else {
  message("No region stats in obj — computing on obj_male")
  obj_male <- RegionStats(obj_male, genome = BSgenome.Mmusculus.UCSC.mm10)
}

# 3c. Variable genes for LinkPeaks (intersection of gene activity and RNA)
rna_tmp    <- FindVariableFeatures(NormalizeData(rna), nfeatures = 5000)
genes_use  <- intersect(rownames(obj_male[["RNA"]]), VariableFeatures(rna_tmp))
message("Genes for LinkPeaks: ", length(genes_use))

if (length(genes_use) < 200) {
  warning(length(genes_use), " shared genes — very low. ",
          "Consider: genes_use <- rownames(obj_male[['RNA']])")
}

# 3d. Run LinkPeaks
DefaultAssay(obj_male) <- "peaks"
obj_male <- LinkPeaks(
  object           = obj_male,
  peak.assay       = "peaks",
  expression.assay = "RNA",
  genes.use        = genes_use,
  distance         = 5e5,
  min.cells        = 10,
  score_cutoff     = 0.05
)

links_df  <- as.data.frame(Links(obj_male))
links_sig <- links_df[!is.na(links_df$score) & links_df$score < 0.05, ]
write.csv(links_df, "peak_gene_links.csv", row.names = FALSE)
message("Links stored: ", nrow(links_df), " | Significant (p<0.05): ", nrow(links_sig))

# Relax threshold if no significant links
if (nrow(links_sig) == 0) {
  warning("No links at p<0.05 — relaxing to p<0.1 for overlap analysis")
  links_sig <- links_df[!is.na(links_df$score) & links_df$score < 0.1, ]
  message("Links at p<0.1: ", nrow(links_sig))
}

saveRDS(obj_male, "scATAC_05_male_linked.rds")

# ══════════════════════════════════════════════════════════════════════════════
# STEP 4: DEGs + DA×DEG OVERLAP + CONCORDANCE SCATTER
# ══════════════════════════════════════════════════════════════════════════════

# 4a. Checkpoint — compare RNA vs ATAC cell type names

# Collapse all oligodendrocyte subtypes → "Oligodendrocytes" to match ATAC
ol_subtypes <- c("Newly formed oligodendrocytes",
                 "Perineuronal oligodendrocytes",
                 "Myelinating oligodendrocytes")

rna$cell_type <- ifelse(
  rna$cell_type %in% ol_subtypes,
  "Oligodendrocytes",
  rna$cell_type
)

# Re-set Idents to the updated cell_type
rna_ct_names  <- sort(unique(as.character(rna$cell_type)))
atac_ct_names <- sort(levels(obj_male$cell_type))
shared_cts    <- intersect(rna_ct_names, atac_ct_names)

message("\n── RNA cell types (", length(rna_ct_names), ") ──")
print(rna_ct_names)
message("── ATAC cell types (", length(atac_ct_names), ") ──")
print(atac_ct_names)
message("── Shared (", length(shared_cts), ") ──")
print(shared_cts)



message("RNA cell types after OL collapse:")
print(sort(table(rna$cell_type), decreasing = TRUE))

#### DIAGNOSTICS
# ── Diagnostic 1: How many significant links do you have? ────────────────────
message("Total links: ",     nrow(links_df))
message("Sig links (p<0.05): ", nrow(links_sig))

# # If links_sig is empty, relax the threshold
# if (nrow(links_sig) == 0) {
#   message("No sig links at p<0.05 — relaxing to p<0.2")
#   links_sig <- links_df[!is.na(links_df$score) & links_df$score < 0.2, ]
#   message("Links at p<0.2: ", nrow(links_sig))
# }

# ── Diagnostic 2: Do name_map keys match your actual DA filenames? ────────────
da_files_male <- list.files("DA_results", pattern = "__KOvCtrl_male\\.csv$",
                            full.names = TRUE)
ct_keys_found <- gsub("__KOvCtrl_male$", "",
                      gsub("\\.csv$", "", basename(da_files_male)))
message("\nDA file keys:")
print(ct_keys_found)
message("\nname_map keys:")
print(names(name_map))
message("\nMatched:")
print(ct_keys_found[ct_keys_found %in% names(name_map)])
message("Unmatched (will be skipped):")
print(ct_keys_found[!ct_keys_found %in% names(name_map)])

# ── Diagnostic 3: Check DEG files exist for matched types ────────────────────
for (k in ct_keys_found[ct_keys_found %in% names(name_map)]) {
  ct      <- name_map[k]
  ct_safe <- gsub("[^A-Za-z0-9]", "_", ct)
  f       <- file.path("DEG_results", paste0(ct_safe, "__KOvCtrl_male_DEG.csv"))
  message(ct, " → DEG file exists: ", file.exists(f))
}

# ── Diagnostic 4: For one cell type, check how many sig DA peaks and DEGs ────
# Change "PV__interneurons" to whatever key is in your da_files
test_key    <- ct_keys_found[1]
test_ct     <- name_map[test_key]
test_ct_safe <- gsub("[^A-Za-z0-9]", "_", test_ct)

test_da  <- read.csv(file.path("DA_results",  paste0(test_key,     "__KOvCtrl_male.csv")),  row.names = 1)
test_deg <- read.csv(file.path("DEG_results", paste0(test_ct_safe, "__KOvCtrl_male_DEG.csv")), row.names = 1)

message("\nTest cell type: ", test_ct)
message("  Sig DA peaks (FDR<0.05): ", sum(test_da$p_val_fdr  < 0.05, na.rm = TRUE))
message("  Sig DEGs     (FDR<0.05): ", sum(test_deg$p_val_fdr < 0.05, na.rm = TRUE))
message("  Linked genes in links_sig: ",
        sum(links_sig$gene %in% rownames(test_deg)[test_deg$p_val_fdr < 0.05]))
message("  Linked peaks that are also DA: ",
        sum(intersect(
          links_sig$peak[links_sig$gene %in% rownames(test_deg)[test_deg$p_val_fdr < 0.05]],
          rownames(test_da)[test_da$p_val_fdr < 0.05]
        ) |> length()))

# WHAT TO LOOK FOR:
#   shared_cts > 0 required for overlap analysis
#   If RNA names differ from ATAC names, uncomment the mapping block below
#   and fill in the mismatched pairs before continuing.

# ── Manual name mapping if RNA and ATAC labels differ ─────────────────────────
# rna_to_atac <- c(
#   "PV interneurons"  = "PV+ interneurons",   # RNA name = ATAC name
#   "SST interneurons" = "SST+ interneurons"
# )
# rna$cell_type <- factor(
#   ifelse(as.character(rna$cell_type) %in% names(rna_to_atac),
#          rna_to_atac[as.character(rna$cell_type)],
#          as.character(rna$cell_type)),
#   levels = levels(rna$cell_type)
# )
# shared_cts <- intersect(unique(as.character(rna$cell_type)),
#                          levels(obj_male$cell_type))

stopifnot("Shared cell types = 0 — fix RNA/ATAC name mismatch above" = length(shared_cts) > 0)

# 4b. DEGs per cell type (RNA, KO vs Ctrl, male)
dir.create("DEG_results", showWarnings = FALSE)

for (ct in shared_cts) {
  ct_safe <- gsub("[^A-Za-z0-9]", "_", ct)
  outfile <- file.path("DEG_results", paste0(ct_safe, "__KOvCtrl_male_DEG.csv"))
  if (file.exists(outfile)) { message("Skipping DEG: ", ct); next }
  
  cells_ko   <- rownames(rna@meta.data[rna$cell_type == ct & rna$genotype == "KO",   ])
  cells_ctrl <- rownames(rna@meta.data[rna$cell_type == ct & rna$genotype == "Ctrl", ])
  
  if (length(cells_ko) < 10 | length(cells_ctrl) < 10) {
    message("  Skip DEG [", ct, "] — too few cells (",
            length(cells_ko), " KO / ", length(cells_ctrl), " Ctrl)")
    next
  }
  
  res <- tryCatch(
    FindMarkers(rna,
                ident.1         = cells_ko,
                ident.2         = cells_ctrl,
                test.use        = "wilcox",
                min.pct         = 0.1,
                logfc.threshold = 0.26),
    error = function(e) { message("  DEG error [", ct, "]: ", e$message); NULL }
  )
  if (!is.null(res)) {
    res$p_val_fdr <- p.adjust(res$p_val, method = "BH")
    write.csv(res, outfile)
    message("  DEG [", ct, "]: ",
            sum(res$p_val_fdr < 0.05, na.rm = TRUE), " sig genes (FDR<0.05) | ",
            sum(res$p_val     < 0.05, na.rm = TRUE), " nominal")
  }
}

# ══════════════════════════════════════════════════════════════════════════════
# 4c. DA × DEG overlap
# NOTE: DA side uses nominal p<0.05 + |log2FC|>0.1 because the severe
# Ctrl/KO cell imbalance (~1:4) makes FDR very conservative on the DA peaks.
# fdr_atac is still stored in the output for downstream filtering.
# DEG side keeps FDR<0.05 as it is well-powered (hundreds of cells per type).
# ══════════════════════════════════════════════════════════════════════════════

da_files_male <- list.files("DA_results", pattern = "__KOvCtrl_male\\.csv$",
                            full.names = TRUE)
message("\nBuilding DA-DEG overlap from ", length(da_files_male), " DA files...")

overlap_results <- list()

for (da_f in da_files_male) {
  nm     <- gsub("\\.csv$", "", basename(da_f))
  ct_key <- gsub("__KOvCtrl_male$", "", nm)
  ct     <- name_map[ct_key]
  if (is.na(ct) | !ct %in% shared_cts) {
    message("  Skip (unmapped or not shared): ", ct_key)
    next
  }
  
  ct_safe  <- gsub("[^A-Za-z0-9]", "_", ct)
  deg_file <- file.path("DEG_results", paste0(ct_safe, "__KOvCtrl_male_DEG.csv"))
  if (!file.exists(deg_file)) {
    message("  Skip (no DEG file): ", ct)
    next
  }
  
  da_df  <- read.csv(da_f,     row.names = 1)
  deg_df <- read.csv(deg_file, row.names = 1)
  
  # DA: nominal p<0.05 + directional effect filter
  da_sig  <- rownames(da_df)[!is.na(da_df$p_val) &
                               da_df$p_val < 0.05 &
                               abs(da_df$avg_log2FC) > 0.1]
  
  # DEG: FDR<0.05 (well-powered)
  deg_sig <- rownames(deg_df)[!is.na(deg_df$p_val) & deg_df$p_val < 0.05]
  
  message("  [", ct, "] DA nominal: ", length(da_sig),
          " | DEG FDR<0.05: ", length(deg_sig))
  
  if (length(da_sig) == 0 | length(deg_sig) == 0) {
    message("    -> no sig features, skipping")
    next
  }
  
  linked_peaks  <- links_sig$peak[links_sig$gene %in% deg_sig]
  overlap_peaks <- intersect(da_sig, linked_peaks)
  message("    -> linked peaks overlapping DA: ", length(overlap_peaks))
  
  if (length(overlap_peaks) == 0) next
  
  ov <- data.frame(peak = overlap_peaks) |>
    left_join(links_sig[links_sig$peak %in% overlap_peaks,
                        c("peak", "gene", "score")],
              by = "peak", relationship = "many-to-many") |>
    left_join(data.frame(gene       = rownames(deg_df),
                         log2FC_rna = deg_df$avg_log2FC,
                         fdr_rna    = deg_df$p_val_fdr),
              by = "gene") |>
    left_join(data.frame(peak        = rownames(da_df),
                         log2FC_atac = da_df$avg_log2FC,
                         p_val_atac  = da_df$p_val,
                         fdr_atac    = da_df$p_val_fdr),
              by = "peak") |>
    filter(!is.na(gene)) |>
    mutate(cell_type  = ct,
           concordant = sign(log2FC_rna) == sign(log2FC_atac))
  
  if (nrow(ov) > 0) {
    overlap_results[[ct]] <- ov
    message("    -> ", nrow(ov), " peak-gene pairs stored (",
            sum(ov$concordant), " concordant)")
  }
}

if (length(overlap_results) > 0) {
  overlap_df <- bind_rows(overlap_results)
  write.csv(overlap_df, "DA_DEG_overlap.csv", row.names = FALSE)
  message("\nOverlapping peak-gene pairs: ", nrow(overlap_df))
  message("Unique genes:                 ", length(unique(overlap_df$gene)))
  message("Concordant pairs:             ", sum(overlap_df$concordant, na.rm = TRUE))
  message("Cell types with overlaps:     ", length(unique(overlap_df$cell_type)))
} else {
  overlap_df <- data.frame()
  message("WARNING: Still no overlaps after threshold relaxation.")
  message("Run this to check links_sig content:")
  message("  head(links_sig)")
  message("  table(links_sig$gene %in% rownames(deg_df))")
}

# ══════════════════════════════════════════════════════════════════════════════
# 4d. Concordance scatter — only if overlaps found
# ══════════════════════════════════════════════════════════════════════════════

if (nrow(overlap_df) > 5) {
  
  conc_df <- overlap_df |>
    filter(!is.na(log2FC_rna), !is.na(log2FC_atac)) |>
    mutate(
      label_gene = ifelse(
        abs(log2FC_rna)  > quantile(abs(log2FC_rna),  0.90, na.rm = TRUE) &
          abs(log2FC_atac) > quantile(abs(log2FC_atac), 0.90, na.rm = TRUE),
        gene, NA),
      hit_class = factor(case_when(
        concordant  & fdr_rna < 0.05 & fdr_atac < 0.05 ~ "Concordant (both FDR<0.05)",
        concordant  & fdr_rna < 0.05                    ~ "Concordant (DEG FDR<0.05)",
        concordant                                       ~ "Concordant (nominal)",
        !concordant & fdr_rna < 0.05                    ~ "Discordant (DEG FDR<0.05)",
        TRUE                                             ~ "Discordant (nominal)"
      ), levels = c(
        "Concordant (both FDR<0.05)",
        "Concordant (DEG FDR<0.05)",
        "Concordant (nominal)",
        "Discordant (DEG FDR<0.05)",
        "Discordant (nominal)"
      ))
    )
  
  conc_colors <- c(
    "Concordant (both FDR<0.05)" = "#009B77",
    "Concordant (DEG FDR<0.05)"  = "#88B04B",
    "Concordant (nominal)"       = "#C8E6C9",
    "Discordant (DEG FDR<0.05)"  = "#C3447A",
    "Discordant (nominal)"       = "#F8BBD0"
  )
  
  # How many panels needed — only plot cell types with enough points
  cts_for_plot <- conc_df |>
    group_by(cell_type) |>
    filter(n() >= 3) |>
    pull(cell_type) |>
    unique()
  
  if (length(cts_for_plot) > 0) {
    p_conc <- ggplot(conc_df |> filter(cell_type %in% cts_for_plot),
                     aes(log2FC_atac, log2FC_rna, color = hit_class)) +
      geom_hline(yintercept = 0, color = "grey55", linewidth = 0.35) +
      geom_vline(xintercept = 0, color = "grey55", linewidth = 0.35) +
      geom_abline(slope = 1, intercept = 0,
                  linetype = "dashed", color = "grey75", linewidth = 0.4) +
      geom_point(size = 2, alpha = 0.8) +
      geom_text_repel(aes(label = label_gene), size = 2.8,
                      max.overlaps = 20, segment.color = "grey65",
                      min.segment.length = 0.2, na.rm = TRUE) +
      scale_color_manual(values = conc_colors, name = NULL) +
      facet_wrap(~ cell_type, scales = "free",
                 ncol = min(4, length(cts_for_plot))) +
      labs(
        title    = "Chromatin accessibility and gene expression concordance at linked loci",
        subtitle = "KO vs Ctrl (male) | DA peaks linked within 500 kb to DEGs",
        x        = "log\u2082FC peak accessibility (snATAC-seq, nominal p<0.05)",
        y        = "log\u2082FC gene expression (snRNA-seq, FDR<0.05)"
      ) +
      theme_minimal(base_size = 11) +
      theme(
        plot.title       = element_text(face = "bold", size = 12),
        plot.subtitle    = element_text(size = 10, color = "grey40",
                                        margin = margin(b = 6)),
        strip.text       = element_text(face = "bold", size = 9),
        legend.position  = "bottom",
        panel.grid.minor = element_blank(),
        axis.line        = element_line(color = "black", linewidth = 0.4),
        axis.ticks       = element_line(color = "black", linewidth = 0.3)
      ) +
      guides(color = guide_legend(nrow = 2, override.aes = list(size = 3.5)))
    
    n_ct_plot <- length(cts_for_plot)
    save_plot(p_conc, "DA_DEG_concordance",
              width  = min(18, 5 * ceiling(n_ct_plot / 2)),
              height = 4.5 * ceiling(n_ct_plot / min(4, n_ct_plot)))
    message("Concordance plot saved: ", n_ct_plot, " cell type panels")
  } else {
    message("Not enough cell types with >=3 overlap pairs to plot.")
  }
  
} else {
  message("Fewer than 5 overlap pairs — skipping concordance scatter.")
  message("Proceed to chromVAR (Step 5) — that does not require peak-gene links.")
}

message("\n=== Step 4 complete. Proceed to Step 5 (chromVAR). ===")

# ══════════════════════════════════════════════════════════════════════════════
# STEP 5: betterChromVAR — PER-CELL TF ACTIVITY DEVIATIONS
# ══════════════════════════════════════════════════════════════════════════════

DefaultAssay(obj_male) <- "peaks"

# Build count matrix — force standard dgCMatrix
peak_counts_mat <- as(
  GetAssayData(obj_male, assay = "peaks", layer = "counts"),
  "CsparseMatrix"
)

# Build GRanges, keep only standard chromosomes
peak_ranges <- StringToGRanges(rownames(peak_counts_mat))
main_chroms  <- c(paste0("chr", 1:19), "chrX", "chrY")
keep_chr     <- as.character(seqnames(peak_ranges)) %in% main_chroms
peak_counts_mat <- peak_counts_mat[keep_chr, ]
peak_ranges     <- peak_ranges[keep_chr]
seqlevels(peak_ranges, pruning.mode = "coarse") <- main_chroms

# Build SE
se <- SummarizedExperiment(
  assays    = list(counts = peak_counts_mat),
  rowRanges = peak_ranges
)

# Add GC bias
se <- addGCBias(se, genome = BSgenome.Mmusculus.UCSC.mm10)

# ── FIX: remove peaks with NA bias before anything else ──────────────────────
bias_vals <- rowData(se)$bias
na_bias   <- is.na(bias_vals)
message("Peaks with NA GC bias: ", sum(na_bias), " of ", nrow(se), " — removing")
se <- se[!na_bias, ]
message("SE after NA removal: ", nrow(se), " peaks x ", ncol(se), " cells")

# ── FIX: also remove peaks with zero variance across cells (all-zero rows) ───
peak_sums <- Matrix::rowSums(assay(se, "counts"))
zero_peaks <- peak_sums == 0
if (any(zero_peaks)) {
  message("Removing ", sum(zero_peaks), " all-zero peaks")
  se <- se[!zero_peaks, ]
}

# Rebuild pwm_filtered — mouse motifs filtered to expressed genes
pwm_set <- getMatrixSet(JASPAR2020,
                        opts = list(species = 10090, all_versions = FALSE))

# Filter to expressed genes using gene activity matrix
rna_counts      <- GetAssayData(rna, assay = "RNA", layer = "counts")
expressed_genes <- rownames(rna_counts)[rowSums(rna_counts > 0) >= 10]
message("Expressed genes: ", length(expressed_genes))

tf_names_all <- sapply(pwm_set@listData, function(x) x@name)
keep_motif   <- sapply(tf_names_all, function(tf) {
  parts <- sub("\\(.*\\)", "", unlist(strsplit(tf, "::"))) |> trimws()
  any(tolower(parts) %in% tolower(expressed_genes))
})
pwm_filtered <- pwm_set[keep_motif]
message("Motifs retained: ", length(pwm_filtered), " of ", length(pwm_set))

motif_matches <- matchMotifs(
  pwms    = pwm_filtered,
  subject = se,
  genome  = BSgenome.Mmusculus.UCSC.mm10
)
message("Motif matches: ", ncol(motif_matches), " motifs")

# Run betterChromVAR one-shot wrapper
dev <- tryCatch(
  betterChromVAR(se, motif_matches),
  error = function(e) {
    message("betterChromVAR failed: ", e$message)
    message("Falling back to standard computeDeviations...")
    bg  <- getBackgroundPeaks(se, niterations = 200)
    computeDeviations(se, annotations = motif_matches, background_peaks = bg)
  }
)

message("Deviations computed: ", nrow(dev), " motifs x ", ncol(dev), " cells")

# Extract z-score matrix
dev_scores <- tryCatch(
  deviationScores(dev),
  error = function(e) assay(dev, "z")
)

# Replace any remaining NaN/Inf with 0
dev_scores[is.nan(dev_scores)]  <- 0
dev_scores[is.infinite(dev_scores)] <- 0

# Align to obj_male cell order
shared_cells <- intersect(colnames(obj_male), colnames(dev_scores))
message("Cells with scores: ", length(shared_cells), " of ", ncol(obj_male))

obj_male[["chromvar"]] <- CreateAssayObject(
  data = dev_scores[, colnames(obj_male), drop = FALSE]
)
message("chromvar assay added: ", nrow(dev_scores), " motifs x ", ncol(obj_male), " cells")

#saveRDS(obj_male, "scATAC_06_male_chromvar.rds")

# ══════════════════════════════════════════════════════════════════════════════
# STEP 5b: DIFFERENTIAL TF ACTIVITY PER CELL TYPE
# ══════════════════════════════════════════════════════════════════════════════

DefaultAssay(obj_male) <- "chromvar"
Idents(obj_male)       <- "cell_type"
dir.create("chromVAR_results", showWarnings = FALSE)

for (ct in levels(obj_male$cell_type)) {
  ct_safe <- gsub("[^A-Za-z0-9]", "_", ct)
  outfile <- file.path("chromVAR_results", paste0(ct_safe, "__TFactivity_KOvCtrl.csv"))
  if (file.exists(outfile)) { message("Skipping chromVAR: ", ct); next }
  
  cells_ko   <- rownames(obj_male@meta.data[obj_male$cell_type == ct & obj_male$genotype == "KO",   ])
  cells_ctrl <- rownames(obj_male@meta.data[obj_male$cell_type == ct & obj_male$genotype == "Ctrl", ])
  
  if (length(cells_ko) < 10 | length(cells_ctrl) < 10) {
    message("  Skip [", ct, "] — too few cells (",
            length(cells_ko), " KO / ", length(cells_ctrl), " Ctrl)")
    next
  }
  
  res <- tryCatch(
    FindMarkers(obj_male,
                ident.1         = cells_ko,
                ident.2         = cells_ctrl,
                only.pos        = FALSE,
                mean.fxn        = rowMeans,
                fc.name         = "avg_diff",
                test.use        = "wilcox",
                min.pct         = 0,
                logfc.threshold = 0),
    error = function(e) { message("  Error [", ct, "]: ", e$message); NULL }
  )
  
  if (!is.null(res)) {
    res$p_val_fdr  <- p.adjust(res$p_val, method = "BH")
    motif_ids      <- rownames(res)
    res$motif_name <- sapply(motif_ids, function(id)
      tryCatch(TFBSTools::name(pwm_filtered[[id]]), error = function(e) id))
    write.csv(res, outfile)
    message("  chromVAR [", ct, "]: ",
            sum(res$p_val_fdr < 0.05, na.rm = TRUE), " sig TFs (FDR<0.05) | ",
            sum(res$p_val     < 0.05, na.rm = TRUE), " nominal")
  }
}

message("\nchromVAR_results files written: ",
        length(list.files("chromVAR_results", pattern = "\\.csv$")))

# ══════════════════════════════════════════════════════════════════════════════
# STEP 6: DOUBLE-HIT TF ANALYSIS
# ══════════════════════════════════════════════════════════════════════════════

chromvar_files  <- list.files("chromVAR_results", pattern = "\\.csv$", full.names = TRUE)
double_hit_list <- list()

for (cv_f in chromvar_files) {
  nm     <- gsub("\\.csv$", "", basename(cv_f))
  ct_key <- gsub("__TFactivity_KOvCtrl$", "", nm)
  ct     <- name_map[ct_key]
  
  if (is.na(ct)) {
    message("Unmapped chromVAR key: '", ct_key, "' — check name_map")
    next
  }
  if (!ct %in% shared_cts) {
    message("  Skip (not in shared_cts): ", ct)
    next
  }
  
  ct_safe  <- gsub("[^A-Za-z0-9]", "_", ct)
  deg_file <- file.path("DEG_results", paste0(ct_safe, "__KOvCtrl_male_DEG.csv"))
  if (!file.exists(deg_file)) {
    message("  Skip (no DEG file): ", ct)
    next
  }
  
  cv_df  <- read.csv(cv_f,     row.names = 1)
  deg_df <- read.csv(deg_file, row.names = 1)
  
  # Nominal TF activity hits
  cv_nom <- cv_df[!is.na(cv_df$p_val) & cv_df$p_val < 0.05, ]
  if (nrow(cv_nom) == 0) {
    message("  Skip [", ct, "] — no nominal chromVAR hits")
    next
  }
  cv_nom$tf_gene <- clean_tf(cv_nom$motif_name)
  
  dh <- cv_nom |>
    left_join(data.frame(tf_gene    = rownames(deg_df),
                         log2FC_rna = deg_df$avg_log2FC,
                         fdr_rna    = deg_df$p_val_fdr),
              by = "tf_gene") |>
    filter(!is.na(log2FC_rna)) |>
    mutate(cell_type  = ct,
           double_hit = !is.na(fdr_rna) & fdr_rna < 0.05,
           concordant = sign(avg_diff) == sign(log2FC_rna))
  
  if (nrow(dh) > 0) {
    double_hit_list[[ct]] <- dh
    message("  [", ct, "] double-hit pairs: ", nrow(dh),
            " (", sum(dh$double_hit & dh$concordant), " concordant FDR<0.05)")
  }
}

# ── Results ───────────────────────────────────────────────────────────────────
if (length(double_hit_list) > 0) {
  
  double_hit_df <- bind_rows(double_hit_list)
  write.csv(double_hit_df, "TF_double_hit.csv", row.names = FALSE)
  
  message("\nTop double-hit TFs (concordant, both FDR<0.05):")
  top_candidates <- double_hit_df |>
    filter(double_hit, concordant) |>
    arrange(p_val_fdr) |>
    select(cell_type, tf_gene, avg_diff, log2FC_rna, p_val_fdr, fdr_rna)
  print(head(top_candidates, 20))
  
  # ── Double-hit scatter ────────────────────────────────────────────────────
  dh_plot_df <- double_hit_df |>
    filter(!is.na(log2FC_rna), !is.na(avg_diff)) |>
    mutate(
      label       = ifelse(double_hit & concordant, tf_gene, NA),
      point_class = factor(case_when(
        double_hit &  concordant ~ "Both sig. \u2014 concordant",
        double_hit & !concordant ~ "Both sig. \u2014 discordant",
        TRUE                     ~ "Chromatin only (nominal)"
      ), levels = c("Both sig. \u2014 concordant",
                    "Both sig. \u2014 discordant",
                    "Chromatin only (nominal)"))
    )
  
  dh_colors <- c("Both sig. \u2014 concordant" = "#009B77",
                 "Both sig. \u2014 discordant" = "#C3447A",
                 "Chromatin only (nominal)"    = "grey75")
  
  cts_to_plot <- dh_plot_df |>
    group_by(cell_type) |>
    filter(n() >= 2) |>
    pull(cell_type) |>
    unique()
  
  if (length(cts_to_plot) > 0) {
    p_dh <- ggplot(dh_plot_df |> filter(cell_type %in% cts_to_plot),
                   aes(avg_diff, log2FC_rna, color = point_class)) +
      geom_hline(yintercept = 0, color = "grey55", linewidth = 0.35) +
      geom_vline(xintercept = 0, color = "grey55", linewidth = 0.35) +
      geom_point(size = 2.2, alpha = 0.8) +
      geom_text_repel(aes(label = label), size = 2.8, max.overlaps = 25,
                      segment.color = "grey65", min.segment.length = 0.2,
                      na.rm = TRUE) +
      scale_color_manual(values = dh_colors, name = NULL) +
      facet_wrap(~ cell_type, scales = "free",
                 ncol = min(4, length(cts_to_plot))) +
      labs(
        title    = "Double-hit TFs: concordant chromatin activity and expression changes",
        subtitle = "KO vs Ctrl (male) | x = chromVAR deviation | y = log\u2082FC expression (RNA)",
        x        = "Differential TF activity (chromVAR avg deviation, KO \u2212 Ctrl)",
        y        = "log\u2082FC TF expression (snRNA-seq)"
      ) +
      theme_minimal(base_size = 11) +
      theme(
        plot.title       = element_text(face = "bold", size = 12),
        plot.subtitle    = element_text(size = 9.5, color = "grey40",
                                        margin = margin(b = 6)),
        strip.text       = element_text(face = "bold", size = 9),
        legend.position  = "bottom",
        panel.grid.minor = element_blank(),
        axis.line        = element_line(color = "black", linewidth = 0.4),
        axis.ticks       = element_line(color = "black", linewidth = 0.3)
      ) +
      guides(color = guide_legend(nrow = 1, override.aes = list(size = 3.5)))
    
    save_plot(p_dh, "TF_double_hit_scatter",
              width  = min(18, 5 * ceiling(length(cts_to_plot) / 2)),
              height = 4.5 * ceiling(length(cts_to_plot) / min(4, length(cts_to_plot))))
    message("Double-hit scatter saved.")
  }
  
} else {
  double_hit_df <- data.frame()
  message("No double-hit TFs found.")
  message("Check: list.files('chromVAR_results') and list.files('DEG_results')")
}

saveRDS(obj_male, "scATAC_06_male_doubleTF.rds")
message("\n=== Steps 5b + 6 complete ===")


# ══════════════════════════════════════════════════════════════════════════════
# STEP 7: COVERAGE PLOTS — interleaved cell type × genotype ordering
# ══════════════════════════════════════════════════════════════════════════════

DefaultAssay(obj_male) <- "peaks"

# Build interleaved ct_geno factor: CT1_Ctrl, CT1_KO, CT2_Ctrl, CT2_KO ...
ct_levels <- levels(factor(obj_male$cell_type,
                           levels = atac_new_order[atac_new_order %in%
                                                     unique(obj_male$cell_type)]))

interleaved_levels <- as.vector(rbind(
  paste(ct_levels, "Ctrl", sep = " \u2014 "),
  paste(ct_levels, "KO",   sep = " \u2014 ")
))

obj_male$ct_geno <- paste(obj_male$cell_type, obj_male$genotype, sep = " \u2014 ")
obj_male$ct_geno <- factor(
  obj_male$ct_geno,
  levels = interleaved_levels[interleaved_levels %in% unique(obj_male$ct_geno)]
)

# Build interleaved color vector — same cell type color, Ctrl lighter, KO full
ct_geno_colors <- unlist(lapply(ct_levels, function(ct) {
  base_col <- atac_colors[ct]
  if (is.na(base_col)) base_col <- "grey60"
  ctrl_col <- scales::alpha(base_col, 0.45)
  ko_col   <- base_col
  setNames(c(ctrl_col, ko_col),
           c(paste(ct, "Ctrl", sep = " \u2014 "),
             paste(ct, "KO",   sep = " \u2014 ")))
}))
# Keep only levels that exist
ct_geno_colors <- ct_geno_colors[names(ct_geno_colors) %in% levels(obj_male$ct_geno)]

# Gene selection — concordant DA+DEG pairs, fallback to markers
if (nrow(overlap_df) > 0) {
  top_genes_coverage <- overlap_df |>
    filter(!is.na(log2FC_rna), !is.na(log2FC_atac),
           sign(log2FC_rna) == sign(log2FC_atac),
           fdr_rna < 0.05) |>
    arrange(fdr_rna) |>
    pull(gene) |>
    unique() |>
    head(12)
} else {
  top_genes_coverage <- c("Pvalb", "Sst", "Olig2", "Mbp",
                          "P2ry12", "Pdgfra", "Gfap", "Slc17a7")
}

message("Coverage plots for: ", paste(top_genes_coverage, collapse = ", "))
dir.create("coverage_plots", showWarnings = FALSE)
n_tracks <- nlevels(obj_male$ct_geno)

for (gene in top_genes_coverage) {
  tryCatch({
    p_cov <- CoveragePlot(
      object            = obj_male,
      region            = gene,
      features          = gene,
      assay             = "peaks",
      expression.assay  = "RNA",
      group.by          = "ct_geno",
      annotation        = TRUE,
      peaks             = TRUE,
      links             = TRUE,
      extend.upstream   = 2000,
      extend.downstream = 2000
    ) &
      scale_fill_manual(values = ct_geno_colors)
    
    ggsave(
      file.path("coverage_plots", paste0(gene, "_coverage.pdf")),
      p_cov,
      width  = 14,
      height = max(8, n_tracks * 0.55)
    )
    message("  Saved: ", gene)
  }, error = function(e) message("  FAILED [", gene, "]: ", e$message))
}

# ══════════════════════════════════════════════════════════════════════════════
# STEP 8: INTEGRATION SUMMARY TABLE + FINAL SAVE
# ══════════════════════════════════════════════════════════════════════════════

summary_integration <- lapply(levels(obj_male$cell_type), function(ct) {
  ct_safe <- gsub("[^A-Za-z0-9]", "_", ct)
  
  n_da  <- tryCatch({
    f <- file.path("DA_results", paste0(ct_safe, "__KOvCtrl_male.csv"))
    if (file.exists(f)) sum(read.csv(f, row.names=1)$p_val_fdr < 0.05, na.rm=TRUE) else NA_integer_
  }, error = function(e) NA_integer_)
  
  n_deg <- tryCatch({
    f <- file.path("DEG_results", paste0(ct_safe, "__KOvCtrl_male_DEG.csv"))
    if (file.exists(f)) sum(read.csv(f, row.names=1)$p_val_fdr < 0.05, na.rm=TRUE) else NA_integer_
  }, error = function(e) NA_integer_)
  
  n_overlap <- if (nrow(overlap_df) > 0)
    sum(overlap_df$cell_type == ct &
          !is.na(overlap_df$log2FC_rna) &
          sign(overlap_df$log2FC_atac) == sign(overlap_df$log2FC_rna), na.rm = TRUE)
  else NA_integer_
  
  n_dh <- if (nrow(double_hit_df) > 0)
    sum(double_hit_df$cell_type == ct &
          double_hit_df$double_hit &
          double_hit_df$concordant, na.rm = TRUE)
  else NA_integer_
  
  data.frame(cell_type         = ct,
             DA_peaks_KOvCtrl  = n_da,
             DEGs_KOvCtrl      = n_deg,
             concordant_linked = n_overlap,
             double_hit_TFs    = n_dh)
}) |> bind_rows() |> arrange(desc(concordant_linked))

print(tibble::as_tibble(summary_integration), n = Inf)
write.csv(summary_integration, "integration_summary.csv", row.names = FALSE)

saveRDS(obj_male, "scATAC_07_male_integrated.rds")
#obj_male <- readRDS("scATAC_07_male_integrated.rds")

################################# Venn diagram + Concordance --------------------
# ══════════════════════════════════════════════════════════════════════════════
# DA peaks × DEGs overlap + directional concordance per cell cluster
#
# DA filter:  nominal p<0.05 + |log2FC|>0.1 (single-replicate design)
# DEG filter: BH FDR<0.05
# ══════════════════════════════════════════════════════════════════════════════

setwd("~/Downloads/Seurat_scATAC-seq")

library(ggvenn)
library(ggplot2)
library(patchwork)
library(dplyr)

# ── 0. Rebuild session variables ─────────────────────────────────────────────
if (!exists("links_sig")) {
  links_df  <- read.csv("peak_gene_links.csv")
  links_sig <- links_df[!is.na(links_df$pvalue) & links_df$pvalue < 0.05, ]
  message("links_sig rows: ", nrow(links_sig))
}

if (!exists("name_map")) {
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
}

if (!exists("shared_cts")) {
  da_keys <- gsub("__KOvCtrl_male\\.csv$", "",
                  basename(list.files("DA_results", pattern = "__KOvCtrl_male\\.csv$")))
  da_cts  <- na.omit(name_map[da_keys])
  shared_cts <- names(da_cts)[sapply(names(da_cts), function(ct) {
    ct_safe  <- gsub("[^A-Za-z0-9]", "_", ct)
    file.exists(file.path("DEG_results", paste0(ct_safe, "__KOvCtrl_male_DEG.csv")))
  })]
  message("Shared cell types found: ", length(shared_cts))
}

if (exists("diag_df")) {
  shared_cts <- diag_df$cell_type[
    !is.na(diag_df$n_deg_sig) & diag_df$n_deg_sig > 0 & diag_df$has_deg_file
  ]
  message("Rebuilt shared_cts from diag_df: ", length(shared_cts))
}

da_files <- list.files("DA_results", pattern = "__KOvCtrl_male\\.csv$", full.names = TRUE)

# ── Shared DA filter (used identically everywhere) ────────────────────────────
get_da_sig_peaks <- function(da_df) {
  rownames(da_df)[
    !is.na(da_df$p_val) & da_df$p_val < 0.05 & abs(da_df$avg_log2FC) > 0.1
  ]
}

# ── 1. Build gene sets + run Venn + overlap summary + concordance ─────────────
venn_list         <- list()
summary_overlap   <- list()
da_direction_list <- list()

for (da_f in da_files) {
  ct_key <- gsub("__KOvCtrl_male\\.csv$", "", basename(da_f))
  ct     <- name_map[ct_key]
  if (is.na(ct) | !ct %in% shared_cts) next
  
  ct_safe  <- gsub("[^A-Za-z0-9]", "_", ct)
  deg_file <- file.path("DEG_results", paste0(ct_safe, "__KOvCtrl_male_DEG.csv"))
  if (!file.exists(deg_file)) next
  
  da_df  <- read.csv(da_f,     row.names = 1)
  deg_df <- read.csv(deg_file, row.names = 1)
  
  # Consistent filter throughout
  da_sig_peaks <- get_da_sig_peaks(da_df)
  da_genes     <- unique(links_sig$gene[links_sig$peak %in% da_sig_peaks])
  deg_genes    <- rownames(deg_df)[!is.na(deg_df$p_val_fdr) & deg_df$p_val_fdr < 0.05]
  overlap      <- intersect(da_genes, deg_genes)
  
  message(sprintf("  [%s]  DA-linked genes: %d | DEGs: %d | overlap: %d",
                  ct, length(da_genes), length(deg_genes), length(overlap)))
  
  # ── Venn ────────────────────────────────────────────────────────────────────
  if (length(da_genes) > 0 & length(deg_genes) > 0) {
    venn_list[[ct]] <- list(
      "DA-linked genes" = da_genes,
      "DEGs"            = deg_genes
    )
  }
  
  # ── Overlap summary ─────────────────────────────────────────────────────────
  if (length(overlap) > 0) {
    summary_overlap[[ct]] <- data.frame(
      gene             = overlap,
      log2FC_rna       = deg_df[overlap, "avg_log2FC"],
      fdr_rna          = deg_df[overlap, "p_val_fdr"],
      stringsAsFactors = FALSE
    ) |>
      mutate(
        direction      = ifelse(log2FC_rna > 0, "Up in KO", "Down in KO"),
        cell_type      = ct,
        n_da_sig_peaks = length(da_sig_peaks),
        n_da_genes     = length(da_genes),
        n_deg_sig      = length(deg_genes),
        n_overlap      = length(overlap),
        pct_da_in_deg  = round(length(overlap) / length(da_genes)  * 100, 1),
        pct_deg_in_da  = round(length(overlap) / length(deg_genes) * 100, 1)
      ) |>
      select(cell_type, gene, direction, log2FC_rna, fdr_rna,
             n_da_sig_peaks, n_da_genes, n_deg_sig, n_overlap,
             pct_da_in_deg, pct_deg_in_da)
  }
  
  # ── Concordance ─────────────────────────────────────────────────────────────
  if (length(da_sig_peaks) > 0) {
    da_peak_df <- data.frame(
      peak         = da_sig_peaks,
      log2FC_atac  = da_df[da_sig_peaks, "avg_log2FC"],
      da_direction = ifelse(da_df[da_sig_peaks, "avg_log2FC"] > 0, "More open", "Less open"),
      stringsAsFactors = FALSE
    )
    
    peak_gene    <- links_sig[links_sig$peak %in% da_sig_peaks, c("peak", "gene")]
    da_peak_gene <- merge(da_peak_df, peak_gene, by = "peak")
    da_peak_gene <- da_peak_gene[da_peak_gene$gene %in% deg_genes, ]
    
    if (nrow(da_peak_gene) > 0) {
      da_peak_gene$log2FC_rna    <- deg_df[da_peak_gene$gene, "avg_log2FC"]
      da_peak_gene$rna_direction <- ifelse(da_peak_gene$log2FC_rna > 0, "Up in KO", "Down in KO")
      da_peak_gene$cell_type     <- ct
      da_peak_gene$concordant    <- (
        (da_peak_gene$da_direction == "More open" & da_peak_gene$rna_direction == "Up in KO") |
          (da_peak_gene$da_direction == "Less open" & da_peak_gene$rna_direction == "Down in KO")
      )
      da_peak_gene$concordance   <- ifelse(da_peak_gene$concordant, "Concordant", "Discordant")
      da_direction_list[[ct]]    <- da_peak_gene
    }
  }
}

# ── 2. Venn diagrams ──────────────────────────────────────────────────────────
message("\nCell types with data for Venn: ", length(venn_list))

venn_plots <- lapply(names(venn_list), function(ct) {
  ggvenn(
    venn_list[[ct]],
    fill_color    = c("#4E84C4", "#DD4124"),
    fill_alpha    = 0.45,
    stroke_color  = "white",
    text_size     = 3.5,
    set_name_size = 3
  ) +
    ggtitle(ct) +
    theme(
      plot.title  = element_text(face = "bold", size = 9, hjust = 0.5),
      plot.margin = margin(4, 4, 4, 4)
    )
})

n     <- length(venn_plots)
ncols <- min(4, n)
nrows <- ceiling(n / ncols)

p_venn <- wrap_plots(venn_plots, ncol = ncols) +
  plot_annotation(
    title    = "DA peaks × DEGs overlap per cell cluster",
    subtitle = "KO vs Ctrl (male) | DA: nominal p<0.05, |log2FC|>0.1 (single-replicate design) | DEGs: BH FDR<0.05",
    theme    = theme(
      plot.title    = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 10, color = "grey40")
    )
  )

ggsave("DA_DEG_venn_per_cluster.pdf", p_venn,
       width = ncols * 4, height = nrows * 4)
message("Saved: DA_DEG_venn_per_cluster.pdf")

# ── 3. Overlap summary export ─────────────────────────────────────────────────
overlap_long <- bind_rows(summary_overlap)

overlap_wide <- overlap_long |>
  group_by(cell_type, n_da_sig_peaks, n_da_genes, n_deg_sig,
           n_overlap, pct_da_in_deg, pct_deg_in_da) |>
  summarise(
    n_up_in_KO    = sum(direction == "Up in KO"),
    n_down_in_KO  = sum(direction == "Down in KO"),
    overlap_genes = paste(sort(gene), collapse = ";"),
    .groups = "drop"
  ) |>
  arrange(desc(n_overlap))

write.csv(overlap_wide, "DA_DEG_overlap_summary.csv", row.names = FALSE)
write.csv(overlap_long, "DA_DEG_overlap_genes.csv",   row.names = FALSE)
message("Saved: DA_DEG_overlap_summary.csv / DA_DEG_overlap_genes.csv")

overlap_wide |>
  select(cell_type, n_da_sig_peaks, n_da_genes, n_deg_sig,
         n_overlap, pct_da_in_deg, pct_deg_in_da,
         n_up_in_KO, n_down_in_KO) |>
  print(n = Inf)

# ── 4. Concordance analysis ───────────────────────────────────────────────────
da_rna_df <- bind_rows(da_direction_list)

if (nrow(da_rna_df) == 0) {
  message("No DA peak–DEG pairs found for concordance analysis.")
} else {
  
  concordance_summary <- da_rna_df |>
    group_by(cell_type, concordance) |>
    summarise(n = n(), .groups = "drop") |>
    tidyr::pivot_wider(names_from = concordance, values_from = n, values_fill = 0) |>
    mutate(
      total          = Concordant + Discordant,
      pct_concordant = round(100 * Concordant / total, 1)
    ) |>
    arrange(desc(pct_concordant))
  
  print(concordance_summary)
  
  direction_crosstab <- da_rna_df |>
    group_by(cell_type, da_direction, rna_direction) |>
    summarise(n_peak_gene_pairs = n(), .groups = "drop") |>
    arrange(cell_type, da_direction, rna_direction)
  
  print(direction_crosstab)
  
  write.csv(da_rna_df,           "DA_DEG_direction_per_peak_gene.csv", row.names = FALSE)
  write.csv(concordance_summary, "DA_DEG_concordance_summary.csv",     row.names = FALSE)
  write.csv(direction_crosstab,  "DA_DEG_direction_crosstab.csv",      row.names = FALSE)
  message("Saved direction concordance files")
  
  p_conc <- ggplot(
    da_rna_df,
    aes(x = reorder(cell_type, concordant), fill = concordance)
  ) +
    geom_bar(position = "fill") +
    scale_y_continuous(labels = scales::percent) +
    scale_fill_manual(values = c("Concordant" = "#2166ac", "Discordant" = "#d73027")) +
    coord_flip() +
    labs(
      title    = "DA–RNA concordance per cell type",
      subtitle = "Concordant = more open + up, or less open + down in KO\nDA: nominal p<0.05 | DEGs: BH FDR<0.05",
      x        = NULL,
      y        = "Fraction of DA peak–DEG pairs",
      fill     = NULL
    ) +
    theme_classic(base_size = 11) +
    theme(plot.title = element_text(face = "bold"))
  
  ggsave("DA_DEG_concordance_barplot.pdf", p_conc, width = 9, height = 6)
  message("Saved: DA_DEG_concordance_barplot.pdf")
}
############################# Check Dnmt1 -------------------------------------------------
# In R — targeted test for Dnmt1 specifically
cells_ko   <- rownames(rna@meta.data[rna$cell_type == "PV+ interneurons" & rna$Genotype == "KO", ])
cells_ctrl <- rownames(rna@meta.data[rna$cell_type == "PV+ interneurons" & rna$Genotype == "Ctrl", ])

# Option B: use group.by + subset.ident to restrict to PV cells only
dnmt1_test <- FindMarkers(
  rna,
  ident.1         = "KO",
  ident.2         = "Ctrl",
  group.by        = "genotype",
  subset.ident    = "PV+ interneurons",   # restricts to PV cells via cell_type ident
  features        = "Dnmt1",
  test.use        = "wilcox",
  min.pct         = 0,
  logfc.threshold = 0
)
dnmt1_test$p_val_fdr <- p.adjust(dnmt1_test$p_val, method = "BH")
print(dnmt1_test)

# Check raw expression directly too
counts  <- GetAssayData(rna, layer = "counts")
meta    <- rna@meta.data
pv_ko   <- meta$cell_type == "PV+ interneurons" & meta$genotype == "KO"
pv_ctrl <- meta$cell_type == "PV+ interneurons" & meta$genotype == "Ctrl"
dnmt1_c <- counts["Dnmt1", ]

message("Dnmt1 pct KO:    ", round(mean(dnmt1_c[pv_ko]   > 0), 3))
message("Dnmt1 pct Ctrl:  ", round(mean(dnmt1_c[pv_ctrl] > 0), 3))
message("Dnmt1 mean KO:   ", round(mean(dnmt1_c[pv_ko]),  4))
message("Dnmt1 mean Ctrl: ", round(mean(dnmt1_c[pv_ctrl]), 4))