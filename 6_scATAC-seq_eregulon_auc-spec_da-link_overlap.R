library(dplyr)
library(tidyr)

setwd("~/Downloads/Seurat_scATAC-seq")

SCENIC_DIR  <- "scenicplus_output/scenicplus/scenicplus_further_analysis"
TF_DEG_FILE <- file.path(SCENIC_DIR, "TF_DEG_overlap/TF_DEG_overlap_permutation_KO_vs_Ctrl.csv")
MULTI_FILE  <- file.path(SCENIC_DIR, "TF_DEG_overlap/TF_DEG_overlap_permutation_KO_vs_Ctrl_multi_eRegulon_spec.csv")
DA_DEG_FILE <- "DA_DEG_direction_per_peak_gene.csv"

# ── Load ──────────────────────────────────────────────────────
da_deg  <- read.csv(DA_DEG_FILE)
tf_deg  <- read.csv(TF_DEG_FILE)
multi   <- read.csv(MULTI_FILE)

# Verify actual cell type values match between files
message("DA_DEG cell types:  ", paste(sort(unique(da_deg$cell_type)), collapse = " | "))
message("TF_DEG cell types:  ", paste(sort(unique(tf_deg$Cell_Type)), collapse = " | "))
message("Multi  cell types:  ", paste(sort(unique(multi$Cell_Type)),  collapse = " | "))

# ── DA/DEG: one row per gene × cell_type (keep peak with strongest atac signal) ──
da_deg_clean <- da_deg %>%
  dplyr::group_by(gene, cell_type) %>%
  dplyr::slice_max(abs(log2FC_atac), n = 1, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::rename(Gene = gene, Cell_Type = cell_type)

# ── eRegulon: expand Overlap_Genes, keep directional where available ──────────
sig_dir <- tf_deg %>%
  dplyr::filter(Direction %in% c("Up", "Down"), Perm_FDR < 0.05)

sig_all <- tf_deg %>%
  dplyr::filter(Direction == "All", Perm_FDR < 0.05) %>%
  dplyr::anti_join(sig_dir, by = c("TF", "Cell_Type"))

eregulon_targets <- dplyr::bind_rows(sig_dir, sig_all) %>%
  dplyr::mutate(Gene = as.character(Overlap_Genes)) %>%
  tidyr::separate_rows(Gene, sep = ",\\s*") %>%
  dplyr::filter(!is.na(Gene), Gene != "", !grepl("ENSMUSG", Gene)) %>%
  dplyr::select(Gene, Cell_Type, TF, Direction,
                fdr_eregulon = Perm_FDR, enrichment = Enrichment_Ratio) %>%
  dplyr::group_by(Gene, Cell_Type, TF) %>%
  dplyr::slice_min(fdr_eregulon, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup()

# ── Co-regulator multi-eRegulon: same logic ───────────────────
sig_dir_m <- multi %>%
  dplyr::filter(Direction %in% c("Up", "Down"), Perm_FDR < 0.05)

sig_all_m <- multi %>%
  dplyr::filter(Direction == "All", Perm_FDR < 0.05) %>%
  dplyr::anti_join(sig_dir_m, by = c("TF", "Cell_Type"))

coreg_targets <- dplyr::bind_rows(sig_dir_m, sig_all_m) %>%
  dplyr::mutate(Gene = as.character(Overlap_Genes)) %>%
  tidyr::separate_rows(Gene, sep = ",\\s*") %>%
  dplyr::filter(!is.na(Gene), Gene != "", !grepl("ENSMUSG", Gene)) %>%
  dplyr::select(Gene, Cell_Type,
                Co_regulator = TF,
                fdr_coreg    = Perm_FDR) %>%
  dplyr::group_by(Gene, Cell_Type, Co_regulator) %>%
  dplyr::slice_min(fdr_coreg, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup()

message("DA/DEG pairs:         ", nrow(da_deg_clean), " (", length(unique(da_deg_clean$Gene)), " genes)")
message("eRegulon targets:     ", nrow(eregulon_targets), " (", length(unique(eregulon_targets$Gene)), " genes)")
message("Co-regulator targets: ", nrow(coreg_targets), " (", length(unique(coreg_targets$Gene)), " genes)")

# ── Triple join ───────────────────────────────────────────────
triple <- da_deg_clean %>%
  dplyr::inner_join(eregulon_targets, by = c("Gene", "Cell_Type"),
                    relationship = "many-to-many") %>%
  dplyr::left_join(coreg_targets,    by = c("Gene", "Cell_Type"),
                   relationship = "many-to-many") %>%
  dplyr::mutate(
    triple_hit       = !is.na(Co_regulator),
    n_evidence_tiers = 2L + as.integer(triple_hit)
  ) %>%
  dplyr::arrange(desc(triple_hit), desc(concordant), fdr_eregulon)

write.csv(triple, "triple_evidence_integration.csv", row.names = FALSE)
message("Written: triple_evidence_integration.csv  (", nrow(triple), " rows)")

# ── Summary ───────────────────────────────────────────────────
cat("\n=== DA + eRegulon (tier 2) ===\n")
cat("Combinations:", nrow(triple), "\n")
cat("Unique genes:", length(unique(triple$Gene)), "\n")
cat("Concordant:  ", sum(triple$concordant, na.rm = TRUE), "\n")

cat("\n=== Triple hits (tier 3) ===\n")
print(triple %>%
        dplyr::filter(triple_hit) %>%
        dplyr::count(Cell_Type, TF, Co_regulator, sort = TRUE))

cat("\n=== Concordant triple hits ===\n")
print(triple %>%
        dplyr::filter(triple_hit, concordant) %>%
        dplyr::select(Gene, Cell_Type, TF, Co_regulator,
                      log2FC_rna, log2FC_atac, concordance, fdr_eregulon) %>%
        dplyr::arrange(fdr_eregulon) %>%
        head(30))

# Full tier-2 concordant hits for supplementary
tier2_concordant <- triple %>%
  dplyr::filter(concordant) %>%
  dplyr::select(Gene, Cell_Type, TF, Co_regulator,
                log2FC_rna, log2FC_atac, da_direction, rna_direction,
                concordance, fdr_eregulon, enrichment) %>%
  dplyr::arrange(Cell_Type, fdr_eregulon)

write.csv(tier2_concordant, "triple_evidence_tier2_concordant.csv", row.names = FALSE)

cat("Tier-2 concordant hits per cell type:\n")
print(table(tier2_concordant$Cell_Type))

#########
library(ggplot2)
library(ggrepel)
library(dplyr)

tier2_plot <- tier2_concordant %>%
  dplyr::mutate(
    triple_hit = !is.na(Co_regulator),
    # All points get gene label; triple hits get TF annotation appended
    label = ifelse(
      triple_hit,
      paste0(Gene, " (", TF, " + ", Co_regulator, ")"),
      Gene
    )
  )

p <- ggplot(tier2_plot,
            aes(x = log2FC_rna, y = log2FC_atac,
                color = Cell_Type, shape = triple_hit)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey70", linewidth = 0.4) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted",
              color = "grey50", linewidth = 0.4) +
  geom_point(size = 3.5, alpha = 0.85, stroke = 0.8) +
  geom_text_repel(
    aes(label = label),
    size          = 2.8,
    fontface      = "italic",
    segment.size  = 0.25,
    segment.color = "grey60",
    box.padding   = 0.5,
    point.padding = 0.3,
    force         = 12,
    force_pull    = 0.5,
    max.overlaps  = Inf,       # label everything
    min.segment.length = 0,
    show.legend   = FALSE
  ) +
  scale_shape_manual(
    values = c("FALSE" = 16, "TRUE" = 18),
    labels = c("DA + eRegulon", "DA + eRegulon + co-regulator"),
    name   = "Evidence tier"
  ) +
  scale_color_manual(values = atac_colors, name = "Cell type") +
  labs(
    title    = "Chromatin–expression concordance: triple-evidence genes",
    subtitle = "Each point = one gene × cell type | diamonds = triple hits (DA + eRegulon + co-regulator)",
    x        = expression(log[2]*"FC RNA (KO vs Ctrl)"),
    y        = expression(log[2]*"FC ATAC (KO vs Ctrl)")
  ) +
  coord_equal(xlim = c(-2.2, 2.2), ylim = c(-2.8, 2.8)) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid     = element_line(color = "grey93"),
    plot.title     = element_text(face = "bold", size = 12),
    plot.subtitle  = element_text(size = 9, color = "grey40"),
    legend.position = "right",
    axis.text      = element_text(size = 9)
  )

ggsave("triple_evidence_concordance_scatter.pdf",
       p, width = 11, height = 9, limitsize = FALSE)
message("Saved.")
