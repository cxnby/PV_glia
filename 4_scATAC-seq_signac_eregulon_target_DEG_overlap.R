# ══════════════════════════════════════════════════════════════════════════════
# Requires in session from Steps 1–8:
#   obj_male, double_hit_df, overlap_df, pwm_filtered
#   atac_colors, atac_new_order, clean_tf, save_plot
# Files on disk:
#   TF_double_hit.csv, DA_DEG_overlap.csv
# ══════════════════════════════════════════════════════════════════════════════
BiocManager::install(c(
  "RSQLite",
  "AnnotationDbi",
  "GO.db",
  "org.Mm.eg.db",
  "clusterProfiler",
  "enrichplot"
), force = TRUE, ask = FALSE)

BiocManager::install(c(
  "TxDb.Mmusculus.UCSC.mm10.knownGene"
))

library(clusterProfiler)
library(org.Mm.eg.db)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tidyr)
library(patchwork)
library(scales)
library(ggtext)
library(GenomicFeatures)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)

setwd("~/Downloads/Seurat_scATAC-seq")

# Reload if starting fresh
if (!exists("double_hit_df")) {
  double_hit_df <- read.csv("TF_double_hit.csv")
}
if (!exists("overlap_df")) {
  overlap_df <- read.csv("DA_DEG_overlap.csv")
  if (nrow(overlap_df) == 0) overlap_df <- data.frame()
}

###### TF motifs sig. differentially opened or closed --------
# ══════════════════════════════════════════════════════════════════════════════
# STEP 10 (FALLBACK): MOTIF ENRICHMENT ON ALL DA PEAKS PER CELL TYPE
# Uses all KO-gained / KO-lost DA peaks, not just those linked to DEGs
# ══════════════════════════════════════════════════════════════════════════════

DefaultAssay(obj_male) <- "peaks"
dir.create("motif_enrichment_per_celltype", showWarnings = FALSE)

motif_results_list <- list()
da_files_male <- list.files("DA_results", pattern = "__KOvCtrl_male\\.csv$",
                            full.names = TRUE)

for (da_f in da_files_male) {
  
  ct_key <- gsub("__KOvCtrl_male\\.csv$", "", basename(da_f))
  ct     <- name_map[ct_key]
  if (is.na(ct)) { message("Unmapped: ", ct_key); next }
  ct_safe <- gsub("[^A-Za-z0-9]", "_", ct)
  
  da_df <- read.csv(da_f, row.names = 1)
  
  # Background = all tested peaks for this cell type that are in the object
  background_peaks <- intersect(rownames(da_df), rownames(obj_male))
  
  if (length(background_peaks) < 100) {
    message("  Skip [", ct, "] — too few background peaks")
    next
  }
  
  for (direction in c("gained", "lost")) {
    
    sign_val <- if (direction == "gained") 1 else -1
    
    # Nominal p<0.05 + directional — same threshold as your DA analysis
    target_peaks <- rownames(da_df)[
      !is.na(da_df$p_val) &
        da_df$p_val < 0.05 &
        sign(da_df$avg_log2FC) == sign_val
    ]
    target_peaks <- intersect(target_peaks, rownames(obj_male))
    
    if (length(target_peaks) < 10) {
      message("  Skip [", ct, " — ", direction, "] — ",
              length(target_peaks), " peaks")
      next
    }
    
    message("  FindMotifs [", ct, " — ", direction, "] — ",
            length(target_peaks), " peaks vs ",
            length(background_peaks), " background")
    
    enr <- tryCatch(
      FindMotifs(object     = obj_male,
                 features   = target_peaks,
                 background = background_peaks),
      error = function(e) { message("  Error: ", e$message); NULL }
    )
    
    if (!is.null(enr) && nrow(enr) > 0) {
      enr <- enr |>
        filter(observed > 0, p.adjust < 0.05) |>
        mutate(cell_type = ct, direction = direction)
      
      if (nrow(enr) > 0) {
        motif_results_list[[paste(ct_safe, direction, sep = "_")]] <- enr
        write.csv(enr,
                  file.path("motif_enrichment_per_celltype",
                            paste0(ct_safe, "_", direction, "_motifs.csv")),
                  row.names = FALSE)
        message("    Sig motifs: ", nrow(enr))
      } else {
        message("    No sig motifs after filtering")
      }
    }
  }
}

# ── Plot ──────────────────────────────────────────────────────────────────────
if (length(motif_results_list) > 0) {
  
  motif_all <- bind_rows(motif_results_list) |>
    mutate(
      neg_log10_fdr = -log10(p.adjust + 1e-300),
      direction     = factor(direction,
                             levels = c("gained", "lost"),
                             labels = c("Gained (KO more open)",
                                        "Lost (KO less open)"))
    )
  
  motif_top <- motif_all |>
    group_by(cell_type, direction) |>
    slice_min(p.adjust, n = 10) |>
    ungroup() |>
    mutate(
      panel_label = paste0(cell_type, "\n", direction),
      motif_label = factor(motif.name, levels = rev(unique(motif.name)))
    )
  
  n_panels    <- length(unique(motif_top$panel_label))
  dir_colors  <- c("Gained (KO more open)" = "#009B77",
                   "Lost (KO less open)"   = "#C3447A")
  
  p_motif <- ggplot(motif_top,
                    aes(x = neg_log10_fdr, y = motif_label,
                        color = direction, size = fold.enrichment)) +
    geom_point(alpha = 0.85) +
    geom_vline(xintercept = -log10(0.05), linetype = "dashed",
               color = "grey60", linewidth = 0.5) +
    scale_color_manual(values = dir_colors, name = NULL) +
    scale_size_continuous(range = c(2, 8), name = "Fold enrichment") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
    facet_wrap(~ panel_label, scales = "free_y",
               ncol = min(3, n_panels)) +
    labs(
      title    = "Motif enrichment in DA peaks (KO vs Ctrl)",
      subtitle = "Per cell type (male) | nominal p<0.05 DA peaks | FDR<0.05 motifs only",
      x        = "−log₁₀(FDR)", y = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title       = element_text(face = "bold", size = 12),
      plot.subtitle    = element_text(size = 9.5, color = "grey40"),
      strip.text       = element_text(face = "bold", size = 9),
      strip.background = element_rect(fill = "grey92", color = NA),
      legend.position  = "right",
      panel.grid.minor = element_blank(),
      axis.text.y      = element_text(size = 8),
      axis.line        = element_line(color = "black", linewidth = 0.4),
      axis.ticks       = element_line(color = "black", linewidth = 0.3)
    )
  
  p_motif <- p_motif +
    labs(x = expression(-log[10](FDR)), y = NULL) +
    scale_size_continuous(range = c(3, 3)) +   # fixed size — removes uninformative legend
    guides(size = "none")                       # drop size legend
  
  save_plot(p_motif, "motif_enrichment_DA_per_celltype_dotplot_v2",
            width  = max(14, 5 * ceiling(n_panels / 3)),
            height = max(8,  4 * ceiling(n_panels / min(3, n_panels))))
  
  # save_plot(p_motif, "motif_enrichment_DA_per_celltype_dotplot",
  #           width  = max(14, 5 * ceiling(n_panels / 3)),
  #           height = max(8,  4 * ceiling(n_panels / min(3, n_panels))))
  message("Dotplot saved.")
  
} else {
  message("No significant motifs found in any cell type.")
}

message("\n=== Step 10 complete. ===")


######### Load differential data for master regulators -------
# ── Load TF-DEG overlap ───────────────────────────────────────────────────────
tf_deg <- read.csv(file.path(further_dir,
                             "TF_DEG_overlap/TF_DEG_overlap_permutation_KO_vs_Ctrl.csv"))
message("TF-DEG overlap loaded: ", nrow(tf_deg), " rows")
message("Columns: ", paste(colnames(tf_deg), collapse = " | "))

# ── Identify gene + cell type + TF columns ────────────────────────────────────
gene_col <- grep("^gene$|^Gene$|^target$|^Target$", colnames(tf_deg),
                 ignore.case = TRUE, value = TRUE)[1]
tf_col   <- grep("^TF$|^tf$", colnames(tf_deg),
                 ignore.case = TRUE, value = TRUE)[1]
ct_col   <- grep("CellType|cell_type|celltype", colnames(tf_deg),
                 ignore.case = TRUE, value = TRUE)[1]
pval_col <- grep("fdr|padj|adj", colnames(tf_deg),
                 ignore.case = TRUE, value = TRUE)[1]

message("\nMapped columns:")
message("  gene  = ", gene_col)
message("  TF    = ", tf_col)
message("  CT    = ", ct_col)
message("  pval  = ", pval_col)

# ── Load DEG direction (up/down in KO) ────────────────────────────────────────
deg_dir_path <- file.path(further_dir, "DEG_KO_vs_Ctrl")
deg_all      <- read.csv(file.path(deg_dir_path,
                                   "DEG_all_celltypes_KO_vs_Ctrl.csv"))
deg_gene_col <- grep("^gene$|^Gene$", colnames(deg_all),
                     ignore.case = TRUE, value = TRUE)[1]
deg_lfc_col  <- grep("log2FC|log2FoldChange|avg_log2FC",
                     colnames(deg_all), ignore.case = TRUE, value = TRUE)[1]
deg_ct_col   <- grep("CellType|cell_type", colnames(deg_all),
                     ignore.case = TRUE, value = TRUE)[1]

message("\nDEG columns:")
message("  gene  = ", deg_gene_col)
message("  LFC   = ", deg_lfc_col)
message("  CT    = ", deg_ct_col)

# ── Build interleaved ct_geno factor and colors ───────────────────────────────
ct_levels <- levels(factor(obj_male$cell_type,
                           levels = atac_new_order[atac_new_order %in%
                                                     unique(obj_male$cell_type)]))

interleaved_levels <- as.vector(rbind(
  paste(ct_levels, "Ctrl", sep = " \u2014 "),
  paste(ct_levels, "KO",   sep = " \u2014 ")
))

obj_male$ct_geno <- factor(
  paste(obj_male$cell_type, obj_male$genotype, sep = " \u2014 "),
  levels = interleaved_levels[interleaved_levels %in% unique(obj_male$ct_geno)]
)

ct_geno_colors <- unlist(lapply(ct_levels, function(ct) {
  base_col <- atac_colors[ct]
  if (is.na(base_col)) base_col <- "grey60"
  setNames(
    c(scales::alpha(base_col, 0.45), base_col),
    c(paste(ct, "Ctrl", sep = " \u2014 "),
      paste(ct, "KO",   sep = " \u2014 "))
  )
}))
ct_geno_colors <- ct_geno_colors[names(ct_geno_colors) %in% levels(obj_male$ct_geno)]

# ── All target genes from significantly enriched TF × CT pairs ───────────────
# Strategy:
#   1. Keep directional rows (Up / Down) if FDR < 0.05 — these are the most
#      informative: the TF's targets are enriched specifically among up- or
#      down-regulated DEGs, indicating coherent activation or repression.
#   2. For TF × CT pairs where only "All" passes FDR < 0.05 (no directional
#      row significant), fall back to "All" to avoid losing the pair entirely.
#   3. After gene-level expansion, deduplicate: if a gene appears in both Up
#      and Down rows for the same TF × CT pair, keep the row with the lower
#      FDR (stronger signal).

sig_directional <- tf_deg |>
  dplyr::filter(Direction %in% c("Up", "Down"), Perm_FDR < 0.05)

sig_all_only <- tf_deg |>
  dplyr::filter(Direction == "All", Perm_FDR < 0.05) |>
  dplyr::anti_join(sig_directional,
                   by = c("TF" = "TF", "Cell_Type" = "Cell_Type"))

top_targets <- dplyr::bind_rows(sig_directional, sig_all_only) |>
  dplyr::mutate(Overlap_Genes = as.character(Overlap_Genes)) |>
  tidyr::separate_rows(Overlap_Genes, sep = ",\\s*") |>
  dplyr::rename(Gene = Overlap_Genes, CellType = Cell_Type, fdr = Perm_FDR) |>
  dplyr::filter(!is.na(Gene), Gene != "", !grepl("^ENSMUSG", Gene)) |>
  # If a gene appears in both Up and Down for the same TF × CT, keep lower FDR
  dplyr::group_by(TF, CellType, Gene) |>
  dplyr::slice_min(fdr, n = 1, with_ties = FALSE) |>
  dplyr::ungroup() |>
  dplyr::left_join(
    deg_all |>
      dplyr::rename(Gene = all_of(deg_gene_col),
                    LFC  = all_of(deg_lfc_col),
                    CT   = all_of(deg_ct_col)) |>
      dplyr::select(Gene, CT, LFC) |>
      dplyr::filter(!is.na(LFC)),
    by = c("Gene" = "Gene", "CellType" = "CT")
  ) |>
  dplyr::mutate(direction = dplyr::case_when(
    LFC >  0 ~ "Up in KO",
    LFC <  0 ~ "Down in KO",
    TRUE     ~ NA_character_
  ))

message("Total plots to generate: ", nrow(top_targets),
        " across ", length(unique(paste(top_targets$TF, top_targets$CellType))),
        " TF \u00d7 CT pairs")
print(top_targets |> dplyr::count(TF, CellType) |> dplyr::arrange(dplyr::desc(n)))

# just verify:
message("Columns in top_targets: ", paste(colnames(top_targets), collapse = " | "))


########## Difference -----
# ══════════════════════════════════════════════════════════════════════════════
# STEP 8D: DIFFERENTIAL COVERAGE — KO minus Ctrl
# ══════════════════════════════════════════════════════════════════════════════

library(ggtext)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(GenomicFeatures)

# ── Download and load ENCODE mm10 cCREs (run once) ───────────────────────────
ccre_bed <- "mm10_cCREs.bed"

if (!file.exists(ccre_bed)) {
  message("Downloading ENCODE mm10 cCREs (V3)...")
  
  # Current working URL via ENCODE portal file download
  download.file(
    url      = "https://www.encodeproject.org/files/ENCFF904ZZH/@@download/ENCFF904ZZH.bed.gz",
    destfile = "mm10_cCREs.bed.gz",
    mode     = "wb"
  )
  R.utils::gunzip("mm10_cCREs.bed.gz", destname = "mm10_cCREs.bed", overwrite = TRUE)
  
  # Verify
  lines <- readLines("mm10_cCREs.bed", n = 3)
  cat(lines, sep = "\n")
}

ccre_raw <- read.table("mm10_cCREs.bed", sep = "\t", header = FALSE,
                       col.names = c("chr", "start", "end", "name",
                                     "score", "strand", "thickStart", "thickEnd",
                                     "itemRgb", "type", "agnostic"),
                       fill = TRUE)

ccre_gr <- GenomicRanges::GRanges(
  seqnames = ccre_raw$chr,
  ranges   = IRanges::IRanges(ccre_raw$start + 1L, ccre_raw$end),
  type     = ccre_raw$type
)

message("cCREs loaded: ", length(ccre_gr), " elements")
message("Types: ", paste(sort(unique(ccre_raw$type)), collapse = ", "))

#Load data and set up folders

DefaultAssay(obj_male) <- "peaks"

dir.create("coverage_plots_targets/differential", showWarnings = FALSE, recursive = TRUE)

diff_colors <- c("UP" = "#b2182b", "DOWN" = "#2166ac")

# ── Helper: resolve gene alias to name present in annotation ─────────────────
resolve_gene_alias <- function(gene, annot) {
  # 1. Direct hit in annotation
  if (any(!is.na(annot$gene_name) & annot$gene_name == gene))
    return(gene)
  
  # 2. Existing orgDb alias lookup (your current code)
  alias_result <- tryCatch({
    AnnotationDbi::select(org.Mm.eg.db,
                          keys    = gene,
                          columns = c("SYMBOL", "ALIAS"),
                          keytype = "ALIAS") |>
      dplyr::pull(SYMBOL) |> unique() |> na.omit()
  }, error = function(e) character(0))
  
  for (sym in alias_result) {
    if (any(!is.na(annot$gene_name) & annot$gene_name == sym))
      return(sym)
  }
  
  # 3. Human → mouse homologue via biomaRt
  mm_sym <- tryCatch({
    human  <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    mouse  <- biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl")
    homol  <- biomaRt::getLDS(
      attributes  = "hgnc_symbol",
      filters     = "hgnc_symbol",
      values      = gene,
      mart        = human,
      attributesL = "mgi_symbol",
      martL       = mouse,
      uniqueRows  = TRUE
    )
    unique(na.omit(homol[[2]]))
  }, error = function(e) character(0))
  
  for (sym in mm_sym) {
    if (any(!is.na(annot$gene_name) & annot$gene_name == sym))
      return(sym)
  }
  
  # 4. Nothing found — return original, TxDb fallback will handle coordinates
  return(gene)
}

# ── Helper: resolve gene to GRanges via TxDb ─────────────────────────────────
gene_to_region <- function(gene) {
  tryCatch({
    txdb   <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
    entrez <- AnnotationDbi::select(org.Mm.eg.db, keys = gene,
                                    keytype = "ALIAS", columns = "ENTREZID")$ENTREZID
    entrez <- unique(entrez[!is.na(entrez)])
    if (length(entrez) == 0) return(NULL)
    gr <- GenomicFeatures::genes(txdb, filter = list(gene_id = entrez[1]))
    if (length(gr) == 0) return(NULL)
    gr
  }, error = function(e) NULL)
}

# ── Helper: extract x,y from a single-group CoveragePlot ─────────────────────
get_cov_df <- function(obj_sub, gene, label, extend_up, extend_dn) {
  p <- tryCatch(
    CoveragePlot(
      object            = obj_sub,
      region            = gene,
      assay             = "peaks",
      group.by          = "genotype",
      annotation        = FALSE,
      peaks             = FALSE,
      links             = FALSE,
      extend.upstream   = extend_up,
      extend.downstream = extend_dn
    ),
    error = function(e) {
      # Fallback: resolve to genomic coordinates via TxDb
      gr <- gene_to_region(gene)
      if (is.null(gr)) return(NULL)
      region_str <- paste0(
        as.character(seqnames(gr)), "-",
        start(gr) - extend_up, "-",
        end(gr)   + extend_dn
      )
      tryCatch(
        CoveragePlot(
          object            = obj_sub,
          region            = region_str,
          assay             = "peaks",
          group.by          = "genotype",
          annotation        = FALSE,
          peaks             = FALSE,
          links             = FALSE,
          extend.upstream   = 0,
          extend.downstream = 0
        ),
        error = function(e2) NULL
      )
    }
  )
  if (is.null(p)) return(NULL)
  built <- tryCatch(ggplot_build(p), error = function(e) NULL)
  if (is.null(built)) return(NULL)
  df <- tryCatch({
    d <- built$data[[1]]
    if (nrow(d) > 0 && all(c("x", "y") %in% colnames(d))) d else NULL
  }, error = function(e) NULL)
  if (is.null(df)) return(NULL)
  data.frame(x = df$x, y = df$y, genotype = label, stringsAsFactors = FALSE)
}

# ── Helper: clean NA coords from annotation ───────────────────────────────────
clean_annotation <- function(obj) {
  annot <- Annotation(obj)
  if (!is.null(annot) && length(annot) > 0) {
    keep <- !is.na(start(annot)) & !is.na(end(annot))
    Annotation(obj) <- annot[keep]
  }
  obj
}

# ── Helper: gene track with composite exon model + cCRE overlay ──────────────
make_gene_track <- function(obj_ct, gene, extend_up, extend_dn,
                            ccre_gr       = NULL,
                            row_height_in = 0.35, label_size = 3) {
  
  annot  <- Annotation(obj_ct)
  target <- if (!is.null(annot) && length(annot) > 0)
    annot[!is.na(annot$gene_name) & annot$gene_name == gene]
  else GenomicRanges::GRanges()
  
  # ── TxDb fallback when gene absent from Signac annotation ──────────────────
  if (length(target) == 0) {
    gr_fb <- gene_to_region(gene)
    if (is.null(gr_fb) || length(gr_fb) == 0) return(list(plot = NULL, height_in = 0))
    chr  <- as.character(seqnames(gr_fb[1]))
    tss  <- ifelse(as.character(strand(gr_fb[1])) == "+",
                   start(gr_fb[1]), end(gr_fb[1]))
    xmin <- tss - extend_up
    xmax <- tss + extend_dn
  } else {
    chr  <- as.character(seqnames(target[1]))
    tss  <- ifelse(as.character(strand(target[1])) == "+",
                   min(start(target)), max(end(target)))
    xmin <- tss - extend_up
    xmax <- tss + extend_dn
  }
  
  # ── From here: unchanged — uses xmin/xmax/chr regardless of how derived ────
  if (is.null(annot) || length(annot) == 0) {
    in_reg <- GenomicRanges::GRanges()
  } else {
    in_reg <- annot[
      !is.na(start(annot)) & !is.na(end(annot)) &
        as.character(seqnames(annot)) == chr &
        start(annot) <= xmax & end(annot) >= xmin
    ]
  }
  
  # If nothing in the region at all (not even neighbours), still render cCREs
  if (length(in_reg) == 0) {
    df          <- data.frame()
    transcripts <- data.frame(gene_name = character(), start_cl = numeric(),
                              end_cl = numeric(), strand = character(),
                              row = integer(), is_target = logical(),
                              stringsAsFactors = FALSE)
    exons       <- data.frame(gene_name = character(), start_cl = numeric(),
                              end_cl = numeric(), row = integer(),
                              fill_id = character(), is_target = logical(),
                              stringsAsFactors = FALSE)
    arrows_df   <- data.frame(gene_name = character(), row = integer(),
                              strand = character(), tick_pos = numeric(),
                              x_from = numeric(), x_to = numeric(),
                              is_target = logical(), stringsAsFactors = FALSE)
    n_rows      <- 1L
  } else {
    df <- as.data.frame(in_reg, row.names = NULL) |>
      dplyr::mutate(
        gene_name = ifelse(is.na(gene_name), "unknown", as.character(gene_name)),
        strand    = as.character(strand),
        start_cl  = pmax(start, xmin),
        end_cl    = pmin(end,   xmax)
      ) |>
      dplyr::filter(start_cl < end_cl)
    
    if (nrow(df) == 0) {
      transcripts <- data.frame(gene_name = character(), start_cl = numeric(),
                                end_cl = numeric(), strand = character(),
                                row = integer(), is_target = logical(),
                                stringsAsFactors = FALSE)
      exons <- arrows_df <- transcripts
      n_rows <- 1L
    } else {
      transcripts <- df |>
        dplyr::group_by(gene_name) |>
        dplyr::summarise(start_cl = min(start_cl), end_cl = max(end_cl),
                         strand = dplyr::first(strand), .groups = "drop") |>
        dplyr::arrange(start_cl)
      
      exon_raw <- df |>
        dplyr::filter(type %in% c("exon", "cds")) |>
        dplyr::select(gene_name, start_cl, end_cl) |>
        dplyr::arrange(gene_name, start_cl)
      
      exons <- exon_raw |>
        dplyr::group_by(gene_name) |>
        dplyr::group_modify(~ {
          d <- dplyr::arrange(.x, start_cl)
          merged <- d[1, , drop = FALSE]
          if (nrow(d) > 1) {
            for (j in 2:nrow(d)) {
              if (d$start_cl[j] <= merged$end_cl[nrow(merged)]) {
                merged$end_cl[nrow(merged)] <- max(merged$end_cl[nrow(merged)], d$end_cl[j])
              } else {
                merged <- dplyr::bind_rows(merged, d[j, , drop = FALSE])
              }
            }
          }
          merged
        }) |>
        dplyr::ungroup()
      
      # Greedy row packer
      transcripts$row <- 1L
      row_ends <- c(transcripts$end_cl[1])
      if (nrow(transcripts) > 1) {
        for (k in 2:nrow(transcripts)) {
          placed <- FALSE
          for (r in seq_along(row_ends)) {
            if (transcripts$start_cl[k] > row_ends[r] + 500) {
              transcripts$row[k] <- as.integer(r)
              row_ends[r] <- transcripts$end_cl[k]
              placed <- TRUE; break
            }
          }
          if (!placed) {
            row_ends <- c(row_ends, transcripts$end_cl[k])
            transcripts$row[k] <- as.integer(length(row_ends))
          }
        }
      }
      
      n_rows <- max(transcripts$row)
      
      row_lookup <- transcripts |>
        dplyr::select(gene_name, row) |>
        dplyr::distinct(gene_name, .keep_all = TRUE)
      
      exons <- exons |>
        dplyr::left_join(row_lookup, by = "gene_name") |>
        dplyr::filter(!is.na(row))
      
      # Arrow ticks
      arrow_spacing <- (xmax - xmin) / 15
      
      arrows_df <- do.call(dplyr::bind_rows, lapply(seq_len(nrow(transcripts)), function(k) {
        tx  <- transcripts[k, ]
        ex  <- exons[exons$gene_name == tx$gene_name, , drop = FALSE]
        from_pos <- tx$start_cl + arrow_spacing / 2
        to_pos   <- tx$end_cl   - arrow_spacing / 2
        
        pos <- if (to_pos >= from_pos) {
          seq(from_pos, to_pos, by = arrow_spacing)
        } else {
          # transcript narrower than one arrow spacing — place single arrow at midpoint
          (tx$start_cl + tx$end_cl) / 2
        }
        
        if (nrow(ex) > 0) {
          exon_buf <- arrow_spacing * 0.5
          in_exon  <- sapply(pos, function(p)
            any(p >= (ex$start_cl - exon_buf) & p <= (ex$end_cl + exon_buf)))
          pos_intron <- pos[!in_exon]
          
          if (length(pos_intron) == 0) {
            ex_sorted  <- ex[order(ex$start_cl), ]
            if (nrow(ex_sorted) >= 2) {
              gap_starts <- ex_sorted$end_cl[-nrow(ex_sorted)]
              gap_ends   <- ex_sorted$start_cl[-1]
              gap_widths <- gap_ends - gap_starts
              best       <- which.max(gap_widths)
              if (gap_widths[best] > (xmax - xmin) * 0.01) {
                pos_intron <- (gap_starts[best] + gap_ends[best]) / 2
              } else {
                return(data.frame(gene_name=character(),row=integer(),strand=character(),
                                  tick_pos=numeric(),x_from=numeric(),x_to=numeric(),
                                  stringsAsFactors=FALSE))
              }
            } else {
              return(data.frame(gene_name=character(),row=integer(),strand=character(),
                                tick_pos=numeric(),x_from=numeric(),x_to=numeric(),
                                stringsAsFactors=FALSE))
            }
          }
          pos <- pos_intron
        }
        
        data.frame(gene_name=tx$gene_name, row=tx$row, strand=tx$strand,
                   tick_pos=pos,
                   x_from=ifelse(tx$strand=="+", pos-arrow_spacing*0.35, pos+arrow_spacing*0.35),
                   x_to  =ifelse(tx$strand=="+", pos+arrow_spacing*0.35, pos-arrow_spacing*0.35),
                   stringsAsFactors=FALSE)
      }))
      
      transcripts$is_target <- transcripts$gene_name == gene
      arrows_df$is_target   <- if (nrow(arrows_df) > 0) arrows_df$gene_name == gene else logical(0)
      exons$is_target       <- exons$gene_name == gene
    }
  }
  
  all_genes   <- unique(c(transcripts$gene_name, gene))
  gene_colors <- setNames(ifelse(all_genes == gene, "#1a6e1a", "#555555"), all_genes)
  
  exons$fill_id <- exons$gene_name
  
  # ── cCRE data ─────────────────────────────────────────────────────────────
  ccre_df   <- NULL
  ccre_ymin <- -0.55
  ccre_ymax <- -0.15
  
  ccre_colors <- c("PLS"="#FF0000","pELS"="#FFA500","dELS"="#FFCD00",
                   "CTCF-only"="#00B0F0","DNase-H3K4me3"="#00B050")
  
  if (!is.null(ccre_gr) && length(ccre_gr) > 0) {
    hits <- ccre_gr[as.character(seqnames(ccre_gr)) == chr &
                      start(ccre_gr) <= xmax & end(ccre_gr) >= xmin]
    if (length(hits) > 0) {
      ccre_df <- data.frame(
        start_cl = pmax(start(hits), xmin),
        end_cl   = pmin(end(hits),   xmax),
        type     = as.character(hits$type),
        stringsAsFactors = FALSE
      ) |>
        dplyr::filter(start_cl < end_cl) |>
        dplyr::mutate(type = dplyr::case_when(
          grepl("PLS",           type) ~ "PLS",
          grepl("pELS",          type) ~ "pELS",
          grepl("dELS",          type) ~ "dELS",
          grepl("CTCF",          type) ~ "CTCF-only",
          grepl("DNase-H3K4me3", type) ~ "DNase-H3K4me3",
          TRUE ~ "other")) |>
        dplyr::filter(type != "other")
    }
  }
  
  has_ccre   <- !is.null(ccre_df) && nrow(ccre_df) > 0
  ccre_extra <- if (has_ccre) 0.55 else 0
  h_in       <- max(0.8, n_rows * row_height_in + 0.7 + ccre_extra)
  y_floor    <- if (has_ccre) ccre_ymin - 0.15 else 0.3
  
  if (has_ccre) ccre_df$fill_id <- ccre_df$type
  unified_fills <- c(gene_colors, ccre_colors)
  
  # ── Build plot ─────────────────────────────────────────────────────────────
  p <- ggplot()
  
  if (nrow(transcripts) > 0) {
    p <- p +
      geom_segment(data = transcripts,
                   aes(x=start_cl, xend=end_cl, y=row, yend=row, color=gene_name),
                   linewidth=0.5)
  }
  
  if (nrow(exons) > 0) {
    p <- p +
      geom_rect(data = exons,
                aes(xmin=start_cl, xmax=end_cl, ymin=row-0.18, ymax=row+0.18, fill=fill_id),
                color=NA)
  }
  
  if (nrow(arrows_df) > 0) {
    p <- p +
      geom_segment(data = arrows_df,
                   aes(x=x_from, xend=x_to, y=row, yend=row, color=gene_name),
                   linewidth=0.3,
                   arrow=arrow(length=unit(0.14,"cm"), type="open", ends="last"))
  }
  
  if (nrow(transcripts) > 0) {
    target_tx  <- dplyr::filter(transcripts,  is_target)
    other_tx   <- dplyr::filter(transcripts, !is_target)
    if (nrow(target_tx) > 0)
      p <- p + geom_text(data=target_tx,
                         aes(x=(start_cl+end_cl)/2, y=row+0.38, label=gene_name),
                         color="#1a6e1a", size=label_size,
                         hjust=0.5, vjust=0, fontface="bold.italic")
    if (nrow(other_tx) > 0)
      p <- p + geom_text(data=other_tx,
                         aes(x=(start_cl+end_cl)/2, y=row+0.38, label=gene_name),
                         color="#888888", size=label_size,
                         hjust=0.5, vjust=0, fontface="italic")
  } else {
    # Gene absent from annotation: draw a label at TSS position
    p <- p + annotate("text", x=(xmin+xmax)/2, y=1, label=gene,
                      color="#1a6e1a", size=label_size,
                      hjust=0.5, vjust=0, fontface="bold.italic")
  }
  
  p <- p + scale_color_manual(values=gene_colors, guide="none")
  
  if (has_ccre) {
    p <- p +
      geom_hline(yintercept=0, color="grey80", linewidth=0.3, linetype="dashed") +
      geom_rect(data=ccre_df,
                aes(xmin=start_cl, xmax=end_cl, ymin=ccre_ymin, ymax=ccre_ymax, fill=fill_id),
                color=NA, inherit.aes=FALSE) +
      geom_text(data=data.frame(x=xmin, y=(ccre_ymin+ccre_ymax)/2, label="cCREs"),
                aes(x=x, y=y, label=label),
                hjust=0, size=2.8, color="grey40", inherit.aes=FALSE) +
      scale_fill_manual(
        values=unified_fills, breaks=names(ccre_colors),
        labels=c("PLS"="Promoter (PLS)","pELS"="Proximal enhancer (pELS)",
                 "dELS"="Distal enhancer (dELS)","CTCF-only"="CTCF insulator",
                 "DNase-H3K4me3"="DNase+H3K4me3"),
        name="ENCODE cCRE",
        guide=guide_legend(override.aes=list(size=3),
                           keyheight=unit(0.35,"cm"), keywidth=unit(0.35,"cm")))
  } else {
    p <- p + scale_fill_manual(values=unified_fills, guide="none")
  }
  
  p <- p +
    scale_x_continuous(limits=c(xmin,xmax), expand=expansion(mult=c(0,0)),
                       labels=scales::label_number(scale=1e-6,suffix=" Mb",accuracy=0.01)) +
    scale_y_continuous(limits=c(y_floor, n_rows+0.9), expand=c(0,0)) +
    labs(x=NULL, y="Genes") +
    theme_minimal(base_size=10) +
    theme(panel.grid=element_blank(),
          axis.line.x=element_line(color="grey70"),
          axis.ticks.x=element_line(color="grey70"),
          axis.text.x=element_text(size=8),
          axis.title.y=element_text(size=9, angle=90),
          axis.text.y=element_blank(), axis.ticks.y=element_blank(),
          legend.position=if(has_ccre)"right" else "none",
          legend.key.size=unit(0.35,"cm"),
          legend.text=element_text(size=7), legend.title=element_text(size=8),
          plot.margin=margin(2,2,6,2))
  
  list(plot=p, height_in=h_in)
}

# ── Plot loop ─────────────────────────────────────────────────────────────────
failed  <- c()
success <- c()

for (i in seq_len(nrow(top_targets))) {
  
  row     <- top_targets[i, ]
  gene    <- row$Gene
  ct      <- row$CellType
  tf      <- row$TF
  dir_lbl <- row$direction
  ct_safe <- gsub("[^A-Za-z0-9]", "_", ct)
  tf_safe <- gsub("[^A-Za-z0-9]", "_", tf)
  
  outfile <- file.path("coverage_plots_targets/differential",
                       paste0(ct_safe, "__", tf_safe,"__", gene, ".pdf"))
  
  if (file.exists(outfile)) {
    message("  Skip (exists): ", gene, " [", tf, " / ", ct, "]")
    success <- c(success, gene)
    next
  }
  
  obj_ct <- tryCatch(
    subset(obj_male, subset = cell_type == ct),
    error = function(e) NULL
  )
  if (is.null(obj_ct) || ncol(obj_ct) < 20) {
    message("  Skip [", gene, " / ", ct, "] — too few cells")
    next
  }
  
  obj_ct$genotype <- factor(obj_ct$genotype, levels = c("Ctrl", "KO"))
  obj_ct <- clean_annotation(obj_ct)
  
  # ── Resolve gene alias to name present in annotation ───────────────────────
  gene_resolved <- resolve_gene_alias(gene, Annotation(obj_ct))
  if (gene_resolved != gene)
    message("  Alias resolved: ", gene, " \u2192 ", gene_resolved)
  
  # ── Skip if not in annotation AND not resolvable via TxDb ──────────────────
  in_annot <- sum(!is.na(Annotation(obj_ct)$gene_name) &
                    Annotation(obj_ct)$gene_name == gene_resolved, na.rm = TRUE)
  if (in_annot == 0) {
    gr_check <- gene_to_region(gene_resolved)
    if (is.null(gr_check)) {
      message("  Skip [", gene, "]: not in annotation or TxDb — likely pseudogene/lncRNA")
      failed <- c(failed, gene)
      next
    }
  }
  
  tryCatch({
    
    obj_ctrl <- clean_annotation(subset(obj_ct, subset = genotype == "Ctrl"))
    obj_ko   <- clean_annotation(subset(obj_ct, subset = genotype == "KO"))
    obj_ctrl$genotype <- factor("Ctrl", levels = "Ctrl")
    obj_ko$genotype   <- factor("KO",   levels = "KO")
    
    df_ctrl <- get_cov_df(obj_ctrl, gene_resolved, "Ctrl", 10000, 10000)
    df_ko   <- get_cov_df(obj_ko,   gene_resolved, "KO",   10000, 10000)
    
    if (is.null(df_ctrl) && is.null(df_ko)) stop("Could not extract coverage data")
    
    n_ctrl <- ncol(obj_ctrl)
    n_ko   <- ncol(obj_ko)
    
    ref_x    <- if (!is.null(df_ko)) df_ko$x else df_ctrl$x
    x_range  <- range(ref_x)
    n_bins   <- length(ref_x)
    common_x <- seq(x_range[1], x_range[2], length.out = n_bins)
    
    y_ctrl <- if (!is.null(df_ctrl)) {
      approx(df_ctrl$x, df_ctrl$y / n_ctrl, xout = common_x, rule = 2)$y
    } else {
      rep(0, n_bins)
    }
    
    y_ko <- if (!is.null(df_ko)) {
      approx(df_ko$x, df_ko$y / n_ko, xout = common_x, rule = 2)$y
    } else {
      rep(0, n_bins)
    }
    
    # y_ctrl <- stats::filter(y_ctrl, rep(1/3, 3), circular = TRUE)
    # y_ko   <- stats::filter(y_ko,   rep(1/3, 3), circular = TRUE)
    
    # ── Gaussian smoother (replaces the 3-point boxcar) ──────────────────────
    gaussian_smooth <- function(y, sigma = 3) {
      hw  <- ceiling(3 * sigma)
      ker <- dnorm(-hw:hw, sd = sigma)
      ker <- ker / sum(ker)
      as.numeric(stats::filter(y, ker, circular = TRUE))
    }
    
    y_ctrl <- gaussian_smooth(as.numeric(y_ctrl), sigma = 8)
    y_ko   <- gaussian_smooth(as.numeric(y_ko),   sigma = 8)
    
    diff_df <- data.frame(
      x    = common_x,
      diff = as.numeric(y_ko) - as.numeric(y_ctrl)
    )
    
    split_segments <- function(x, y) {
      segs  <- list()
      cur_x <- c()
      cur_y <- c()
      for (k in seq_len(length(x))) {
        if (k > 1) {
          y0 <- y[k - 1]; y1 <- y[k]
          if (!is.na(y0) && !is.na(y1) && sign(y0) != sign(y1) && (y1 - y0) != 0) {
            x_zero <- x[k - 1] + (x[k] - x[k - 1]) * (-y0 / (y1 - y0))
            cur_x  <- c(cur_x, x_zero)
            cur_y  <- c(cur_y, 0)
            segs[[length(segs) + 1]] <- data.frame(x = cur_x, y = cur_y)
            cur_x <- c(x_zero)
            cur_y <- c(0)
          }
        }
        cur_x <- c(cur_x, x[k])
        cur_y <- c(cur_y, y[k])
      }
      if (length(cur_x) > 0) segs[[length(segs) + 1]] <- data.frame(x = cur_x, y = cur_y)
      segs
    }
    
    segments   <- split_segments(diff_df$x, diff_df$diff)
    seg_pos_df <- dplyr::bind_rows(lapply(segments, function(s) if (mean(s$y, na.rm = TRUE) >= 0) s else NULL))
    seg_neg_df <- dplyr::bind_rows(lapply(segments, function(s) if (mean(s$y, na.rm = TRUE) <  0) s else NULL))
    
    yabs <- max(abs(diff_df$diff), na.rm = TRUE) * 1.05
    if (!is.finite(yabs) || yabs == 0) yabs <- 1
    
    fill_col <- dplyr::case_when(
      grepl("up",   dir_lbl, ignore.case = TRUE) ~ diff_colors["UP"],
      grepl("down", dir_lbl, ignore.case = TRUE) ~ diff_colors["DOWN"],
      TRUE ~ "#555555"
    )
    
    dir_color_hex <- dplyr::case_when(
      grepl("up",   dir_lbl, ignore.case = TRUE) ~ "#b2182b",
      grepl("down", dir_lbl, ignore.case = TRUE) ~ "#2166ac",
      TRUE ~ "#555555"
    )
    
    p_cov <- ggplot(diff_df, aes(x = x, y = diff)) +
      geom_hline(yintercept = 0, color = "grey50", linewidth = 0.4) +
      geom_area(data = seg_pos_df, aes(x = x, y = y),
                fill = diff_colors["UP"],   alpha = 0.75, color = NA) +
      geom_area(data = seg_neg_df, aes(x = x, y = y),
                fill = diff_colors["DOWN"], alpha = 0.75, color = NA) +
      scale_x_continuous(expand = expansion(mult = c(0, 0))) +
      scale_y_continuous(
        limits = c(-yabs, yabs),
        expand = expansion(mult = c(0, 0)),
        labels = function(x) as.character(abs(x))
      ) +
      labs(
        y     = "\u0394 Coverage per cell (KO \u2212 Ctrl)",
        x     = NULL,
        title = paste0(
          "**", gene, "**  |  Regulated by: ", tf,
          "  |  <span style='color:", dir_color_hex, "'>**Expression: ",
          dir_lbl, "**</span>"
        ),
        subtitle = paste0(
          "Cell type: ", ct,
          "  |  Differential ATAC (KO \u2212 Ctrl, male)  |  \u00b110kb",
          "  |  <span style='color:#b2182b'>**\u25a0 KO > Ctrl**</span>",
          "  <span style='color:#2166ac'>**\u25a0 KO < Ctrl**</span>"
        )
      ) +
      theme_minimal(base_size = 11) +
      theme(
        plot.title      = ggtext::element_markdown(size = 12),
        plot.subtitle   = ggtext::element_markdown(size = 9, color = "grey40"),
        panel.grid      = element_blank(),
        axis.line.y     = element_line(color = "grey70"),
        axis.ticks.y    = element_line(color = "grey70"),
        axis.text.x     = element_blank(),
        axis.ticks.x    = element_blank(),
        legend.position = "none",
        plot.margin     = margin(2, 2, 2, 2)
      )
    
    # ── Peaks panel ──────────────────────────────────────────────────────────
    p_peaks <- tryCatch(
      PeakPlot(
        object            = obj_ct,
        region            = gene_resolved,
        assay             = "peaks",
        extend.upstream   = 10000,
        extend.downstream = 10000
      ) + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()),
      error = function(e) NULL
    )
    
    # ── Gene track ───────────────────────────────────────────────────────────
    gene_track <- make_gene_track(obj_ct, gene_resolved, 10000, 10000,
                                  ccre_gr       = ccre_gr,
                                  row_height_in = 0.35,
                                  label_size    = 3)
    p_annot <- gene_track$plot
    annot_h <- gene_track$height_in
    
    peaks_h <- 0.5
    cov_h   <- 4.0
    title_h <- 0.8
    
    if (!is.null(p_peaks) && !is.null(p_annot)) {
      p_final <- p_cov / p_peaks / p_annot +
        plot_layout(heights = unit(c(cov_h, peaks_h, annot_h), "in"))
      total_h <- title_h + cov_h + peaks_h + annot_h + 0.5
    } else if (!is.null(p_annot)) {
      p_final <- p_cov / p_annot +
        plot_layout(heights = unit(c(cov_h, annot_h), "in"))
      total_h <- title_h + cov_h + annot_h + 0.5
    } else if (!is.null(p_peaks)) {
      p_final <- p_cov / p_peaks +
        plot_layout(heights = unit(c(cov_h, peaks_h), "in"))
      total_h <- title_h + cov_h + peaks_h + 0.5
    } else {
      p_final <- p_cov
      total_h <- title_h + cov_h + 0.5
    }
    
    ggsave(outfile, p_final, width = 10, height = total_h, limitsize = FALSE)
    success <- c(success, gene)
    message("  Saved: ", gene, " [", tf, " \u2192 ", ct, ", ", dir_lbl, "]")
    
  }, error = function(e) {
    failed <<- c(failed, gene)
    message("  FAILED [", gene, "]: ", e$message)
  })
}

message("\n\u2500\u2500 Step 8D summary \u2500\u2500")
message("  Saved:  ", length(success))
message("  Failed: ", length(failed), " \u2014 ",
        if (length(failed) > 0) paste(failed, collapse = ", ") else "none")
message("\n=== Step 8D complete. ===")


######### Load differential data for co-regulators ------
# ── Load TF-DEG overlap (multi-eRegulon) ─────────────────────────────────────
tf_deg <- read.csv(file.path(further_dir,
                             "TF_DEG_overlap/TF_DEG_overlap_permutation_KO_vs_Ctrl_multi_eRegulon_spec.csv"))
message("TF-DEG overlap loaded: ", nrow(tf_deg), " rows")
message("Columns: ", paste(colnames(tf_deg), collapse = " | "))

# ── Identify gene + cell type + TF columns ────────────────────────────────────
gene_col <- grep("^gene$|^Gene$|^target$|^Target$", colnames(tf_deg),
                 ignore.case = TRUE, value = TRUE)[1]
tf_col   <- grep("^TF$|^tf$", colnames(tf_deg),
                 ignore.case = TRUE, value = TRUE)[1]
ct_col   <- grep("CellType|cell_type|celltype", colnames(tf_deg),
                 ignore.case = TRUE, value = TRUE)[1]
pval_col <- grep("fdr|padj|adj", colnames(tf_deg),
                 ignore.case = TRUE, value = TRUE)[1]

message("\nMapped columns:")
message("  gene  = ", gene_col)
message("  TF    = ", tf_col)
message("  CT    = ", ct_col)
message("  pval  = ", pval_col)

# ── Load DEG direction (up/down in KO) ────────────────────────────────────────
deg_dir_path <- file.path(further_dir, "DEG_KO_vs_Ctrl")
deg_all      <- read.csv(file.path(deg_dir_path,
                                   "DEG_all_celltypes_KO_vs_Ctrl.csv"))
deg_gene_col <- grep("^gene$|^Gene$", colnames(deg_all),
                     ignore.case = TRUE, value = TRUE)[1]
deg_lfc_col  <- grep("log2FC|log2FoldChange|avg_log2FC",
                     colnames(deg_all), ignore.case = TRUE, value = TRUE)[1]
deg_ct_col   <- grep("CellType|cell_type", colnames(deg_all),
                     ignore.case = TRUE, value = TRUE)[1]

message("\nDEG columns:")
message("  gene  = ", deg_gene_col)
message("  LFC   = ", deg_lfc_col)
message("  CT    = ", deg_ct_col)

# ── Output directory ──────────────────────────────────────────────────────────
multi_diff_dir <- "/Users/cnbr/Downloads/Seurat_scATAC-seq/coverage_plots_targets/differential_multi"
dir.create(multi_diff_dir, showWarnings = FALSE, recursive = TRUE)
message("Output directory: ", multi_diff_dir)

# ── Build top_targets from multi-eRegulon CSV ─────────────────────────────────
sig_directional <- tf_deg |>
  dplyr::filter(Direction %in% c("Up", "Down"), Perm_FDR < 0.05)

sig_all_only <- tf_deg |>
  dplyr::filter(Direction == "All", Perm_FDR < 0.05) |>
  dplyr::anti_join(sig_directional,
                   by = c("TF" = "TF", "Cell_Type" = "Cell_Type"))

top_targets <- dplyr::bind_rows(sig_directional, sig_all_only) |>
  dplyr::mutate(Overlap_Genes = as.character(Overlap_Genes)) |>
  tidyr::separate_rows(Overlap_Genes, sep = ",\\s*") |>
  dplyr::rename(Gene = Overlap_Genes, CellType = Cell_Type, fdr = Perm_FDR) |>
  dplyr::filter(!is.na(Gene), Gene != "", !grepl("^ENSMUSG", Gene)) |>
  dplyr::group_by(TF, CellType, Gene) |>
  dplyr::slice_min(fdr, n = 1, with_ties = FALSE) |>
  dplyr::ungroup() |>
  dplyr::left_join(
    deg_all |>
      dplyr::rename(Gene = all_of(deg_gene_col),
                    LFC  = all_of(deg_lfc_col),
                    CT   = all_of(deg_ct_col)) |>
      dplyr::select(Gene, CT, LFC) |>
      dplyr::filter(!is.na(LFC)),
    by = c("Gene" = "Gene", "CellType" = "CT")
  ) |>
  dplyr::mutate(direction = dplyr::case_when(
    LFC >  0 ~ "Up in KO",
    LFC <  0 ~ "Down in KO",
    TRUE     ~ NA_character_
  ))

message("Total plots to generate: ", nrow(top_targets),
        " across ", length(unique(paste(top_targets$TF, top_targets$CellType))),
        " TF \u00d7 CT pairs")
print(top_targets |> dplyr::count(TF, CellType) |> dplyr::arrange(dplyr::desc(n)))
message("Columns in top_targets: ", paste(colnames(top_targets), collapse = " | "))


########## Difference multi -----
# ══════════════════════════════════════════════════════════════════════════════
# STEP 8D: DIFFERENTIAL COVERAGE — KO minus Ctrl
# ══════════════════════════════════════════════════════════════════════════════

library(ggtext)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(GenomicFeatures)

# ── Download and load ENCODE mm10 cCREs (run once) ───────────────────────────
ccre_bed <- "mm10_cCREs.bed"

if (!file.exists(ccre_bed)) {
  message("Downloading ENCODE mm10 cCREs (V3)...")
  
  # Current working URL via ENCODE portal file download
  download.file(
    url      = "https://www.encodeproject.org/files/ENCFF904ZZH/@@download/ENCFF904ZZH.bed.gz",
    destfile = "mm10_cCREs.bed.gz",
    mode     = "wb"
  )
  R.utils::gunzip("mm10_cCREs.bed.gz", destname = "mm10_cCREs.bed", overwrite = TRUE)
  
  # Verify
  lines <- readLines("mm10_cCREs.bed", n = 3)
  cat(lines, sep = "\n")
}

ccre_raw <- read.table("mm10_cCREs.bed", sep = "\t", header = FALSE,
                       col.names = c("chr", "start", "end", "name",
                                     "score", "strand", "thickStart", "thickEnd",
                                     "itemRgb", "type", "agnostic"),
                       fill = TRUE)

ccre_gr <- GenomicRanges::GRanges(
  seqnames = ccre_raw$chr,
  ranges   = IRanges::IRanges(ccre_raw$start + 1L, ccre_raw$end),
  type     = ccre_raw$type
)

message("cCREs loaded: ", length(ccre_gr), " elements")
message("Types: ", paste(sort(unique(ccre_raw$type)), collapse = ", "))

#Load data and set up folders

DefaultAssay(obj_male) <- "peaks"

dir.create("coverage_plots_targets/differential", showWarnings = FALSE, recursive = TRUE)

diff_colors <- c("UP" = "#b2182b", "DOWN" = "#2166ac")

# ── Helper: resolve gene alias to name present in annotation ─────────────────
resolve_gene_alias <- function(gene, annot) {
  # 1. Direct hit in annotation
  if (any(!is.na(annot$gene_name) & annot$gene_name == gene))
    return(gene)
  
  # 2. Existing orgDb alias lookup (your current code)
  alias_result <- tryCatch({
    AnnotationDbi::select(org.Mm.eg.db,
                          keys    = gene,
                          columns = c("SYMBOL", "ALIAS"),
                          keytype = "ALIAS") |>
      dplyr::pull(SYMBOL) |> unique() |> na.omit()
  }, error = function(e) character(0))
  
  for (sym in alias_result) {
    if (any(!is.na(annot$gene_name) & annot$gene_name == sym))
      return(sym)
  }
  
  # 3. Human → mouse homologue via biomaRt
  mm_sym <- tryCatch({
    human  <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    mouse  <- biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl")
    homol  <- biomaRt::getLDS(
      attributes  = "hgnc_symbol",
      filters     = "hgnc_symbol",
      values      = gene,
      mart        = human,
      attributesL = "mgi_symbol",
      martL       = mouse,
      uniqueRows  = TRUE
    )
    unique(na.omit(homol[[2]]))
  }, error = function(e) character(0))
  
  for (sym in mm_sym) {
    if (any(!is.na(annot$gene_name) & annot$gene_name == sym))
      return(sym)
  }
  
  # 4. Nothing found — return original, TxDb fallback will handle coordinates
  return(gene)
}

# ── Helper: resolve gene to GRanges via TxDb ─────────────────────────────────
gene_to_region <- function(gene) {
  tryCatch({
    txdb   <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
    entrez <- AnnotationDbi::select(org.Mm.eg.db, keys = gene,
                                    keytype = "ALIAS", columns = "ENTREZID")$ENTREZID
    entrez <- unique(entrez[!is.na(entrez)])
    if (length(entrez) == 0) return(NULL)
    gr <- GenomicFeatures::genes(txdb, filter = list(gene_id = entrez[1]))
    if (length(gr) == 0) return(NULL)
    gr
  }, error = function(e) NULL)
}

# ── Helper: extract x,y from a single-group CoveragePlot ─────────────────────
get_cov_df <- function(obj_sub, gene, label, extend_up, extend_dn) {
  p <- tryCatch(
    CoveragePlot(
      object            = obj_sub,
      region            = gene,
      assay             = "peaks",
      group.by          = "genotype",
      annotation        = FALSE,
      peaks             = FALSE,
      links             = FALSE,
      extend.upstream   = extend_up,
      extend.downstream = extend_dn
    ),
    error = function(e) {
      # Fallback: resolve to genomic coordinates via TxDb
      gr <- gene_to_region(gene)
      if (is.null(gr)) return(NULL)
      region_str <- paste0(
        as.character(seqnames(gr)), "-",
        start(gr) - extend_up, "-",
        end(gr)   + extend_dn
      )
      tryCatch(
        CoveragePlot(
          object            = obj_sub,
          region            = region_str,
          assay             = "peaks",
          group.by          = "genotype",
          annotation        = FALSE,
          peaks             = FALSE,
          links             = FALSE,
          extend.upstream   = 0,
          extend.downstream = 0
        ),
        error = function(e2) NULL
      )
    }
  )
  if (is.null(p)) return(NULL)
  built <- tryCatch(ggplot_build(p), error = function(e) NULL)
  if (is.null(built)) return(NULL)
  df <- tryCatch({
    d <- built$data[[1]]
    if (nrow(d) > 0 && all(c("x", "y") %in% colnames(d))) d else NULL
  }, error = function(e) NULL)
  if (is.null(df)) return(NULL)
  data.frame(x = df$x, y = df$y, genotype = label, stringsAsFactors = FALSE)
}

# ── Helper: clean NA coords from annotation ───────────────────────────────────
clean_annotation <- function(obj) {
  annot <- Annotation(obj)
  if (!is.null(annot) && length(annot) > 0) {
    keep <- !is.na(start(annot)) & !is.na(end(annot))
    Annotation(obj) <- annot[keep]
  }
  obj
}

# ── Helper: gene track with composite exon model + cCRE overlay ──────────────
make_gene_track <- function(obj_ct, gene, extend_up, extend_dn,
                            ccre_gr       = NULL,
                            row_height_in = 0.35, label_size = 3) {
  
  annot  <- Annotation(obj_ct)
  target <- if (!is.null(annot) && length(annot) > 0)
    annot[!is.na(annot$gene_name) & annot$gene_name == gene]
  else GenomicRanges::GRanges()
  
  # ── TxDb fallback when gene absent from Signac annotation ──────────────────
  if (length(target) == 0) {
    gr_fb <- gene_to_region(gene)
    if (is.null(gr_fb) || length(gr_fb) == 0) return(list(plot = NULL, height_in = 0))
    chr  <- as.character(seqnames(gr_fb[1]))
    tss  <- ifelse(as.character(strand(gr_fb[1])) == "+",
                   start(gr_fb[1]), end(gr_fb[1]))
    xmin <- tss - extend_up
    xmax <- tss + extend_dn
  } else {
    # chr  <- as.character(seqnames(target[1]))
    # tss  <- ifelse(as.character(strand(target[1])) == "+",
    #                min(start(target)), max(end(target)))
    # xmin <- tss - extend_up
    # xmax <- tss + extend_dn
  #}
    chr  <- as.character(seqnames(target)[1])
    # ── KEY FIX: window is centred on gene midpoint, not just TSS ──────────
    gene_start <- min(start(target))
    gene_end   <- max(end(target))
    gene_mid   <- (gene_start + gene_end) / 2
    xmin       <- gene_mid - extend_up
    xmax       <- gene_mid + extend_dn
  
  }
  
  # ── From here: unchanged — uses xmin/xmax/chr regardless of how derived ────
  if (is.null(annot) || length(annot) == 0) {
    in_reg <- GenomicRanges::GRanges()
  } else {
    in_reg <- annot[
      !is.na(start(annot)) & !is.na(end(annot)) &
        as.character(seqnames(annot)) == chr &
        start(annot) <= xmax & end(annot) >= xmin
    ]
  }
  
  # If nothing in the region at all (not even neighbours), still render cCREs
  if (length(in_reg) == 0) {
    df          <- data.frame()
    transcripts <- data.frame(gene_name = character(), start_cl = numeric(),
                              end_cl = numeric(), strand = character(),
                              row = integer(), is_target = logical(),
                              stringsAsFactors = FALSE)
    exons       <- data.frame(gene_name = character(), start_cl = numeric(),
                              end_cl = numeric(), row = integer(),
                              fill_id = character(), is_target = logical(),
                              stringsAsFactors = FALSE)
    arrows_df   <- data.frame(gene_name = character(), row = integer(),
                              strand = character(), tick_pos = numeric(),
                              x_from = numeric(), x_to = numeric(),
                              is_target = logical(), stringsAsFactors = FALSE)
    n_rows      <- 1L
  } else {
    df <- as.data.frame(in_reg, row.names = NULL) |>
      dplyr::mutate(
        gene_name = ifelse(is.na(gene_name), "unknown", as.character(gene_name)),
        strand    = as.character(strand),
        start_cl  = pmax(start, xmin),
        end_cl    = pmin(end,   xmax)
      ) |>
      dplyr::filter(start_cl < end_cl)
    
    if (nrow(df) == 0) {
      transcripts <- data.frame(gene_name = character(), start_cl = numeric(),
                                end_cl = numeric(), strand = character(),
                                row = integer(), is_target = logical(),
                                stringsAsFactors = FALSE)
      exons <- arrows_df <- transcripts
      n_rows <- 1L
    } else {
      transcripts <- df |>
        dplyr::group_by(gene_name) |>
        dplyr::summarise(start_cl = min(start_cl), end_cl = max(end_cl),
                         strand = dplyr::first(strand), .groups = "drop") |>
        dplyr::arrange(start_cl)
      
      exon_raw <- df |>
        dplyr::filter(type %in% c("exon", "cds")) |>
        dplyr::select(gene_name, start_cl, end_cl) |>
        dplyr::arrange(gene_name, start_cl)
      
      exons <- exon_raw |>
        dplyr::group_by(gene_name) |>
        dplyr::group_modify(~ {
          d <- dplyr::arrange(.x, start_cl)
          merged <- d[1, , drop = FALSE]
          if (nrow(d) > 1) {
            for (j in 2:nrow(d)) {
              if (d$start_cl[j] <= merged$end_cl[nrow(merged)]) {
                merged$end_cl[nrow(merged)] <- max(merged$end_cl[nrow(merged)], d$end_cl[j])
              } else {
                merged <- dplyr::bind_rows(merged, d[j, , drop = FALSE])
              }
            }
          }
          merged
        }) |>
        dplyr::ungroup()
      
      # Greedy row packer
      transcripts$row <- 1L
      row_ends <- c(transcripts$end_cl[1])
      if (nrow(transcripts) > 1) {
        for (k in 2:nrow(transcripts)) {
          placed <- FALSE
          for (r in seq_along(row_ends)) {
            if (transcripts$start_cl[k] > row_ends[r] + 500) {
              transcripts$row[k] <- as.integer(r)
              row_ends[r] <- transcripts$end_cl[k]
              placed <- TRUE; break
            }
          }
          if (!placed) {
            row_ends <- c(row_ends, transcripts$end_cl[k])
            transcripts$row[k] <- as.integer(length(row_ends))
          }
        }
      }
      
      n_rows <- max(transcripts$row)
      
      row_lookup <- transcripts |>
        dplyr::select(gene_name, row) |>
        dplyr::distinct(gene_name, .keep_all = TRUE)
      
      exons <- exons |>
        dplyr::left_join(row_lookup, by = "gene_name") |>
        dplyr::filter(!is.na(row))
      
      # Arrow ticks
      arrow_spacing <- (xmax - xmin) / 15
      
      arrows_df <- do.call(dplyr::bind_rows, lapply(seq_len(nrow(transcripts)), function(k) {
        tx  <- transcripts[k, ]
        ex  <- exons[exons$gene_name == tx$gene_name, , drop = FALSE]
        from_pos <- tx$start_cl + arrow_spacing / 2
        to_pos   <- tx$end_cl   - arrow_spacing / 2
        
        pos <- if (to_pos >= from_pos) {
          seq(from_pos, to_pos, by = arrow_spacing)
        } else {
          # transcript narrower than one arrow spacing — place single arrow at midpoint
          (tx$start_cl + tx$end_cl) / 2
        }
        
        if (nrow(ex) > 0) {
          exon_buf <- arrow_spacing * 0.5
          in_exon  <- sapply(pos, function(p)
            any(p >= (ex$start_cl - exon_buf) & p <= (ex$end_cl + exon_buf)))
          pos_intron <- pos[!in_exon]
          
          if (length(pos_intron) == 0) {
            ex_sorted  <- ex[order(ex$start_cl), ]
            if (nrow(ex_sorted) >= 2) {
              gap_starts <- ex_sorted$end_cl[-nrow(ex_sorted)]
              gap_ends   <- ex_sorted$start_cl[-1]
              gap_widths <- gap_ends - gap_starts
              best       <- which.max(gap_widths)
              if (gap_widths[best] > (xmax - xmin) * 0.01) {
                pos_intron <- (gap_starts[best] + gap_ends[best]) / 2
              } else {
                return(data.frame(gene_name=character(),row=integer(),strand=character(),
                                  tick_pos=numeric(),x_from=numeric(),x_to=numeric(),
                                  stringsAsFactors=FALSE))
              }
            } else {
              return(data.frame(gene_name=character(),row=integer(),strand=character(),
                                tick_pos=numeric(),x_from=numeric(),x_to=numeric(),
                                stringsAsFactors=FALSE))
            }
          }
          pos <- pos_intron
        }
        
        data.frame(gene_name=tx$gene_name, row=tx$row, strand=tx$strand,
                   tick_pos=pos,
                   x_from=ifelse(tx$strand=="+", pos-arrow_spacing*0.35, pos+arrow_spacing*0.35),
                   x_to  =ifelse(tx$strand=="+", pos+arrow_spacing*0.35, pos-arrow_spacing*0.35),
                   stringsAsFactors=FALSE)
      }))
      
      transcripts$is_target <- transcripts$gene_name == gene
      arrows_df$is_target   <- if (nrow(arrows_df) > 0) arrows_df$gene_name == gene else logical(0)
      exons$is_target       <- exons$gene_name == gene
    }
  }
  
  all_genes   <- unique(c(transcripts$gene_name, gene))
  gene_colors <- setNames(ifelse(all_genes == gene, "#1a6e1a", "#555555"), all_genes)
  
  exons$fill_id <- exons$gene_name
  
  # ── cCRE data ─────────────────────────────────────────────────────────────
  ccre_df   <- NULL
  ccre_ymin <- -0.55
  ccre_ymax <- -0.15
  
  ccre_colors <- c("PLS"="#FF0000","pELS"="#FFA500","dELS"="#FFCD00",
                   "CTCF-only"="#00B0F0","DNase-H3K4me3"="#00B050")
  
  if (!is.null(ccre_gr) && length(ccre_gr) > 0) {
    hits <- ccre_gr[as.character(seqnames(ccre_gr)) == chr &
                      start(ccre_gr) <= xmax & end(ccre_gr) >= xmin]
    if (length(hits) > 0) {
      ccre_df <- data.frame(
        start_cl = pmax(start(hits), xmin),
        end_cl   = pmin(end(hits),   xmax),
        type     = as.character(hits$type),
        stringsAsFactors = FALSE
      ) |>
        dplyr::filter(start_cl < end_cl) |>
        dplyr::mutate(type = dplyr::case_when(
          grepl("PLS",           type) ~ "PLS",
          grepl("pELS",          type) ~ "pELS",
          grepl("dELS",          type) ~ "dELS",
          grepl("CTCF",          type) ~ "CTCF-only",
          grepl("DNase-H3K4me3", type) ~ "DNase-H3K4me3",
          TRUE ~ "other")) |>
        dplyr::filter(type != "other")
    }
  }
  
  has_ccre   <- !is.null(ccre_df) && nrow(ccre_df) > 0
  ccre_extra <- if (has_ccre) 0.55 else 0
  h_in       <- max(0.8, n_rows * row_height_in + 0.7 + ccre_extra)
  y_floor    <- if (has_ccre) ccre_ymin - 0.15 else 0.3
  
  if (has_ccre) ccre_df$fill_id <- ccre_df$type
  unified_fills <- c(gene_colors, ccre_colors)
  
  # ── Build plot ─────────────────────────────────────────────────────────────
  p <- ggplot()
  
  if (nrow(transcripts) > 0) {
    p <- p +
      geom_segment(data = transcripts,
                   aes(x=start_cl, xend=end_cl, y=row, yend=row, color=gene_name),
                   linewidth=0.5)
  }
  
  if (nrow(exons) > 0) {
    p <- p +
      geom_rect(data = exons,
                aes(xmin=start_cl, xmax=end_cl, ymin=row-0.18, ymax=row+0.18, fill=fill_id),
                color=NA)
  }
  
  if (nrow(arrows_df) > 0) {
    p <- p +
      geom_segment(data = arrows_df,
                   aes(x=x_from, xend=x_to, y=row, yend=row, color=gene_name),
                   linewidth=0.3,
                   arrow=arrow(length=unit(0.14,"cm"), type="open", ends="last"))
  }
  
  if (nrow(transcripts) > 0) {
    target_tx  <- dplyr::filter(transcripts,  is_target)
    other_tx   <- dplyr::filter(transcripts, !is_target)
    if (nrow(target_tx) > 0)
      p <- p + geom_text(data=target_tx,
                         aes(x=(start_cl+end_cl)/2, y=row+0.38, label=gene_name),
                         color="#1a6e1a", size=label_size,
                         hjust=0.5, vjust=0, fontface="bold.italic")
    if (nrow(other_tx) > 0)
      p <- p + geom_text(data=other_tx,
                         aes(x=(start_cl+end_cl)/2, y=row+0.38, label=gene_name),
                         color="#888888", size=label_size,
                         hjust=0.5, vjust=0, fontface="italic")
  } else {
    # Gene absent from annotation: draw a label at TSS position
    p <- p + annotate("text", x=(xmin+xmax)/2, y=1, label=gene,
                      color="#1a6e1a", size=label_size,
                      hjust=0.5, vjust=0, fontface="bold.italic")
  }
  
  p <- p + scale_color_manual(values=gene_colors, guide="none")
  
  if (has_ccre) {
    p <- p +
      geom_hline(yintercept=0, color="grey80", linewidth=0.3, linetype="dashed") +
      geom_rect(data=ccre_df,
                aes(xmin=start_cl, xmax=end_cl, ymin=ccre_ymin, ymax=ccre_ymax, fill=fill_id),
                color=NA, inherit.aes=FALSE) +
      geom_text(data=data.frame(x=xmin, y=(ccre_ymin+ccre_ymax)/2, label="cCREs"),
                aes(x=x, y=y, label=label),
                hjust=0, size=2.8, color="grey40", inherit.aes=FALSE) +
      scale_fill_manual(
        values=unified_fills, breaks=names(ccre_colors),
        labels=c("PLS"="Promoter (PLS)","pELS"="Proximal enhancer (pELS)",
                 "dELS"="Distal enhancer (dELS)","CTCF-only"="CTCF insulator",
                 "DNase-H3K4me3"="DNase+H3K4me3"),
        name="ENCODE cCRE",
        guide=guide_legend(override.aes=list(size=3),
                           keyheight=unit(0.35,"cm"), keywidth=unit(0.35,"cm")))
  } else {
    p <- p + scale_fill_manual(values=unified_fills, guide="none")
  }
  
  p <- p +
    scale_x_continuous(limits=c(xmin,xmax), expand=expansion(mult=c(0,0)),
                       labels=scales::label_number(scale=1e-6,suffix=" Mb",accuracy=0.01)) +
    scale_y_continuous(limits=c(y_floor, n_rows+0.9), expand=c(0,0)) +
    labs(x=NULL, y="Genes") +
    theme_minimal(base_size=10) +
    theme(panel.grid=element_blank(),
          axis.line.x=element_line(color="grey70"),
          axis.ticks.x=element_line(color="grey70"),
          axis.text.x=element_text(size=8),
          axis.title.y=element_text(size=9, angle=90),
          axis.text.y=element_blank(), axis.ticks.y=element_blank(),
          legend.position=if(has_ccre)"right" else "none",
          legend.key.size=unit(0.35,"cm"),
          legend.text=element_text(size=7), legend.title=element_text(size=8),
          plot.margin=margin(2,2,6,2))
  
  list(plot=p, height_in=h_in)
}

# ── Plot loop ─────────────────────────────────────────────────────────────────
failed  <- c()
success <- c()

for (i in seq_len(nrow(top_targets))) {
  
  row     <- top_targets[i, ]
  gene    <- row$Gene
  ct      <- row$CellType
  tf      <- row$TF
  dir_lbl <- row$direction
  ct_safe <- gsub("[^A-Za-z0-9]", "_", ct)
  tf_safe <- gsub("[^A-Za-z0-9]", "_", tf)
  
  outfile <- file.path(multi_diff_dir,
                       paste0(ct_safe,"__", tf_safe,"__", gene,".pdf"))
  
  if (file.exists(outfile)) {
    message("  Skip (exists): ", gene, " [", tf, " / ", ct, "]")
    success <- c(success, gene)
    next
  }
  
  obj_ct <- tryCatch(
    subset(obj_male, subset = cell_type == ct),
    error = function(e) NULL
  )
  if (is.null(obj_ct) || ncol(obj_ct) < 20) {
    message("  Skip [", gene, " / ", ct, "] — too few cells")
    next
  }
  
  obj_ct$genotype <- factor(obj_ct$genotype, levels = c("Ctrl", "KO"))
  obj_ct <- clean_annotation(obj_ct)
  
  gene_resolved <- resolve_gene_alias(gene, Annotation(obj_ct))
  if (gene_resolved != gene)
    message("  Alias resolved: ", gene, " \u2192 ", gene_resolved)
  
  in_annot <- sum(!is.na(Annotation(obj_ct)$gene_name) &
                    Annotation(obj_ct)$gene_name == gene_resolved, na.rm = TRUE)
  if (in_annot == 0) {
    gr_check <- gene_to_region(gene_resolved)
    if (is.null(gr_check)) {
      message("  Skip [", gene, "]: not in annotation or TxDb — likely pseudogene/lncRNA")
      failed <- c(failed, gene)
      next
    }
  }
  
  tryCatch({
    
    obj_ctrl <- clean_annotation(subset(obj_ct, subset = genotype == "Ctrl"))
    obj_ko   <- clean_annotation(subset(obj_ct, subset = genotype == "KO"))
    obj_ctrl$genotype <- factor("Ctrl", levels = "Ctrl")
    obj_ko$genotype   <- factor("KO",   levels = "KO")
    
    df_ctrl <- get_cov_df(obj_ctrl, gene_resolved, "Ctrl", 10000, 10000)
    df_ko   <- get_cov_df(obj_ko,   gene_resolved, "KO",   10000, 10000)
    
    if (is.null(df_ctrl) && is.null(df_ko)) stop("Could not extract coverage data")
    
    n_ctrl <- ncol(obj_ctrl)
    n_ko   <- ncol(obj_ko)
    
    ref_x    <- if (!is.null(df_ko)) df_ko$x else df_ctrl$x
    x_range  <- range(ref_x)
    n_bins   <- length(ref_x)
    common_x <- seq(x_range[1], x_range[2], length.out = n_bins)
    
    y_ctrl <- if (!is.null(df_ctrl)) {
      approx(df_ctrl$x, df_ctrl$y / n_ctrl, xout = common_x, rule = 2)$y
    } else {
      rep(0, n_bins)
    }
    
    y_ko <- if (!is.null(df_ko)) {
      approx(df_ko$x, df_ko$y / n_ko, xout = common_x, rule = 2)$y
    } else {
      rep(0, n_bins)
    }
    
    gaussian_smooth <- function(y, sigma = 3) {
      hw  <- ceiling(3 * sigma)
      ker <- dnorm(-hw:hw, sd = sigma)
      ker <- ker / sum(ker)
      as.numeric(stats::filter(y, ker, circular = TRUE))
    }
    
    y_ctrl <- gaussian_smooth(as.numeric(y_ctrl), sigma = 8)
    y_ko   <- gaussian_smooth(as.numeric(y_ko),   sigma = 8)
    
    diff_df <- data.frame(
      x    = common_x,
      diff = as.numeric(y_ko) - as.numeric(y_ctrl)
    )
    
    split_segments <- function(x, y) {
      segs  <- list()
      cur_x <- c()
      cur_y <- c()
      for (k in seq_len(length(x))) {
        if (k > 1) {
          y0 <- y[k - 1]; y1 <- y[k]
          if (!is.na(y0) && !is.na(y1) && sign(y0) != sign(y1) && (y1 - y0) != 0) {
            x_zero <- x[k - 1] + (x[k] - x[k - 1]) * (-y0 / (y1 - y0))
            cur_x  <- c(cur_x, x_zero)
            cur_y  <- c(cur_y, 0)
            segs[[length(segs) + 1]] <- data.frame(x = cur_x, y = cur_y)
            cur_x <- c(x_zero)
            cur_y <- c(0)
          }
        }
        cur_x <- c(cur_x, x[k])
        cur_y <- c(cur_y, y[k])
      }
      if (length(cur_x) > 0) segs[[length(segs) + 1]] <- data.frame(x = cur_x, y = cur_y)
      segs
    }
    
    segments   <- split_segments(diff_df$x, diff_df$diff)
    seg_pos_df <- dplyr::bind_rows(lapply(segments, function(s) if (mean(s$y, na.rm = TRUE) >= 0) s else NULL))
    seg_neg_df <- dplyr::bind_rows(lapply(segments, function(s) if (mean(s$y, na.rm = TRUE) <  0) s else NULL))
    
    yabs <- max(abs(diff_df$diff), na.rm = TRUE) * 1.05
    if (!is.finite(yabs) || yabs == 0) yabs <- 1
    
    dir_color_hex <- dplyr::case_when(
      grepl("up",   dir_lbl, ignore.case = TRUE) ~ "#b2182b",
      grepl("down", dir_lbl, ignore.case = TRUE) ~ "#2166ac",
      TRUE ~ "#555555"
    )
    
    p_cov <- ggplot(diff_df, aes(x = x, y = diff)) +
      geom_hline(yintercept = 0, color = "grey50", linewidth = 0.4) +
      geom_area(data = seg_pos_df, aes(x = x, y = y),
                fill = diff_colors["UP"],   alpha = 0.75, color = NA) +
      geom_area(data = seg_neg_df, aes(x = x, y = y),
                fill = diff_colors["DOWN"], alpha = 0.75, color = NA) +
      scale_x_continuous(expand = expansion(mult = c(0, 0))) +
      scale_y_continuous(
        limits = c(-yabs, yabs),
        expand = expansion(mult = c(0, 0)),
        labels = function(x) as.character(abs(x))
      ) +
      labs(
        y     = "\u0394 Coverage per cell (KO \u2212 Ctrl)",
        x     = NULL,
        title = paste0(
          "**", gene, "**  |  Regulated by: ", tf,
          "  |  <span style='color:", dir_color_hex, "'>**Expression: ",
          dir_lbl, "**</span>"
        ),
        subtitle = paste0(
          "Cell type: ", ct,
          "  |  Differential ATAC (KO \u2212 Ctrl, male)  |  \u00b110kb",
          "  |  <span style='color:#b2182b'>**\u25a0 KO > Ctrl**</span>",
          "  <span style='color:#2166ac'>**\u25a0 KO < Ctrl**</span>"
        )
      ) +
      theme_minimal(base_size = 11) +
      theme(
        plot.title      = ggtext::element_markdown(size = 12),
        plot.subtitle   = ggtext::element_markdown(size = 9, color = "grey40"),
        panel.grid      = element_blank(),
        axis.line.y     = element_line(color = "grey70"),
        axis.ticks.y    = element_line(color = "grey70"),
        axis.text.x     = element_blank(),
        axis.ticks.x    = element_blank(),
        legend.position = "none",
        plot.margin     = margin(2, 2, 2, 2)
      )
    
    p_peaks <- tryCatch(
      PeakPlot(
        object            = obj_ct,
        region            = gene_resolved,
        assay             = "peaks",
        extend.upstream   = 10000,
        extend.downstream = 10000
      ) + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()),
      error = function(e) NULL
    )
    
    gene_track <- make_gene_track(obj_ct, gene_resolved, 10000, 10000,
                                  ccre_gr       = ccre_gr,
                                  row_height_in = 0.35,
                                  label_size    = 3)
    p_annot <- gene_track$plot
    annot_h <- gene_track$height_in
    
    peaks_h <- 0.5
    cov_h   <- 4.0
    title_h <- 0.8
    
    if (!is.null(p_peaks) && !is.null(p_annot)) {
      p_final <- p_cov / p_peaks / p_annot +
        plot_layout(heights = unit(c(cov_h, peaks_h, annot_h), "in"))
      total_h <- title_h + cov_h + peaks_h + annot_h + 0.5
    } else if (!is.null(p_annot)) {
      p_final <- p_cov / p_annot +
        plot_layout(heights = unit(c(cov_h, annot_h), "in"))
      total_h <- title_h + cov_h + annot_h + 0.5
    } else if (!is.null(p_peaks)) {
      p_final <- p_cov / p_peaks +
        plot_layout(heights = unit(c(cov_h, peaks_h), "in"))
      total_h <- title_h + cov_h + peaks_h + 0.5
    } else {
      p_final <- p_cov
      total_h <- title_h + cov_h + 0.5
    }
    
    ggsave(outfile, p_final, width = 10, height = total_h, limitsize = FALSE)
    success <- c(success, gene)
    message("  Saved: ", gene, " [", tf, " \u2192 ", ct, ", ", dir_lbl, "]")
    
  }, error = function(e) {
    failed <<- c(failed, gene)
    message("  FAILED [", gene, "]: ", e$message)
  })
}

message("\n\u2500\u2500 Step 8D multi-eRegulon summary \u2500\u2500")
message("  Saved:  ", length(success))
message("  Failed: ", length(failed), " \u2014 ",
        if (length(failed) > 0) paste(failed, collapse = ", ") else "none")
message("\n=== Step 8D multi-eRegulon complete. ===")





###################################### ------
###################################### -----



####### INDIVIDUAL COVERAGE PLOTS #######
# ══════════════════════════════════════════════════════════════════════════════
# STEP 7 (REVISED): COVERAGE PLOTS AT BIOLOGICALLY MEANINGFUL LOCI
# ══════════════════════════════════════════════════════════════════════════════

DefaultAssay(obj_male) <- "peaks"

# ── Build interleaved ct_geno factor ─────────────────────────────────────────
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

# ── Interleaved color vector — Ctrl = lighter, KO = full ─────────────────────
ct_geno_colors <- unlist(lapply(ct_levels, function(ct) {
  base_col <- atac_colors[ct]
  if (is.na(base_col)) base_col <- "grey60"
  setNames(
    c(scales::alpha(base_col, 0.45), base_col),
    c(paste(ct, "Ctrl", sep = " \u2014 "),
      paste(ct, "KO",   sep = " \u2014 "))
  )
}))
ct_geno_colors <- ct_geno_colors[names(ct_geno_colors) %in% levels(obj_male$ct_geno)]

# ── Gene sets ─────────────────────────────────────────────────────────────────
# Group 1: SCENIC+ master regulators — plot at their own locus
scenic_tf_genes <- c("Msi2", "Satb2", "Cux2", "Zfp536", "Celf4", "Celf5",
                     "Srrm3", "Mef2c", "Bcl11b", "Rora", "Zmat4", "Zfp831",
                     "Esrrg", "Pparg", "Zfhx3", "Zbtb20", "Fli1")

# Group 2: Top DEGs from your SCENIC+ eRegulon target genes that are also
# biologically meaningful — pull from your DEG results
# These are the genes *regulated by* the master TFs, not the TFs themselves
deg_loci <- c("Pvalb", "Sst", "Vip",        # interneuron markers
              "Mbp", "Plp1",                  # OL/myelin
              "Gfap", "Aldh1l1",              # astrocyte
              "Slc17a7", "Satb2",             # excitatory neuron
              "P2ry12", "Cx3cr1")             # microglia

all_genes <- unique(c(scenic_tf_genes, deg_loci))
message("Coverage plots planned: ", length(all_genes), " loci")

# ── Output directory ──────────────────────────────────────────────────────────
dir.create("coverage_plots", showWarnings = FALSE)
n_tracks <- nlevels(obj_male$ct_geno)

# ── Plot loop ─────────────────────────────────────────────────────────────────
failed  <- c()
success <- c()

for (gene in all_genes) {
  outfile <- file.path("coverage_plots", paste0(gene, "_coverage.pdf"))
  if (file.exists(outfile)) { message("  Skip (exists): ", gene); success <- c(success, gene); next }
  
  tryCatch({
    p_cov <- CoveragePlot(
      object            = obj_male,
      region            = gene,
      assay             = "peaks",
      group.by          = "ct_geno",
      annotation        = TRUE,
      peaks             = TRUE,
      links             = TRUE,
      extend.upstream   = 5000,
      extend.downstream = 5000
    ) & scale_fill_manual(values = ct_geno_colors)
    
    ggsave(outfile, p_cov, width = 14, height = max(8, n_tracks * 0.5))
    success <- c(success, gene)
    message("  Saved: ", gene)
    
  }, error = function(e) { failed <<- c(failed, gene); message("  FAILED [", gene, "]: ", e$message) })
}

message("\n── Coverage plot summary ──")
message("  Saved:  ", length(success), " — ", paste(success,  collapse = ", "))
message("  Failed: ", length(failed),  " — ", paste(failed, collapse = ", "))

# ── Optional: focused single-cell-type plots for key TFs ─────────────────────
# For a paper figure you likely want just the relevant cell type per TF,
# not all 22 tracks. This loops over TF + its cell type from SCENIC+.

tf_ct_pairs <- list(
  Msi2    = "Astrocytes",
  Satb2   = "Layer 5/6 IT neurons",
  Cux2    = "Layer 2/3 IT neurons",
  Zfp536  = "Oligodendrocytes",
  Celf4   = "VIP+ interneurons",
  Celf5   = "Deep-layer extratelencephalic neurons",
  Srrm3   = "PV+ interneurons",
  Mef2c   = "SST+ interneurons",
  Bcl11b  = "Layer 5b PT neurons",
  Rora    = "Oligodendrocyte precursor cells",
  Zmat4   = "Layer 4 sensory neurons"
)

dir.create("coverage_plots/focused", showWarnings = FALSE)

for (tf in names(tf_ct_pairs)) {
  ct <- tf_ct_pairs[[tf]]
  outfile <- file.path("coverage_plots/focused",
                       paste0(tf, "_", gsub("[^A-Za-z0-9]", "_", ct), "_coverage.pdf"))
  if (file.exists(outfile)) { message("  Skip (exists): ", tf); next }
  
  obj_ct <- tryCatch(subset(obj_male, subset = cell_type == ct), error = function(e) NULL)
  if (is.null(obj_ct) || ncol(obj_ct) < 20) { message("  Skip [", tf, "] — too few cells"); next }
  
  obj_ct$ct_geno <- factor(
    paste(obj_ct$cell_type, obj_ct$genotype, sep = " \u2014 "),
    levels = c(paste(ct, "Ctrl", sep = " \u2014 "),
               paste(ct, "KO",   sep = " \u2014 "))
  )
  
  focus_colors <- setNames(
    c(scales::alpha(atac_colors[ct], 0.45), atac_colors[ct]),
    c(paste(ct, "Ctrl", sep = " \u2014 "),
      paste(ct, "KO",   sep = " \u2014 "))
  )
  
  tryCatch({
    p_focus <- CoveragePlot(
      object            = obj_ct,
      region            = tf,
      assay             = "peaks",
      group.by          = "ct_geno",
      annotation        = TRUE,
      peaks             = TRUE,
      links             = TRUE,
      extend.upstream   = 5000,
      extend.downstream = 5000
    ) & scale_fill_manual(values = focus_colors)
    
    ggsave(outfile, p_focus, width = 10, height = 6)
    message("  Focused plot saved: ", tf, " in ", ct)
    
  }, error = function(e) message("  FAILED [", tf, "]: ", e$message))
}

message("\n── Coverage plot summary ──")
message("  Saved:  ", length(success), " — ", paste(success, collapse = ", "))
message("  Failed: ", length(failed),  " — ", paste(failed,  collapse = ", "))
message("\n=== Step 7 complete. ===")




########## Diagnostics #1 ------
# ── Diagnose failed genes ─────────────────────────────────────────────────────
failed_genes <- c("Gm11266", "Gm12239", "Gm37679", "Gm37885", "Hdhd5", "Mfsd14a",
                  "Myorg", "Prxl2a", "Zcchc14", "Atosa", "Epb41l4a", "Nectin1",
                  "Ppp4r1", "Prag1", "Prox1os", "Tprn", "Zfhx2os", "Adgrb2", "Bbln",
                  "Dusp14", "Fam131a", "Gm10635", "Gm31356", "Gm34583", "Gm35161",
                  "Gm35256", "Gm48530", "Gm5468.1", "Gm61537", "Gm64463", "Lemd1",
                  "Mtarc2", "Pnma8a", "Rasal1", "Smim43", "Tmem254", "Calm1", "Cracdl",
                  "Gm56999", "Gm57076", "Gm64463", "Josd1", "Marchf1", "Rundc3a",
                  "Tuba4a", "Vmn2r86", "Gm13629", "Gm30382", "Gucy1a1", "Inka2",
                  "Nectin1", "Nectin3", "9630014M24Rik", "Adgra1", "Akain1", "Ark2c",
                  "C130073E24Rik", "Caln1", "Gm13629", "Gm15246.1", "Gm28175",
                  "Gm29674", "Gm32200", "Gm36736", "Gm38604", "Gm42616", "Gm43606",
                  "Gm45837", "Gm49575", "Gm5820", "Gm64463", "Pcsk2os2", "Pnma8a",
                  "Ptchd1as", "Smim43")

failed_unique <- unique(failed_genes)

annot <- Annotation(obj_male)

for (g in failed_unique) {
  # Get cell type from top_targets
  ct <- top_targets$CellType[top_targets$Gene == g][1]
  
  obj_ct <- tryCatch(subset(obj_male, subset = cell_type == ct), error = function(e) NULL)
  n_cells <- if (!is.null(obj_ct)) ncol(obj_ct) else NA
  n_ctrl  <- if (!is.null(obj_ct)) sum(obj_ct$genotype == "Ctrl") else NA
  n_ko    <- if (!is.null(obj_ct)) sum(obj_ct$genotype == "KO")   else NA
  
  # Check annotation
  in_annot <- sum(!is.na(annot$gene_name) & annot$gene_name == g)
  
  # Check if gene resolves as a region in Signac
  region_ok <- tryCatch({
    StringToGRanges(g, sep = c("-", "-"))
    TRUE
  }, error = function(e) {
    tryCatch({
      LookupGeneCoords(obj_male, gene = g)
      TRUE
    }, error = function(e2) FALSE)
  })
  
  message(sprintf("%-20s | ct: %-30s | cells: %3s (ctrl:%3s ko:%3s) | annot: %3s | region_ok: %s",
                  g, ct, n_cells, n_ctrl, n_ko, in_annot, region_ok))
}
########## Diagnostics #2 ------
# ══════════════════════════════════════════════════════════════════════════════
# STEP 8D — RETRY: failed genes only
# ══════════════════════════════════════════════════════════════════════════════

still_failed <- unique(c(
  "Gm12239", "Gm37679", "Gm37885", "Zcchc14", "Ppp4r1", "Tprn", "Dusp14",
  "Fam131a", "Gm31356", "Gm34583", "Gm35161", "Gm35256", "Gm48530",
  "Gm5468.1", "Gm61537", "Gm64463", "Lemd1", "Rasal1", "Calm1", "Gm56999",
  "Gm57076", "Josd1", "Rundc3a", "Tuba4a", "Vmn2r86", "Gm30382",
  "9630014M24Rik", "C130073E24Rik", "Caln1", "Gm15246.1", "Gm28175",
  "Gm29674", "Gm36736", "Gm38604", "Gm42616", "Gm43606", "Gm45837", "Gm49575"
))

retry_targets <- top_targets[top_targets$Gene %in% still_failed, ]
failed2  <- c()
success2 <- c()

for (i in seq_len(nrow(retry_targets))) {
  
  row     <- retry_targets[i, ]
  gene    <- row$Gene
  ct      <- row$CellType
  tf      <- row$TF
  dir_lbl <- row$direction
  ct_safe <- gsub("[^A-Za-z0-9]", "_", ct)
  tf_safe <- gsub("[^A-Za-z0-9]", "_", tf)
  
  outfile <- file.path("coverage_plots_targets/differential",
                       paste0(gene, "__", tf_safe, "__", ct_safe, ".pdf"))
  
  if (file.exists(outfile)) {
    message("  Skip (exists): ", gene)
    success2 <- c(success2, gene)
    next
  }
  
  obj_ct <- tryCatch(subset(obj_male, subset = cell_type == ct), error = function(e) NULL)
  if (is.null(obj_ct) || ncol(obj_ct) < 20) {
    message("  Skip [", gene, "] — too few cells"); next
  }
  
  obj_ct$genotype <- factor(obj_ct$genotype, levels = c("Ctrl", "KO"))
  obj_ct <- clean_annotation(obj_ct)
  
  gene_resolved <- resolve_gene_alias(gene, Annotation(obj_ct))
  if (gene_resolved != gene)
    message("  Alias resolved: ", gene, " \u2192 ", gene_resolved)
  
  in_annot <- sum(!is.na(Annotation(obj_ct)$gene_name) &
                    Annotation(obj_ct)$gene_name == gene_resolved, na.rm = TRUE)
  if (in_annot == 0) {
    gr_check <- gene_to_region(gene_resolved)
    if (is.null(gr_check)) {
      message("  Skip [", gene, "]: not in annotation or TxDb")
      failed2 <- c(failed2, gene)
      next
    }
  }
  
  last_ok <- "init"
  tryCatch({
    
    obj_ctrl <- clean_annotation(subset(obj_ct, subset = genotype == "Ctrl"))
    obj_ko   <- clean_annotation(subset(obj_ct, subset = genotype == "KO"))
    obj_ctrl$genotype <- factor("Ctrl", levels = "Ctrl")
    obj_ko$genotype   <- factor("KO",   levels = "KO")
    last_ok <- "subset"
    
    df_ctrl <- get_cov_df(obj_ctrl, gene_resolved, "Ctrl", 10000, 10000)
    df_ko   <- get_cov_df(obj_ko,   gene_resolved, "KO",   10000, 10000)
    if (is.null(df_ctrl) && is.null(df_ko)) stop("no coverage data")
    last_ok <- "get_cov"
    
    n_ctrl <- ncol(obj_ctrl); n_ko <- ncol(obj_ko)
    ref_x    <- if (!is.null(df_ko)) df_ko$x else df_ctrl$x
    common_x <- seq(range(ref_x)[1], range(ref_x)[2], length.out = length(ref_x))
    y_ctrl <- if (!is.null(df_ctrl)) approx(df_ctrl$x, df_ctrl$y / n_ctrl, xout = common_x, rule = 2)$y else rep(0, length(common_x))
    y_ko   <- if (!is.null(df_ko))   approx(df_ko$x,   df_ko$y   / n_ko,   xout = common_x, rule = 2)$y else rep(0, length(common_x))
    y_ctrl <- stats::filter(y_ctrl, rep(1/3, 3), circular = TRUE)
    y_ko   <- stats::filter(y_ko,   rep(1/3, 3), circular = TRUE)
    last_ok <- "grid"
    
    diff_df <- data.frame(x = common_x, diff = as.numeric(y_ko) - as.numeric(y_ctrl))
    last_ok <- "diff_df"
    
    split_segments <- function(x, y) {
      segs <- list(); cur_x <- c(); cur_y <- c()
      for (k in seq_len(length(x))) {
        if (k > 1) {
          y0 <- y[k-1]; y1 <- y[k]
          if (!is.na(y0) && !is.na(y1) && sign(y0) != sign(y1) && (y1-y0) != 0) {
            x_zero <- x[k-1] + (x[k]-x[k-1]) * (-y0/(y1-y0))
            cur_x <- c(cur_x, x_zero); cur_y <- c(cur_y, 0)
            segs[[length(segs)+1]] <- data.frame(x=cur_x, y=cur_y)
            cur_x <- c(x_zero); cur_y <- c(0) } }
        cur_x <- c(cur_x, x[k]); cur_y <- c(cur_y, y[k]) }
      if (length(cur_x) > 0) segs[[length(segs)+1]] <- data.frame(x=cur_x, y=cur_y)
      segs
    }
    
    segments   <- split_segments(diff_df$x, diff_df$diff)
    seg_pos_df <- dplyr::bind_rows(lapply(segments, function(s) if (mean(s$y, na.rm=TRUE) >= 0) s else NULL))
    seg_neg_df <- dplyr::bind_rows(lapply(segments, function(s) if (mean(s$y, na.rm=TRUE) <  0) s else NULL))
    last_ok <- "segments"
    
    yabs <- max(abs(diff_df$diff), na.rm = TRUE) * 1.05
    if (!is.finite(yabs) || yabs == 0) yabs <- 1
    last_ok <- "yabs"
    
    dir_color_hex <- dplyr::case_when(
      grepl("up",   dir_lbl, ignore.case = TRUE) ~ "#b2182b",
      grepl("down", dir_lbl, ignore.case = TRUE) ~ "#2166ac",
      TRUE ~ "#555555"
    )
    
    p_cov <- ggplot(diff_df, aes(x = x, y = diff)) +
      geom_hline(yintercept = 0, color = "grey50", linewidth = 0.4) +
      geom_area(data = seg_pos_df, aes(x = x, y = y),
                fill = diff_colors["UP"],   alpha = 0.75, color = NA) +
      geom_area(data = seg_neg_df, aes(x = x, y = y),
                fill = diff_colors["DOWN"], alpha = 0.75, color = NA) +
      scale_x_continuous(expand = expansion(mult = c(0, 0))) +
      scale_y_continuous(
        limits = c(-yabs, yabs), expand = expansion(mult = c(0, 0)),
        labels = function(x) as.character(abs(x))
      ) +
      labs(
        y     = "\u0394 Coverage per cell (KO \u2212 Ctrl)",
        x     = NULL,
        title = paste0("**", gene, "**  |  Regulated by: ", tf,
                       "  |  <span style='color:", dir_color_hex,
                       "'>**Expression: ", dir_lbl, "**</span>"),
        subtitle = paste0("Cell type: ", ct,
                          "  |  Differential ATAC (KO \u2212 Ctrl, male)  |  \u00b110kb",
                          "  |  <span style='color:#b2182b'>**\u25a0 KO > Ctrl**</span>",
                          "  <span style='color:#2166ac'>**\u25a0 KO < Ctrl**</span>")
      ) +
      theme_minimal(base_size = 11) +
      theme(
        plot.title      = ggtext::element_markdown(size = 12),
        plot.subtitle   = ggtext::element_markdown(size = 9, color = "grey40"),
        panel.grid      = element_blank(),
        axis.line.y     = element_line(color = "grey70"),
        axis.ticks.y    = element_line(color = "grey70"),
        axis.text.x     = element_blank(),
        axis.ticks.x    = element_blank(),
        legend.position = "none",
        plot.margin     = margin(2, 2, 2, 2)
      )
    last_ok <- "p_cov"
    
    p_peaks <- tryCatch(
      PeakPlot(obj_ct, region = gene_resolved, assay = "peaks",
               extend.upstream = 10000, extend.downstream = 10000) +
        theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()),
      error = function(e) NULL
    )
    last_ok <- "p_peaks"
    
    gene_track <- make_gene_track(obj_ct, gene_resolved, 10000, 10000,
                                  row_height_in = 0.35, label_size = 3)
    p_annot <- gene_track$plot
    annot_h <- gene_track$height_in
    last_ok <- "gene_track"
    
    annot_gene <- Annotation(obj_ct)
    target_gr  <- annot_gene[!is.na(annot_gene$gene_name) &
                               annot_gene$gene_name == gene_resolved]
    p_ccre <- NULL; ccre_h <- 0
    if (length(target_gr) > 0) {
      chr_plot  <- as.character(seqnames(target_gr[1]))
      tss_plot  <- ifelse(as.character(strand(target_gr[1])) == "+",
                          min(start(target_gr)), max(end(target_gr)))
      ct_res    <- make_ccre_track(ccre_gr, tss_plot - 10000, tss_plot + 10000, chr_plot)
      p_ccre    <- ct_res$plot
      ccre_h    <- ct_res$height_in
    }
    last_ok <- "ccre_track"
    
    peaks_h <- 0.5; cov_h <- 4.0; title_h <- 0.8
    panels  <- list(p_cov, p_peaks, p_ccre, p_annot)
    heights <- c(cov_h, peaks_h, ccre_h, annot_h)
    keep    <- !sapply(panels, is.null) & heights > 0
    panels  <- panels[keep]
    heights <- heights[keep]
    
    p_final <- Reduce(`/`, panels) +
      plot_layout(heights = unit(heights, "in"))
    total_h <- title_h + sum(heights) + 0.5
    last_ok <- "assemble"
    
    ggsave(outfile, p_final, width = 10, height = total_h, limitsize = FALSE)
    last_ok <- "ggsave"
    
    success2 <- c(success2, gene)
    message("  Saved: ", gene, " [last_ok: ", last_ok, "]")
    
  }, error = function(e) {
    failed2 <<- c(failed2, gene)
    message("  FAILED [", gene, "] at step '", last_ok, "': ", e$message)
  })
}

message("\n\u2500\u2500 Retry summary \u2500\u2500")
message("  Saved:  ", length(success2))
message("  Failed: ", length(failed2), " \u2014 ",
        if (length(failed2) > 0) paste(failed2, collapse = ", ") else "none")
########## Diagnostics #3 ------
# # Check if links exist
# Links(obj_ct)
########## Diagnostics #4 ------
head(ccre_raw$type, 20)
table(ccre_raw$type)

# Inspect raw file
lines <- readLines(ccre_bed, n = 5)
cat(lines, sep = "\n")

# Count columns in first row
strsplit(lines[2], "\t")[[1]]

########## Diagnostics #5 -----
obj_ct <- subset(obj_male, subset = cell_type == "Astrocytes")
obj_ct <- clean_annotation(obj_ct)
annot  <- Annotation(obj_ct)

gene_resolved <- resolve_gene_alias("Myorg", annot)
message("Resolved: ", gene_resolved)
message("In annot: ", sum(!is.na(annot$gene_name) & annot$gene_name == gene_resolved))
gene_to_region("Myorg")
########## Difference #2 ------
# ══════════════════════════════════════════════════════════════════════════════
# STEP 8D: DIFFERENTIAL COVERAGE — KO minus Ctrl
# ══════════════════════════════════════════════════════════════════════════════

DefaultAssay(obj_male) <- "peaks"

further_dir <- "~/Downloads/Seurat_scATAC-seq/scenicplus_output/scenicplus/scenicplus_further_analysis"
dir.create("coverage_plots_targets/differential", showWarnings = FALSE, recursive = TRUE)

diff_colors <- c("UP" = "#b2182b", "DOWN" = "#2166ac")

# ── Helper: extract x,y from a single-group CoveragePlot ─────────────────────
get_cov_df <- function(obj_sub, gene, label, extend_up, extend_dn) {
  p <- CoveragePlot(
    object            = obj_sub,
    region            = gene,
    assay             = "peaks",
    group.by          = "genotype",
    annotation        = FALSE,
    peaks             = FALSE,
    links             = FALSE,
    extend.upstream   = extend_up,
    extend.downstream = extend_dn
  )
  built <- tryCatch(ggplot_build(p), error = function(e) NULL)
  if (is.null(built)) return(NULL)
  df <- tryCatch({
    d <- built$data[[1]]
    if (nrow(d) > 0 && all(c("x", "y") %in% colnames(d))) d else NULL
  }, error = function(e) NULL)
  if (is.null(df)) return(NULL)
  data.frame(x = df$x, y = df$y, genotype = label, stringsAsFactors = FALSE)
}

# ── Helper: clean NA coords from annotation ───────────────────────────────────
clean_annotation <- function(obj) {
  annot <- Annotation(obj)
  if (!is.null(annot) && length(annot) > 0) {
    keep <- !is.na(start(annot)) & !is.na(end(annot))
    Annotation(obj) <- annot[keep]
  }
  obj
}

# ── Helper: gene track with composite exon model ─────────────────────────────
make_gene_track <- function(obj_ct, gene, extend_up, extend_dn,
                            row_height_in = 0.35, label_size = 3) {
  
  annot <- Annotation(obj_ct)
  if (is.null(annot) || length(annot) == 0) return(list(plot = NULL, height_in = 0))
  
  target <- annot[!is.na(annot$gene_name) & annot$gene_name == gene]
  if (length(target) == 0) return(list(plot = NULL, height_in = 0))
  
  chr  <- as.character(seqnames(target[1]))
  tss  <- ifelse(as.character(strand(target[1])) == "+",
                 min(start(target)), max(end(target)))
  xmin <- tss - extend_up
  xmax <- tss + extend_dn
  
  in_reg <- annot[
    !is.na(start(annot)) & !is.na(end(annot)) &
      as.character(seqnames(annot)) == chr &
      start(annot) <= xmax & end(annot) >= xmin
  ]
  if (length(in_reg) == 0) return(list(plot = NULL, height_in = 0))
  
  df <- as.data.frame(in_reg, row.names = NULL) |>
    dplyr::mutate(
      gene_name = ifelse(is.na(gene_name), "unknown", as.character(gene_name)),
      strand    = as.character(strand),
      start_cl  = pmax(start, xmin),
      end_cl    = pmin(end,   xmax)
    ) |>
    dplyr::filter(start_cl < end_cl)
  
  if (nrow(df) == 0) return(list(plot = NULL, height_in = 0))
  
  # One backbone per gene: full extent
  transcripts <- df |>
    dplyr::group_by(gene_name) |>
    dplyr::summarise(
      start_cl = min(start_cl),
      end_cl   = max(end_cl),
      strand   = dplyr::first(strand),
      .groups  = "drop"
    ) |>
    dplyr::arrange(start_cl)
  
  # Union of exon/cds intervals per gene (merge overlapping)
  exon_raw <- df |>
    dplyr::filter(type %in% c("exon", "cds")) |>
    dplyr::select(gene_name, start_cl, end_cl) |>
    dplyr::arrange(gene_name, start_cl)
  
  exons <- exon_raw |>
    dplyr::group_by(gene_name) |>
    dplyr::group_modify(~ {
      d <- dplyr::arrange(.x, start_cl)
      merged <- d[1, , drop = FALSE]
      if (nrow(d) > 1) {
        for (j in 2:nrow(d)) {
          if (d$start_cl[j] <= merged$end_cl[nrow(merged)]) {
            merged$end_cl[nrow(merged)] <- max(merged$end_cl[nrow(merged)], d$end_cl[j])
          } else {
            merged <- dplyr::bind_rows(merged, d[j, , drop = FALSE])
          }
        }
      }
      merged
    }) |>
    dplyr::ungroup()
  
  if (nrow(transcripts) == 0) return(list(plot = NULL, height_in = 0))
  
  # Greedy row packer
  transcripts$row <- 1L
  row_ends <- c(transcripts$end_cl[1])
  if (nrow(transcripts) > 1) {
    for (k in 2:nrow(transcripts)) {
      placed <- FALSE
      for (r in seq_along(row_ends)) {
        if (transcripts$start_cl[k] > row_ends[r] + 500) {
          transcripts$row[k] <- as.integer(r)
          row_ends[r] <- transcripts$end_cl[k]
          placed <- TRUE; break
        }
      }
      if (!placed) {
        row_ends <- c(row_ends, transcripts$end_cl[k])
        transcripts$row[k] <- as.integer(length(row_ends))
      }
    }
  }
  
  n_rows <- max(transcripts$row)
  h_in   <- max(0.8, n_rows * row_height_in + 0.7)
  
  row_lookup <- transcripts |>
    dplyr::select(gene_name, row) |>
    dplyr::distinct(gene_name, .keep_all = TRUE)
  
  exons <- exons |>
    dplyr::left_join(row_lookup, by = "gene_name") |>
    dplyr::filter(!is.na(row))
  
  # Direction arrow ticks
  arrow_spacing <- (xmax - xmin) / 20
  arrows_df <- transcripts |>
    dplyr::rowwise() |>
    dplyr::mutate(
      tick_pos = list({
        pos <- seq(start_cl + arrow_spacing / 2,
                   end_cl  - arrow_spacing / 2,
                   by = arrow_spacing)
        if (length(pos) == 0) pos <- (start_cl + end_cl) / 2
        pos
      })
    ) |>
    tidyr::unnest(tick_pos) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      x_from = ifelse(strand == "+", tick_pos - arrow_spacing * 0.3,
                      tick_pos + arrow_spacing * 0.3),
      x_to   = ifelse(strand == "+", tick_pos + arrow_spacing * 0.3,
                      tick_pos - arrow_spacing * 0.3)
    )
  
  all_genes   <- unique(transcripts$gene_name)
  gene_colors <- setNames(
    ifelse(all_genes == gene, "#1a6e1a", "#3333aa"),
    all_genes
  )
  
  p <- ggplot() +
    # Backbone
    geom_segment(data = transcripts,
                 aes(x = start_cl, xend = end_cl, y = row, yend = row,
                     color = gene_name),
                 linewidth = 0.5) +
    # Direction arrows
    geom_segment(data = arrows_df,
                 aes(x = x_from, xend = x_to, y = row, yend = row,
                     color = gene_name),
                 linewidth = 0.4,
                 arrow = arrow(length = unit(0.12, "cm"), type = "open",
                               ends = "last")) +
    # Exon boxes — height reduced to ±0.15 (was ±0.25/0.30)
    geom_rect(data = exons,
              aes(xmin = start_cl, xmax = end_cl,
                  ymin = row - 0.15, ymax = row + 0.15,
                  fill = gene_name),
              color = NA) +
    # Gene labels
    geom_text(data = transcripts,
              aes(x = (start_cl + end_cl) / 2, y = row + 0.35,
                  label = gene_name, color = gene_name),
              size = label_size, hjust = 0.5, vjust = 0, fontface = "italic") +
    scale_color_manual(values = gene_colors, guide = "none") +
    scale_fill_manual(values  = gene_colors, guide = "none") +
    scale_x_continuous(
      limits = c(xmin, xmax),
      expand = expansion(mult = c(0, 0)),
      labels = scales::label_number(scale = 1e-6, suffix = " Mb", accuracy = 0.01)
    ) +
    scale_y_continuous(limits = c(0.3, n_rows + 0.9), expand = c(0, 0)) +
    labs(x = NULL, y = "Genes") +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid   = element_blank(),
      axis.line.x  = element_line(color = "grey70"),
      axis.ticks.x = element_line(color = "grey70"),
      axis.text.x  = element_text(size = 8),
      axis.title.y = element_text(size = 9, angle = 90),
      axis.text.y  = element_blank(),
      axis.ticks.y = element_blank(),
      plot.margin  = margin(2, 2, 6, 2)
    )
  
  list(plot = p, height_in = h_in)
}

# ── Plot loop ──────────────────────────────────────────────────────────────────
failed  <- c()
success <- c()

for (i in seq_len(nrow(top_targets))) {
  
  row     <- top_targets[i, ]
  gene    <- row$Gene
  ct      <- row$CellType
  tf      <- row$TF
  dir_lbl <- row$direction
  ct_safe <- gsub("[^A-Za-z0-9]", "_", ct)
  tf_safe <- gsub("[^A-Za-z0-9]", "_", tf)
  
  outfile <- file.path("coverage_plots_targets/differential",
                       paste0(gene, "__", tf_safe, "__", ct_safe, ".pdf"))
  
  if (file.exists(outfile)) {
    message("  Skip (exists): ", gene, " [", tf, " / ", ct, "]")
    success <- c(success, gene)
    next
  }
  
  obj_ct <- tryCatch(
    subset(obj_male, subset = cell_type == ct),
    error = function(e) NULL
  )
  if (is.null(obj_ct) || ncol(obj_ct) < 20) {
    message("  Skip [", gene, " / ", ct, "] — too few cells")
    next
  }
  
  obj_ct$genotype <- factor(obj_ct$genotype, levels = c("Ctrl", "KO"))
  obj_ct <- clean_annotation(obj_ct)
  
  tryCatch({
    
    obj_ctrl <- clean_annotation(subset(obj_ct, subset = genotype == "Ctrl"))
    obj_ko   <- clean_annotation(subset(obj_ct, subset = genotype == "KO"))
    obj_ctrl$genotype <- factor("Ctrl", levels = "Ctrl")
    obj_ko$genotype   <- factor("KO",   levels = "KO")
    
    df_ctrl <- get_cov_df(obj_ctrl, gene, "Ctrl", 10000, 10000)
    df_ko   <- get_cov_df(obj_ko,   gene, "KO",   10000, 10000)
    
    if (is.null(df_ctrl) && is.null(df_ko)) stop("Could not extract coverage data")
    
    # ── Cell counts for per-cell normalization ────────────────────────────────
    n_ctrl <- ncol(obj_ctrl)
    n_ko   <- ncol(obj_ko)
    
    # ── Align both signals to a shared uniform grid before subtracting ────────
    ref_x <- if (!is.null(df_ko)) df_ko$x else df_ctrl$x
    
    # Build a common grid with the same resolution as the data
    x_range  <- range(ref_x)
    n_bins   <- length(ref_x)
    common_x <- seq(x_range[1], x_range[2], length.out = n_bins)
    
    y_ctrl <- if (!is.null(df_ctrl)) {
      approx(df_ctrl$x, df_ctrl$y / n_ctrl, xout = common_x, rule = 2)$y
    } else {
      rep(0, n_bins)
    }
    
    y_ko <- if (!is.null(df_ko)) {
      approx(df_ko$x, df_ko$y / n_ko, xout = common_x, rule = 2)$y
    } else {
      rep(0, n_bins)
    }
    
    # Optional: smooth slightly to reduce interpolation noise
    y_ctrl <- stats::filter(y_ctrl, rep(1/3, 3), circular = TRUE)
    y_ko   <- stats::filter(y_ko,   rep(1/3, 3), circular = TRUE)
    
    diff_df <- data.frame(
      x    = common_x,
      diff = as.numeric(y_ko) - as.numeric(y_ctrl)
    )
    
    # ── Split into contiguous positive/negative segments ─────────────────────
    # Insert zero-crossing interpolated points so polygons close cleanly at zero
    split_segments <- function(x, y) {
      segs <- list()
      n <- length(x)
      cur_x <- c()
      cur_y <- c()
      for (k in seq_len(n)) {
        if (k > 1) {
          y0 <- y[k - 1]; y1 <- y[k]
          # Sign change — interpolate exact zero crossing
          if (!is.na(y0) && !is.na(y1) && sign(y0) != sign(y1) && (y1 - y0) != 0) {
            x_zero <- x[k - 1] + (x[k] - x[k - 1]) * (-y0 / (y1 - y0))
            cur_x <- c(cur_x, x_zero)
            cur_y <- c(cur_y, 0)
            segs[[length(segs) + 1]] <- data.frame(x = cur_x, y = cur_y)
            cur_x <- c(x_zero)
            cur_y <- c(0)
          }
        }
        cur_x <- c(cur_x, x[k])
        cur_y <- c(cur_y, y[k])
      }
      if (length(cur_x) > 0) segs[[length(segs) + 1]] <- data.frame(x = cur_x, y = cur_y)
      segs
    }
    
    segments     <- split_segments(diff_df$x, diff_df$diff)
    seg_pos_df   <- dplyr::bind_rows(lapply(segments, function(s) if (mean(s$y, na.rm=TRUE) >= 0) s else NULL))
    seg_neg_df   <- dplyr::bind_rows(lapply(segments, function(s) if (mean(s$y, na.rm=TRUE) <  0) s else NULL))
    
    yabs <- max(abs(diff_df$diff), na.rm = TRUE) * 1.05
    if (!is.finite(yabs) || yabs == 0) yabs <- 1
    
    # Normalize direction label for color lookup
    fill_col <- dplyr::case_when(
      grepl("up", dir_lbl, ignore.case = TRUE)   ~ diff_colors["UP"],
      grepl("down", dir_lbl, ignore.case = TRUE) ~ diff_colors["DOWN"],
      TRUE ~ "#555555"
    )
    
    dir_color_hex <- dplyr::case_when(
      grepl("up",   dir_lbl, ignore.case = TRUE) ~ "#b2182b",
      grepl("down", dir_lbl, ignore.case = TRUE) ~ "#2166ac",
      TRUE ~ "#555555"
    )
    
    p_cov <- ggplot(diff_df, aes(x = x, y = diff)) +
      geom_hline(yintercept = 0, color = "grey50", linewidth = 0.4) +
      geom_area(data = seg_pos_df, aes(x = x, y = y),
                fill = diff_colors["UP"],   alpha = 0.75, color = NA) +
      geom_area(data = seg_neg_df, aes(x = x, y = y),
                fill = diff_colors["DOWN"], alpha = 0.75, color = NA) +
      scale_x_continuous(expand = expansion(mult = c(0, 0))) +
      scale_y_continuous(
        limits = c(-yabs, yabs),
        expand = expansion(mult = c(0, 0)),
        labels = function(x) as.character(abs(x))
      ) +
      labs(
        y     = "\u0394 Coverage (KO \u2212 Ctrl)",
        x     = NULL,
        title = paste0(
          "**", gene, "**  |  Regulated by: ", tf,
          "  |  <span style='color:", dir_color_hex, "'>**Expression: ",
          dir_lbl, "**</span>"
        ),
        subtitle = paste0(
          "Cell type: ", ct,
          "  |  Differential ATAC (KO \u2212 Ctrl, male)  |  \u00b110kb",
          "  |  <span style='color:#b2182b'>**\u25a0 KO > Ctrl**</span>",
          "  <span style='color:#2166ac'>**\u25a0 KO < Ctrl**</span>"
        )
      ) +
      theme_minimal(base_size = 11) +
      theme(
        plot.title      = ggtext::element_markdown(size = 12),
        plot.subtitle   = ggtext::element_markdown(size = 9, color = "grey40"),
        panel.grid      = element_blank(),
        axis.line.y     = element_line(color = "grey70"),
        axis.ticks.y    = element_line(color = "grey70"),
        axis.text.x     = element_blank(),
        axis.ticks.x    = element_blank(),
        legend.position = "none",
        plot.margin     = margin(2, 2, 2, 2)
      )
    
    # ── Peaks panel ───────────────────────────────────────────────────────────
    p_peaks <- tryCatch(
      PeakPlot(
        object            = obj_ct,
        region            = gene,
        assay             = "peaks",
        extend.upstream   = 10000,
        extend.downstream = 10000
      ) + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()),
      error = function(e) NULL
    )
    
    # ── Gene track ────────────────────────────────────────────────────────────
    gene_track <- make_gene_track(obj_ct, gene, 10000, 10000,
                                  row_height_in = 0.35, label_size = 3)
    p_annot <- gene_track$plot
    annot_h <- gene_track$height_in
    
    peaks_h <- 0.5
    cov_h   <- 4.0
    title_h <- 0.8
    
    if (!is.null(p_peaks) && !is.null(p_annot)) {
      p_final <- p_cov / p_peaks / p_annot +
        plot_layout(heights = unit(c(cov_h, peaks_h, annot_h), "in"))
      total_h <- title_h + cov_h + peaks_h + annot_h + 0.5
    } else if (!is.null(p_annot)) {
      p_final <- p_cov / p_annot +
        plot_layout(heights = unit(c(cov_h, annot_h), "in"))
      total_h <- title_h + cov_h + annot_h + 0.5
    } else if (!is.null(p_peaks)) {
      p_final <- p_cov / p_peaks +
        plot_layout(heights = unit(c(cov_h, peaks_h), "in"))
      total_h <- title_h + cov_h + peaks_h + 0.5
    } else {
      p_final <- p_cov
      total_h <- title_h + cov_h + 0.5
    }
    
    ggsave(outfile, p_final, width = 10, height = total_h, limitsize = FALSE)
    success <- c(success, gene)
    message("  Saved: ", gene, " [", tf, " \u2192 ", ct, ", ", dir_lbl, "]")
    
  }, error = function(e) {
    failed <<- c(failed, gene)
    message("  FAILED [", gene, "]: ", e$message)
  })
}

message("\n\u2500\u2500 Step 8D summary \u2500\u2500")
message("  Saved:  ", length(success))
message("  Failed: ", length(failed), " \u2014 ",
        if (length(failed) > 0) paste(failed, collapse = ", ") else "none")
message("\n=== Step 8D complete. ===")


######## Coverage plots (initial version) -------
# ══════════════════════════════════════════════════════════════════════════════
# STEP 8: COVERAGE PLOTS AT TF TARGET GENE LOCI
# Shows chromatin openness at DEG loci regulated by SCENIC+ master TFs
# Uses TF_DEG_overlap to identify which target genes to plot per cell type
# ══════════════════════════════════════════════════════════════════════════════

DefaultAssay(obj_male) <- "peaks"

further_dir <- "~/Downloads/Seurat_scATAC-seq/scenicplus_output/scenicplus/scenicplus_further_analysis"
dir.create("coverage_plots_targets", showWarnings = FALSE)
dir.create("coverage_plots_targets/focused", showWarnings = FALSE)

# ── Load TF-DEG overlap ───────────────────────────────────────────────────────
tf_deg <- read.csv(file.path(further_dir,
                             "TF_DEG_overlap/TF_DEG_overlap_permutation_KO_vs_Ctrl.csv"))
message("TF-DEG overlap loaded: ", nrow(tf_deg), " rows")
message("Columns: ", paste(colnames(tf_deg), collapse = " | "))

# ── Identify gene + cell type + TF columns (adjust if names differ) ───────────
# Expected columns based on your Python script:
#   Gene, TF, CellType, Overlap, PercentOverlap, Permutation_pval, fdr_adj_pval
gene_col  <- grep("^gene$|^Gene$|^target$|^Target$", colnames(tf_deg),
                  ignore.case = TRUE, value = TRUE)[1]
tf_col    <- grep("^TF$|^tf$", colnames(tf_deg),
                  ignore.case = TRUE, value = TRUE)[1]
ct_col    <- grep("CellType|cell_type|celltype", colnames(tf_deg),
                  ignore.case = TRUE, value = TRUE)[1]
pval_col  <- grep("fdr|padj|adj", colnames(tf_deg),
                  ignore.case = TRUE, value = TRUE)[1]

message("\nMapped columns:")
message("  gene  = ", gene_col)
message("  TF    = ", tf_col)
message("  CT    = ", ct_col)
message("  pval  = ", pval_col)

# ── Load DEG direction (up/down in KO) ───────────────────────────────────────
deg_dir_path <- file.path(further_dir, "DEG_KO_vs_Ctrl")
deg_all      <- read.csv(file.path(deg_dir_path,
                                   "DEG_all_celltypes_KO_vs_Ctrl.csv"))
# Expected: gene, avg_log2FC (or log2FoldChange), cell_type
deg_gene_col <- grep("^gene$|^Gene$", colnames(deg_all),
                     ignore.case = TRUE, value = TRUE)[1]
deg_lfc_col  <- grep("log2FC|log2FoldChange|avg_log2FC",
                     colnames(deg_all), ignore.case = TRUE, value = TRUE)[1]
deg_ct_col   <- grep("CellType|cell_type", colnames(deg_all),
                     ignore.case = TRUE, value = TRUE)[1]

message("\nDEG columns:")
message("  gene  = ", deg_gene_col)
message("  LFC   = ", deg_lfc_col)
message("  CT    = ", deg_ct_col)

# ── Build interleaved ct_geno factor and colors ───────────────────────────────
ct_levels <- levels(factor(obj_male$cell_type,
                           levels = atac_new_order[atac_new_order %in%
                                                     unique(obj_male$cell_type)]))

interleaved_levels <- as.vector(rbind(
  paste(ct_levels, "Ctrl", sep = " \u2014 "),
  paste(ct_levels, "KO",   sep = " \u2014 ")
))

obj_male$ct_geno <- factor(
  paste(obj_male$cell_type, obj_male$genotype, sep = " \u2014 "),
  levels = interleaved_levels[interleaved_levels %in% unique(obj_male$ct_geno)]
)

ct_geno_colors <- unlist(lapply(ct_levels, function(ct) {
  base_col <- atac_colors[ct]
  if (is.na(base_col)) base_col <- "grey60"
  setNames(
    c(scales::alpha(base_col, 0.45), base_col),
    c(paste(ct, "Ctrl", sep = " \u2014 "),
      paste(ct, "KO",   sep = " \u2014 "))
  )
}))
ct_geno_colors <- ct_geno_colors[names(ct_geno_colors) %in% levels(obj_male$ct_geno)]

# ── All target genes from significantly enriched TF × CT pairs ───────────────
top_targets <- tf_deg |>
  rename(TF = TF, CellType = Cell_Type, fdr = Perm_FDR) |>
  filter(Direction == "All", fdr < 0.05) |>
  mutate(Overlap_Genes = as.character(Overlap_Genes)) |>
  separate_rows(Overlap_Genes, sep = ",\\s*") |>
  rename(Gene = Overlap_Genes) |>
  filter(!is.na(Gene), Gene != "", !grepl("^ENSMUSG", Gene)) |>
  left_join(
    deg_all |>
      rename(Gene = all_of(deg_gene_col),
             LFC  = all_of(deg_lfc_col),
             CT   = all_of(deg_ct_col)) |>
      select(Gene, CT, LFC) |>
      filter(!is.na(LFC)),
    by = c("Gene" = "Gene", "CellType" = "CT")
  ) |>
  mutate(direction = case_when(
    LFC >  0 ~ "Up in KO",
    LFC <  0 ~ "Down in KO",
    TRUE     ~ "NA"
  ))

message("Total plots to generate: ", nrow(top_targets),
        " across ", length(unique(paste(top_targets$TF, top_targets$CellType))),
        " TF × CT pairs")
print(top_targets |> count(TF, CellType) |> arrange(desc(n)))

# ── Add DEG direction label ───────────────────────────────────────────────────
deg_lfc <- deg_all |>
  rename(Gene = all_of(deg_gene_col),
         LFC  = all_of(deg_lfc_col),
         CT   = all_of(deg_ct_col)) |>
  select(Gene, CT, LFC) |>
  distinct()

# top_targets already has LFC and direction from the build block above
# just verify:
message("Columns in top_targets: ", paste(colnames(top_targets), collapse = " | "))

# ── Coverage plot loop — focused (one CT per TF × target gene) ───────────────
failed  <- c()
success <- c()

for (i in seq_len(nrow(top_targets))) {
  
  row     <- top_targets[i, ]
  gene    <- row$Gene
  ct      <- row$CellType
  tf      <- row$TF
  dir_lbl <- row$direction
  ct_safe <- gsub("[^A-Za-z0-9]", "_", ct)
  tf_safe <- gsub("[^A-Za-z0-9]", "_", tf)
  
  outfile <- file.path("coverage_plots_targets/focused",
                       paste0(gene, "__", tf_safe, "__", ct_safe, ".pdf"))
  
  if (file.exists(outfile)) {
    message("  Skip (exists): ", gene, " [", tf, " / ", ct, "]")
    success <- c(success, gene)
    next
  }
  
  # Subset to relevant cell type only
  obj_ct <- tryCatch(
    subset(obj_male, subset = cell_type == ct),
    error = function(e) NULL
  )
  
  if (is.null(obj_ct) || ncol(obj_ct) < 20) {
    message("  Skip [", gene, " / ", ct, "] — too few cells: ", ncol(obj_ct))
    next
  }
  
  obj_ct$ct_geno <- factor(
    paste(obj_ct$cell_type, obj_ct$genotype, sep = " \u2014 "),
    levels = c(paste(ct, "Ctrl", sep = " \u2014 "),
               paste(ct, "KO",   sep = " \u2014 "))
  )
  
  focus_colors <- setNames(
    c(scales::alpha(atac_colors[ct], 0.45), atac_colors[ct]),
    c(paste(ct, "Ctrl", sep = " \u2014 "),
      paste(ct, "KO",   sep = " \u2014 "))
  )
  
  tryCatch({
    p <- CoveragePlot(
      object            = obj_ct,
      region            = gene,
      assay             = "peaks",
      group.by          = "ct_geno",
      annotation        = TRUE,
      peaks             = TRUE,
      links             = TRUE,
      extend.upstream   = 10000,
      extend.downstream = 10000
    ) & scale_fill_manual(values = focus_colors)
    
    # Annotate with TF and direction
    p_final <- patchwork::wrap_elements(p) +
      patchwork::plot_annotation(
        title    = paste0(gene, "  |  Regulated by: ", tf, "  |  ", dir_lbl),
        subtitle = paste0("Cell type: ", ct, "  |  KO vs Ctrl (male)  |  ±10kb window"),
        theme    = theme(
          plot.title    = element_text(face = "bold", size = 12),
          plot.subtitle = element_text(size = 9, color = "grey40")
        )
      )
    
    ggsave(outfile, p_final, width = 10, height = 7)
    success <- c(success, gene)
    message("  Saved: ", gene, " [", tf, " → ", ct, ", ", dir_lbl, "]")
    
  }, error = function(e) {
    failed <<- c(failed, gene)
    message("  FAILED [", gene, "]: ", e$message)
  })
}

# ── Summary ───────────────────────────────────────────────────────────────────
message("\n── Step 8 coverage plots summary ──")
message("  Saved:  ", length(success), " plots")
message("  Failed: ", length(failed),  " — ",
        if (length(failed) > 0) paste(failed, collapse = ", ") else "none")
message("\n=== Step 8 complete. ===")

# ══════════════════════════════════════════════════════════════════════════════
# STEP 8: COVERAGE PLOTS — (Ctrl vs KO, 2 tracks, per gene)
# Uses CoveragePlot() directly — same approach as Step 8 that was working
# ══════════════════════════════════════════════════════════════════════════════

DefaultAssay(obj_male) <- "peaks"

further_dir <- "~/Downloads/Seurat_scATAC-seq/scenicplus_output/scenicplus/scenicplus_further_analysis"
dir.create("coverage_plots_targets/focused", showWarnings = FALSE, recursive = TRUE)

overlay_colors <- c("Ctrl" = "#c6c6c6", "KO" = "#7FCDCD")

# ── Plot loop ─────────────────────────────────────────────────────────────────
failed  <- c()
success <- c()

for (i in seq_len(nrow(top_targets))) {
  
  row     <- top_targets[i, ]
  gene    <- row$Gene
  ct      <- row$CellType
  tf      <- row$TF
  dir_lbl <- row$direction
  ct_safe <- gsub("[^A-Za-z0-9]", "_", ct)
  tf_safe <- gsub("[^A-Za-z0-9]", "_", tf)
  
  outfile <- file.path("coverage_plots_targets/focused",
                       paste0(gene, "__", tf_safe, "__", ct_safe, ".pdf"))
  
  if (file.exists(outfile)) {
    message("  Skip (exists): ", gene, " [", tf, " / ", ct, "]")
    success <- c(success, gene)
    next
  }
  
  # ── Subset to relevant cell type only ────────────────────────────────────
  obj_ct <- tryCatch(
    subset(obj_male, subset = cell_type == ct),
    error = function(e) NULL
  )
  if (is.null(obj_ct) || ncol(obj_ct) < 20) {
    message("  Skip [", gene, " / ", ct, "] — too few cells")
    next
  }
  
  # group.by genotype directly — gives exactly 2 tracks: Ctrl and KO
  obj_ct$genotype <- factor(obj_ct$genotype, levels = c("Ctrl", "KO"))
  
  tryCatch({
    
    # ── Pass 1: build without ymax to discover signal range ────────────────
    p_raw <- CoveragePlot(
      object            = obj_ct,
      region            = gene,
      assay             = "peaks",
      group.by          = "genotype",
      annotation        = TRUE,
      peaks             = TRUE,
      links             = TRUE,
      extend.upstream   = 10000,
      extend.downstream = 10000
    )
    
    # ── Extract shared ymax across all panels and layers ───────────────────
    ymax_vals <- c()
    for (j in seq_len(length(p_raw))) {
      built <- tryCatch(ggplot_build(p_raw[[j]]), error = function(e) NULL)
      if (is.null(built)) next
      for (layer_data in built$data) {
        if (is.data.frame(layer_data) && "y" %in% colnames(layer_data)) {
          vals <- layer_data$y
          vals <- vals[is.finite(vals) & vals > 0]
          if (length(vals) > 0) ymax_vals <- c(ymax_vals, max(vals))
        }
      }
    }
    
    ymax <- if (length(ymax_vals) > 0 && max(ymax_vals) > 0)
      max(ymax_vals) * 1.05
    else NULL
    
    # ── Pass 2: rebuild with shared ymax and correct colors ────────────────
    p <- CoveragePlot(
      object            = obj_ct,
      region            = gene,
      assay             = "peaks",
      group.by          = "genotype",
      annotation        = TRUE,
      peaks             = TRUE,
      links             = TRUE,
      extend.upstream   = 10000,
      extend.downstream = 10000,
      ymax              = ymax
    ) & scale_fill_manual(values = overlay_colors)
    
    ymax_label <- if (!is.null(ymax)) round(ymax, 1) else "auto"
    
    p_final <- patchwork::wrap_elements(p) +
      patchwork::plot_annotation(
        title    = paste0(gene, "  |  Regulated by: ", tf, "  |  ", dir_lbl),
        subtitle = paste0("Cell type: ", ct,
                          "  |  KO vs Ctrl (male)  |  \u00b110kb  |  shared ymax = ",
                          ymax_label),
        theme    = theme(
          plot.title    = element_text(face = "bold", size = 12),
          plot.subtitle = element_text(size = 9, color = "grey40")
        )
      )
    
    ggsave(outfile, p_final, width = 10, height = 7)
    success <- c(success, gene)
    message("  Saved: ", gene, " [", tf, " \u2192 ", ct, ", ", dir_lbl,
            ", ymax=", ymax_label, "]")
    
  }, error = function(e) {
    failed <<- c(failed, gene)
    message("  FAILED [", gene, "]: ", e$message)
  })
}

message("\n\u2500\u2500 Step 8B summary \u2500\u2500")
message("  Saved:  ", length(success))
message("  Failed: ", length(failed), " \u2014 ",
        if (length(failed) > 0) paste(failed, collapse = ", ") else "none")
message("\n=== Step 8 complete. ===")

########## Overlaid ------
# ══════════════════════════════════════════════════════════════════════════════
# STEP 8B: OVERLAID COVERAGE — Ctrl vs KO same panel
# ══════════════════════════════════════════════════════════════════════════════

DefaultAssay(obj_male) <- "peaks"

further_dir <- "~/Downloads/Seurat_scATAC-seq/scenicplus_output/scenicplus/scenicplus_further_analysis"
dir.create("coverage_plots_targets/overlaid", showWarnings = FALSE, recursive = TRUE)

overlay_colors <- c("Ctrl" = "#c6c6c6", "KO" = "#7FCDCD")
overlay_alpha   <- c("Ctrl" = 0.45,      "KO" = 0.75)     # KO more opaque
overlay_linecolor <- c("Ctrl" = NA,      "KO" = "#1a8a8a") # KO outline only

# ── Helper: extract x,y from a single-group CoveragePlot ─────────────────────
get_cov_df <- function(obj_sub, gene, label, extend_up, extend_dn) {
  p <- CoveragePlot(
    object            = obj_sub,
    region            = gene,
    assay             = "peaks",
    group.by          = "genotype",
    annotation        = FALSE,
    peaks             = FALSE,
    links             = FALSE,
    extend.upstream   = extend_up,
    extend.downstream = extend_dn
  )
  built <- tryCatch(ggplot_build(p), error = function(e) NULL)
  if (is.null(built)) return(NULL)
  df <- tryCatch({
    d <- built$data[[1]]
    if (nrow(d) > 0 && all(c("x", "y") %in% colnames(d))) d else NULL
  }, error = function(e) NULL)
  if (is.null(df)) return(NULL)
  data.frame(x = df$x, y = df$y, genotype = label, stringsAsFactors = FALSE)
}

# ── Helper: clean NA coords from annotation before plotting ──────────────────
clean_annotation <- function(obj) {
  annot <- Annotation(obj)
  if (!is.null(annot) && length(annot) > 0) {
    keep <- !is.na(start(annot)) & !is.na(end(annot))
    Annotation(obj) <- annot[keep]
  }
  obj
}

# ── Plot loop ─────────────────────────────────────────────────────────────────
failed  <- c()
success <- c()

for (i in seq_len(nrow(top_targets))) {
  
  row     <- top_targets[i, ]
  gene    <- row$Gene
  ct      <- row$CellType
  tf      <- row$TF
  dir_lbl <- row$direction
  ct_safe <- gsub("[^A-Za-z0-9]", "_", ct)
  tf_safe <- gsub("[^A-Za-z0-9]", "_", tf)
  
  outfile <- file.path("coverage_plots_targets/overlaid",
                       paste0(gene, "__", tf_safe, "__", ct_safe, ".pdf"))
  
  if (file.exists(outfile)) {
    message("  Skip (exists): ", gene, " [", tf, " / ", ct, "]")
    success <- c(success, gene)
    next
  }
  
  obj_ct <- tryCatch(
    subset(obj_male, subset = cell_type == ct),
    error = function(e) NULL
  )
  if (is.null(obj_ct) || ncol(obj_ct) < 20) {
    message("  Skip [", gene, " / ", ct, "] — too few cells")
    next
  }
  
  obj_ct$genotype <- factor(obj_ct$genotype, levels = c("Ctrl", "KO"))
  
  # Clean NA coords from annotation once per cell type
  obj_ct <- clean_annotation(obj_ct)
  
  tryCatch({
    
    # ── Subset per genotype and extract coverage data ─────────────────────────
    obj_ctrl <- clean_annotation(subset(obj_ct, subset = genotype == "Ctrl"))
    obj_ko   <- clean_annotation(subset(obj_ct, subset = genotype == "KO"))
    
    obj_ctrl$genotype <- factor("Ctrl", levels = "Ctrl")
    obj_ko$genotype   <- factor("KO",   levels = "KO")
    
    df_ctrl <- get_cov_df(obj_ctrl, gene, "Ctrl", 10000, 10000)
    df_ko   <- get_cov_df(obj_ko,   gene, "KO",   10000, 10000)
    
    if (is.null(df_ctrl) && is.null(df_ko)) stop("Could not extract coverage data")
    
    ref_x <- if (!is.null(df_ctrl)) df_ctrl$x else df_ko$x
    if (is.null(df_ctrl)) df_ctrl <- data.frame(x = ref_x, y = 0, genotype = "Ctrl")
    if (is.null(df_ko))   df_ko   <- data.frame(x = ref_x, y = 0, genotype = "KO")
    
    cov_df <- bind_rows(df_ctrl, df_ko) |>
      mutate(genotype = factor(genotype, levels = c("Ctrl", "KO")))
    
    ymax <- max(cov_df$y, na.rm = TRUE) * 1.05
    if (!is.finite(ymax) || ymax == 0) ymax <- 1
    
    # ── Overlaid coverage panel ───────────────────────────────────────────────
    p_cov <- ggplot() +
      # KO: filled area as base
      geom_area(data = dplyr::filter(cov_df, genotype == "KO"),
                aes(x = x, y = y),
                fill = "#2ab3b3", alpha = 0.55, color = NA) +
      # Ctrl: line only — always visible on top regardless of amplitude
      geom_line(data = dplyr::filter(cov_df, genotype == "Ctrl"),
                aes(x = x, y = y),
                color = "#e05c00", linewidth = 0.5, alpha = 0.9) +
      scale_x_continuous(expand = expansion(mult = c(0, 0))) +
      scale_y_continuous(limits = c(0, ymax),
                         expand = expansion(mult = c(0, 0.02))) +
      labs(y = "Coverage", x = NULL,
           title    = paste0(gene, "  |  Regulated by: ", tf, "  |  ", dir_lbl),
           subtitle = paste0("Cell type: ", ct,
                             "  |  KO vs Ctrl (male)  |  \u00b110kb")) +
      theme_minimal(base_size = 11) +
      theme(
        plot.title      = element_text(face = "bold", size = 12),
        plot.subtitle   = element_text(size = 9, color = "grey40"),
        panel.grid      = element_blank(),
        axis.line       = element_line(color = "grey70"),
        axis.ticks      = element_line(color = "grey70"),
        axis.text.x     = element_blank(),
        axis.ticks.x    = element_blank(),
        legend.position = "none"
      ) +
      # Inline labels with colored swatches using annotate + point trick
      annotate("point", x = -Inf, y = ymax * 0.94, shape = 15,
               size = 3, color = "#2ab3b3") +
      annotate("text",  x = -Inf, y = ymax * 0.94,
               label = " KO (filled)", hjust = -0.15, size = 3.5,
               color = "#1a8a8a", fontface = "bold") +
      annotate("segment", x = -Inf, xend = -Inf, y = ymax * 0.82, yend = ymax * 0.82,
               color = "#e05c00", linewidth = 1.2) +
      annotate("text",  x = -Inf, y = ymax * 0.82,
               label = " Ctrl (line)", hjust = -0.15, size = 3.5,
               color = "#e05c00", fontface = "bold")
    
    # ── Peaks panel (tryCatch in case no peaks in region) ────────────────────
    p_peaks <- tryCatch(
      PeakPlot(
        object            = obj_ct,
        region            = gene,
        assay             = "peaks",
        extend.upstream   = 10000,
        extend.downstream = 10000
      ) + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()),
      error = function(e) NULL
    )
    
    # ── Annotation panel (tryCatch in case gene not in annotation) ────────────
    p_annot <- tryCatch(
      AnnotationPlot(
        object            = obj_ct,
        region            = gene,
        extend.upstream   = 10000,
        extend.downstream = 10000
      ),
      error = function(e) NULL
    )
    
    # ── Combine — adapt layout to which panels succeeded ─────────────────────
    if (!is.null(p_peaks) && !is.null(p_annot)) {
      p_final <- p_cov / p_peaks / p_annot +
        plot_layout(heights = c(4, 0.5, 1))
    } else if (!is.null(p_annot)) {
      p_final <- p_cov / p_annot +
        plot_layout(heights = c(4, 1))
    } else if (!is.null(p_peaks)) {
      p_final <- p_cov / p_peaks +
        plot_layout(heights = c(4, 0.5))
    } else {
      p_final <- p_cov
    }
    
    ggsave(outfile, p_final, width = 10, height = 5)
    success <- c(success, gene)
    message("  Saved: ", gene, " [", tf, " \u2192 ", ct, ", ", dir_lbl, "]")
    
  }, error = function(e) {
    failed <<- c(failed, gene)
    message("  FAILED [", gene, "]: ", e$message)
  })
}

message("\n\u2500\u2500 Step 8B summary \u2500\u2500")
message("  Saved:  ", length(success))
message("  Failed: ", length(failed), " \u2014 ",
        if (length(failed) > 0) paste(failed, collapse = ", ") else "none")
message("\n=== Step 8B complete. ===")

########## Mirrored #1 ------
# ══════════════════════════════════════════════════════════════════════════════
# STEP 8C: MIRRORED COVERAGE — Ctrl up / KO down
# ══════════════════════════════════════════════════════════════════════════════

DefaultAssay(obj_male) <- "peaks"

further_dir <- "~/Downloads/Seurat_scATAC-seq/scenicplus_output/scenicplus/scenicplus_further_analysis"
dir.create("coverage_plots_targets/mirrored", showWarnings = FALSE, recursive = TRUE)

overlay_colors <- c("Ctrl" = "#c6c6c6", "KO" = "#7FCDCD")

failed  <- c()
success <- c()

for (i in seq_len(nrow(top_targets))) {
  
  row     <- top_targets[i, ]
  gene    <- row$Gene
  ct      <- row$CellType
  tf      <- row$TF
  dir_lbl <- row$direction
  ct_safe <- gsub("[^A-Za-z0-9]", "_", ct)
  tf_safe <- gsub("[^A-Za-z0-9]", "_", tf)
  
  outfile <- file.path("coverage_plots_targets/mirrored",
                       paste0(gene, "__", tf_safe, "__", ct_safe, ".pdf"))
  
  if (file.exists(outfile)) {
    message("  Skip (exists): ", gene, " [", tf, " / ", ct, "]")
    success <- c(success, gene)
    next
  }
  
  obj_ct <- tryCatch(
    subset(obj_male, subset = cell_type == ct),
    error = function(e) NULL
  )
  if (is.null(obj_ct) || ncol(obj_ct) < 20) {
    message("  Skip [", gene, " / ", ct, "] — too few cells")
    next
  }
  
  obj_ct$genotype <- factor(obj_ct$genotype, levels = c("Ctrl", "KO"))
  obj_ct <- clean_annotation(obj_ct)
  
  tryCatch({
    
    obj_ctrl <- clean_annotation(subset(obj_ct, subset = genotype == "Ctrl"))
    obj_ko   <- clean_annotation(subset(obj_ct, subset = genotype == "KO"))
    obj_ctrl$genotype <- factor("Ctrl", levels = "Ctrl")
    obj_ko$genotype   <- factor("KO",   levels = "KO")
    
    df_ctrl <- get_cov_df(obj_ctrl, gene, "Ctrl", 10000, 10000)
    df_ko   <- get_cov_df(obj_ko,   gene, "KO",   10000, 10000)
    
    if (is.null(df_ctrl) && is.null(df_ko)) stop("Could not extract coverage data")
    
    ref_x <- if (!is.null(df_ctrl)) df_ctrl$x else df_ko$x
    if (is.null(df_ctrl)) df_ctrl <- data.frame(x = ref_x, y = 0, genotype = "Ctrl")
    if (is.null(df_ko))   df_ko   <- data.frame(x = ref_x, y = 0, genotype = "KO")
    
    # ── Mirror: Ctrl positive, KO negative ───────────────────────────────────
    df_ctrl_plot <- df_ctrl |> mutate(y_plot =  y, genotype = "Ctrl")
    df_ko_plot   <- df_ko   |> mutate(y_plot = -y, genotype = "KO")
    
    cov_df <- bind_rows(df_ctrl_plot, df_ko_plot) |>
      mutate(genotype = factor(genotype, levels = c("Ctrl", "KO")))
    
    ymax <- max(df_ctrl$y, na.rm = TRUE) * 1.05
    ymin <- -max(df_ko$y,  na.rm = TRUE) * 1.05
    if (!is.finite(ymax) || ymax == 0) ymax <- 1
    if (!is.finite(ymin) || ymin == 0) ymin <- -1
    
    y_breaks <- pretty(c(ymin, ymax), n = 4)
    y_labels <- as.character(abs(y_breaks))
    
    p_cov <- ggplot(cov_df, aes(x = x, y = y_plot, fill = genotype)) +
      geom_area(alpha = 0.7, position = "identity") +
      geom_hline(yintercept = 0, color = "grey40", linewidth = 0.4) +
      scale_fill_manual(values = overlay_colors, name = NULL) +
      scale_x_continuous(expand = expansion(mult = c(0, 0))) +
      scale_y_continuous(
        limits = c(ymin, ymax),
        breaks = y_breaks,
        labels = y_labels,
        expand = expansion(mult = c(0.02, 0.02))
      ) +
      annotate("text", x = -Inf, y =  ymax * 0.85, label = "Ctrl",
               hjust = -0.2, size = 3.5, color = "#999999", fontface = "bold") +
      annotate("text", x = -Inf, y =  ymin * 0.85, label = "KO",
               hjust = -0.2, size = 3.5, color = "#7FCDCD", fontface = "bold") +
      labs(y = "Coverage", x = NULL,
           title    = paste0(gene, "  |  Regulated by: ", tf, "  |  ", dir_lbl),
           subtitle = paste0("Cell type: ", ct,
                             "  |  KO \u25bc vs Ctrl \u25b2 (male)  |  \u00b110kb")) +
      theme_minimal(base_size = 11) +
      theme(
        plot.title           = element_text(face = "bold", size = 12),
        plot.subtitle        = element_text(size = 9, color = "grey40"),
        panel.grid           = element_blank(),
        axis.line.y          = element_line(color = "grey70"),
        axis.ticks.y         = element_line(color = "grey70"),
        axis.text.x          = element_blank(),
        axis.ticks.x         = element_blank(),
        legend.position      = "none"
      )
    
    # ── Peaks + annotation ────────────────────────────────────────────────────
    # ── Peaks panel ───────────────────────────────────────────────────────────
    p_peaks <- tryCatch(
      PeakPlot(
        object            = obj_ct,
        region            = gene,
        assay             = "peaks",
        extend.upstream   = 10000,
        extend.downstream = 10000
      ) + theme(
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin  = margin(2, 2, 2, 2)
      ),
      error = function(e) NULL
    )
    
    # ── Annotation panel ──────────────────────────────────────────────────────
    p_annot <- tryCatch(
      AnnotationPlot(
        object            = obj_ct,
        region            = gene,
        extend.upstream   = 10000,
        extend.downstream = 10000
      ) + theme(plot.margin = margin(4, 2, 4, 2)),
      error = function(e) NULL
    )
    
    # ── Measure annotation panel's natural height by rendering to grob ────────
    get_natural_height_in <- function(p, width_in = 10) {
      tryCatch({
        gt   <- ggplot_gtable(ggplot_build(p))
        # grob heights are in absolute units where possible; sum them
        h    <- sum(grid::convertHeight(gt$heights, "in", valueOnly = TRUE),
                    na.rm = TRUE)
        # clamp: minimum 0.8in, maximum 3in
        min(3, max(0.8, h))
      }, error = function(e) 1.5)
    }
    
    annot_h_in <- if (!is.null(p_annot)) get_natural_height_in(p_annot) else 0
    peaks_h_in <- 0.4
    
    # ── Combine with exact inch heights ───────────────────────────────────────
    cov_h_in <- 4
    
    if (!is.null(p_peaks) && !is.null(p_annot)) {
      p_final <- p_cov / p_peaks / p_annot +
        plot_layout(heights = unit(c(cov_h_in, peaks_h_in, annot_h_in), "in"))
      total_h <- cov_h_in + peaks_h_in + annot_h_in
    } else if (!is.null(p_annot)) {
      p_final <- p_cov / p_annot +
        plot_layout(heights = unit(c(cov_h_in, annot_h_in), "in"))
      total_h <- cov_h_in + annot_h_in
    } else if (!is.null(p_peaks)) {
      p_final <- p_cov / p_peaks +
        plot_layout(heights = unit(c(cov_h_in, peaks_h_in), "in"))
      total_h <- cov_h_in + peaks_h_in
    } else {
      p_final <- p_cov
      total_h <- cov_h_in
    }
    
    # Add small padding for title/subtitle at top
    total_h <- total_h + 0.6
    
    ggsave(outfile, p_final, width = 10, height = total_h, limitsize = FALSE)
    success <- c(success, gene)
    message("  Saved: ", gene, " [", tf, " \u2192 ", ct, ", ", dir_lbl, "]")
    
  }, error = function(e) {
    failed <<- c(failed, gene)
    message("  FAILED [", gene, "]: ", e$message)
  })
}

message("\n\u2500\u2500 Step 8C summary \u2500\u2500")
message("  Saved:  ", length(success))
message("  Failed: ", length(failed), " \u2014 ",
        if (length(failed) > 0) paste(failed, collapse = ", ") else "none")
message("\n=== Step 8C complete. ===")


########## Mirrored #2 ------
# ══════════════════════════════════════════════════════════════════════════════
# STEP 8C: MIRRORED COVERAGE — Ctrl up / KO down
# ══════════════════════════════════════════════════════════════════════════════
DefaultAssay(obj_male) <- "peaks"

further_dir <- "~/Downloads/Seurat_scATAC-seq/scenicplus_output/scenicplus/scenicplus_further_analysis"
dir.create("coverage_plots_targets/mirrored", showWarnings = FALSE, recursive = TRUE)

overlay_colors <- c("Ctrl" = "#c6c6c6", "KO" = "#7FCDCD")

failed  <- c()
success <- c()

# ── Helper: build gene track manually with fixed row height ──────────────────
make_gene_track <- function(obj_ct, gene, extend_up, extend_dn,
                            row_height_in = 0.35, label_size = 3) {
  
  annot <- Annotation(obj_ct)
  if (is.null(annot) || length(annot) == 0) return(list(plot = NULL, height_in = 0))
  
  target <- annot[!is.na(annot$gene_name) & annot$gene_name == gene]
  if (length(target) == 0) return(list(plot = NULL, height_in = 0))
  
  chr  <- as.character(seqnames(target[1]))
  tss  <- ifelse(as.character(strand(target[1])) == "+",
                 min(start(target)), max(end(target)))
  xmin <- tss - extend_up
  xmax <- tss + extend_dn
  
  in_reg <- annot[
    !is.na(start(annot)) & !is.na(end(annot)) &
      as.character(seqnames(annot)) == chr &
      start(annot) <= xmax & end(annot) >= xmin
  ]
  if (length(in_reg) == 0) return(list(plot = NULL, height_in = 0))
  
  df <- as.data.frame(in_reg, row.names = NULL) |>
    dplyr::mutate(
      gene_name = ifelse(is.na(gene_name), "unknown", as.character(gene_name)),
      strand    = as.character(strand),
      start_cl  = pmax(start, xmin),
      end_cl    = pmin(end,   xmax)
    ) |>
    dplyr::filter(start_cl < end_cl)
  
  if (nrow(df) == 0) return(list(plot = NULL, height_in = 0))
  
  # ── ONE row per gene: full extent as backbone ─────────────────────────────
  transcripts <- df |>
    dplyr::group_by(gene_name) |>
    dplyr::summarise(
      start_cl = min(start_cl),
      end_cl   = max(end_cl),
      strand   = dplyr::first(strand),
      .groups  = "drop"
    ) |>
    dplyr::arrange(start_cl)
  
  # ── Union of exons/CDS per gene (merge overlapping ranges) ────────────────
  exon_raw <- df |>
    dplyr::filter(type %in% c("exon", "cds")) |>
    dplyr::select(gene_name, start_cl, end_cl) |>
    dplyr::arrange(gene_name, start_cl)
  
  # Merge overlapping exon intervals per gene
  exons <- exon_raw |>
    dplyr::group_by(gene_name) |>
    dplyr::group_modify(~ {
      d <- dplyr::arrange(.x, start_cl)
      merged <- d[1, , drop = FALSE]
      if (nrow(d) > 1) {
        for (j in 2:nrow(d)) {
          if (d$start_cl[j] <= merged$end_cl[nrow(merged)]) {
            merged$end_cl[nrow(merged)] <- max(merged$end_cl[nrow(merged)], d$end_cl[j])
          } else {
            merged <- dplyr::bind_rows(merged, d[j, , drop = FALSE])
          }
        }
      }
      merged
    }) |>
    dplyr::ungroup()
  
  if (nrow(transcripts) == 0) return(list(plot = NULL, height_in = 0))
  
  # Greedy row packer (one row per gene now, so usually very few rows)
  transcripts$row <- 1L
  row_ends <- c(transcripts$end_cl[1])
  if (nrow(transcripts) > 1) {
    for (k in 2:nrow(transcripts)) {
      placed <- FALSE
      for (r in seq_along(row_ends)) {
        if (transcripts$start_cl[k] > row_ends[r] + 500) {
          transcripts$row[k] <- as.integer(r)
          row_ends[r] <- transcripts$end_cl[k]
          placed <- TRUE; break
        }
      }
      if (!placed) {
        row_ends <- c(row_ends, transcripts$end_cl[k])
        transcripts$row[k] <- as.integer(length(row_ends))
      }
    }
  }
  
  n_rows <- max(transcripts$row)
  h_in   <- max(0.8, n_rows * row_height_in + 0.7)
  
  row_lookup <- transcripts |>
    dplyr::select(gene_name, row) |>
    dplyr::distinct(gene_name, .keep_all = TRUE)
  
  exons <- exons |>
    dplyr::left_join(row_lookup, by = "gene_name") |>
    dplyr::filter(!is.na(row))
  
  # Direction arrow ticks along backbone
  arrow_spacing <- (xmax - xmin) / 20
  arrows_df <- transcripts |>
    dplyr::rowwise() |>
    dplyr::mutate(
      tick_pos = list({
        pos <- seq(start_cl + arrow_spacing / 2,
                   end_cl  - arrow_spacing / 2,
                   by = arrow_spacing)
        if (length(pos) == 0) pos <- (start_cl + end_cl) / 2
        pos
      })
    ) |>
    tidyr::unnest(tick_pos) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      x_from = ifelse(strand == "+", tick_pos - arrow_spacing * 0.3,
                      tick_pos + arrow_spacing * 0.3),
      x_to   = ifelse(strand == "+", tick_pos + arrow_spacing * 0.3,
                      tick_pos - arrow_spacing * 0.3)
    )
  
  all_genes   <- unique(transcripts$gene_name)
  gene_colors <- setNames(
    ifelse(all_genes == gene, "#1a6e1a", "#3333aa"),
    all_genes
  )
  
  p <- ggplot() +
    geom_segment(data = transcripts,
                 aes(x = start_cl, xend = end_cl, y = row, yend = row,
                     color = gene_name),
                 linewidth = 0.5) +
    geom_segment(data = arrows_df,
                 aes(x = x_from, xend = x_to, y = row, yend = row,
                     color = gene_name),
                 linewidth = 0.4,
                 arrow = arrow(length = unit(0.12, "cm"), type = "open",
                               ends = "last")) +
    geom_rect(data = exons,
              aes(xmin = start_cl, xmax = end_cl,
                  ymin = row - 0.3, ymax = row + 0.3,
                  fill = gene_name),
              color = NA) +
    geom_text(data = transcripts,
              aes(x = (start_cl + end_cl) / 2, y = row + 0.44,
                  label = gene_name, color = gene_name),
              size = label_size, hjust = 0.5, vjust = 0, fontface = "italic") +
    scale_color_manual(values = gene_colors, guide = "none") +
    scale_fill_manual(values  = gene_colors, guide = "none") +
    scale_x_continuous(
      limits = c(xmin, xmax),
      expand = expansion(mult = c(0, 0)),
      labels = scales::label_number(scale = 1e-6, suffix = " Mb", accuracy = 0.01)
    ) +
    scale_y_continuous(limits = c(0.3, n_rows + 0.9), expand = c(0, 0)) +
    labs(x = NULL, y = "Genes") +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid   = element_blank(),
      axis.line.x  = element_line(color = "grey70"),
      axis.ticks.x = element_line(color = "grey70"),
      axis.text.x  = element_text(size = 8),
      axis.title.y = element_text(size = 9, angle = 90),
      axis.text.y  = element_blank(),
      axis.ticks.y = element_blank(),
      plot.margin  = margin(2, 2, 6, 2)
    )
  
  list(plot = p, height_in = h_in)
}

for (i in seq_len(nrow(top_targets))) {
  
  row     <- top_targets[i, ]
  gene    <- row$Gene
  ct      <- row$CellType
  tf      <- row$TF
  dir_lbl <- row$direction
  ct_safe <- gsub("[^A-Za-z0-9]", "_", ct)
  tf_safe <- gsub("[^A-Za-z0-9]", "_", tf)
  
  outfile <- file.path("coverage_plots_targets/mirrored",
                       paste0(gene, "__", tf_safe, "__", ct_safe, ".pdf"))
  
  if (file.exists(outfile)) {
    message("  Skip (exists): ", gene, " [", tf, " / ", ct, "]")
    success <- c(success, gene)
    next
  }
  
  obj_ct <- tryCatch(
    subset(obj_male, subset = cell_type == ct),
    error = function(e) NULL
  )
  if (is.null(obj_ct) || ncol(obj_ct) < 20) {
    message("  Skip [", gene, " / ", ct, "] — too few cells")
    next
  }
  
  obj_ct$genotype <- factor(obj_ct$genotype, levels = c("Ctrl", "KO"))
  obj_ct <- clean_annotation(obj_ct)
  
  tryCatch({
    
    obj_ctrl <- clean_annotation(subset(obj_ct, subset = genotype == "Ctrl"))
    obj_ko   <- clean_annotation(subset(obj_ct, subset = genotype == "KO"))
    obj_ctrl$genotype <- factor("Ctrl", levels = "Ctrl")
    obj_ko$genotype   <- factor("KO",   levels = "KO")
    
    df_ctrl <- get_cov_df(obj_ctrl, gene, "Ctrl", 10000, 10000)
    df_ko   <- get_cov_df(obj_ko,   gene, "KO",   10000, 10000)
    
    if (is.null(df_ctrl) && is.null(df_ko)) stop("Could not extract coverage data")
    
    ref_x <- if (!is.null(df_ctrl)) df_ctrl$x else df_ko$x
    if (is.null(df_ctrl)) df_ctrl <- data.frame(x = ref_x, y = 0, genotype = "Ctrl")
    if (is.null(df_ko))   df_ko   <- data.frame(x = ref_x, y = 0, genotype = "KO")
    
    # ── Mirror: Ctrl positive, KO negative ───────────────────────────────────
    df_ctrl_plot <- df_ctrl |> mutate(y_plot =  y, genotype = "Ctrl")
    df_ko_plot   <- df_ko   |> mutate(y_plot = -y, genotype = "KO")
    
    cov_df <- bind_rows(df_ctrl_plot, df_ko_plot) |>
      mutate(genotype = factor(genotype, levels = c("Ctrl", "KO")))
    
    ymax <- max(df_ctrl$y, na.rm = TRUE) * 1.05
    ymin <- -max(df_ko$y,  na.rm = TRUE) * 1.05
    if (!is.finite(ymax) || ymax == 0) ymax <- 1
    if (!is.finite(ymin) || ymin == 0) ymin <- -1
    
    y_breaks <- pretty(c(ymin, ymax), n = 4)
    y_labels <- as.character(abs(y_breaks))
    
    p_cov <- ggplot(cov_df, aes(x = x, y = y_plot, fill = genotype)) +
      geom_area(alpha = 0.7, position = "identity") +
      geom_hline(yintercept = 0, color = "grey40", linewidth = 0.4) +
      scale_fill_manual(values = overlay_colors, name = NULL) +
      scale_x_continuous(expand = expansion(mult = c(0, 0))) +
      scale_y_continuous(
        limits = c(ymin, ymax),
        breaks = y_breaks,
        labels = y_labels,
        expand = expansion(mult = c(0.02, 0.02))
      ) +
      annotate("text", x = -Inf, y =  ymax * 0.85, label = "Ctrl",
               hjust = -0.2, size = 3.5, color = "#999999", fontface = "bold") +
      annotate("text", x = -Inf, y =  ymin * 0.85, label = "KO",
               hjust = -0.2, size = 3.5, color = "#7FCDCD", fontface = "bold") +
      labs(y = "Coverage", x = NULL,
           title    = paste0(gene, "  |  Regulated by: ", tf, "  |  ", dir_lbl),
           subtitle = paste0("Cell type: ", ct,
                             "  |  KO \u25bc vs Ctrl \u25b2 (male)  |  \u00b110kb")) +
      theme_minimal(base_size = 11) +
      theme(
        plot.title           = element_text(face = "bold", size = 12),
        plot.subtitle        = element_text(size = 9, color = "grey40"),
        panel.grid           = element_blank(),
        axis.line.y          = element_line(color = "grey70"),
        axis.ticks.y         = element_line(color = "grey70"),
        axis.text.x          = element_blank(),
        axis.ticks.x         = element_blank(),
        legend.position      = "none"
      )
    
    # ── Peaks + annotation ────────────────────────────────────────────────────
    # ── Peaks panel ───────────────────────────────────────────────────────────
    p_peaks <- tryCatch(
      PeakPlot(
        object            = obj_ct,
        region            = gene,
        assay             = "peaks",
        extend.upstream   = 10000,
        extend.downstream = 10000
      ) + theme(
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin  = margin(2, 2, 2, 2)
      ),
      error = function(e) NULL
    )
    
    # ── Annotation: manual gene track ─────────────────────────────────────────
    gene_track <- make_gene_track(obj_ct, gene, 10000, 10000,
                                  row_height_in = 0.35, label_size = 3)
    p_annot  <- gene_track$plot
    annot_h  <- gene_track$height_in   # exact inches this track needs
    
    peaks_h  <- 0.5
    cov_h    <- 4.0
    title_h  <- 0.8   # increased — accounts for title + subtitle + patchwork top padding
    
    p_cov <- p_cov + theme(plot.margin = margin(2, 2, 2, 2))
    
    if (!is.null(p_peaks) && !is.null(p_annot)) {
      p_final <- p_cov / p_peaks / p_annot +
        plot_layout(heights = unit(c(cov_h, peaks_h, annot_h), "in"))
      total_h <- title_h + cov_h + peaks_h + annot_h + 0.5   # ← +0.5 buffer
    } else if (!is.null(p_annot)) {
      p_final <- p_cov / p_annot +
        plot_layout(heights = unit(c(cov_h, annot_h), "in"))
      total_h <- title_h + cov_h + annot_h + 0.5
    } else if (!is.null(p_peaks)) {
      p_final <- p_cov / p_peaks +
        plot_layout(heights = unit(c(cov_h, peaks_h), "in"))
      total_h <- title_h + cov_h + peaks_h + 0.5
    } else {
      p_final <- p_cov
      total_h <- title_h + cov_h + 0.5
    }
    
    ggsave(outfile, p_final, width = 10, height = total_h, limitsize = FALSE)
    success <- c(success, gene)
    message("  Saved: ", gene, " [", tf, " \u2192 ", ct, ", ", dir_lbl, "]")
    
  }, error = function(e) {
    failed <<- c(failed, gene)
    message("  FAILED [", gene, "]: ", e$message)
  })
}

message("\n\u2500\u2500 Step 8C summary \u2500\u2500")
message("  Saved:  ", length(success))
message("  Failed: ", length(failed), " \u2014 ",
        if (length(failed) > 0) paste(failed, collapse = ", ") else "none")
message("\n=== Step 8C complete. ===")
