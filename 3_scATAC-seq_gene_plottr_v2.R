# ══════════════════════════════════════════════════════════════════════════════
# plot_gene_differential()
# Interactive wrapper: give a gene name (+ optional cell type), get the
# differential ATAC coverage plot with peak track, gene model, and cCREs.
#
# Prerequisites (must already be loaded in session):
#   obj_male, ccre_gr, diff_colors, atac_colors
#   helpers: resolve_gene_alias, gene_to_region, get_cov_df,
#            clean_annotation, make_gene_track
# ══════════════════════════════════════════════════════════════════════════════

###### Set up -------

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

# Reload if starting fresh
if (!exists("double_hit_df")) {
  double_hit_df <- read.csv("TF_double_hit.csv")
}
if (!exists("overlap_df")) {
  overlap_df <- read.csv("DA_DEG_overlap.csv")
  if (nrow(overlap_df) == 0) overlap_df <- data.frame()
}

obj_male <- readRDS("scATAC_07_male_integrated.rds")

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

DefaultAssay(obj_male) <- "peaks"

diff_colors <- c("UP" = "#b2182b", "DOWN" = "#2166ac")

###### Helper functions ------------

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

###### Make gene tracks -------

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
    chr        <- as.character(seqnames(target)[1])
    gene_start <- min(start(target))
    gene_end   <- max(end(target))
    gene_mid   <- (gene_start + gene_end) / 2
    xmin       <- gene_mid - extend_up
    xmax       <- gene_mid + extend_dn
  }
  
  if (is.null(annot) || length(annot) == 0) {
    in_reg <- GenomicRanges::GRanges()
  } else {
    in_reg <- annot[
      !is.na(start(annot)) & !is.na(end(annot)) &
        as.character(seqnames(annot)) == chr &
        start(annot) <= xmax & end(annot) >= xmin
    ]
  }
  
  if (length(in_reg) == 0) {
    transcripts <- data.frame(gene_name = character(), start_cl = numeric(),
                              end_cl = numeric(), strand = character(),
                              row = integer(), is_target = logical(),
                              stringsAsFactors = FALSE)
    exons     <- data.frame(gene_name = character(), start_cl = numeric(),
                            end_cl = numeric(), row = integer(),
                            fill_id = character(), stringsAsFactors = FALSE)
    arrows_df <- data.frame(gene_name = character(), row = integer(),
                            strand = character(), x_from = numeric(),
                            x_to = numeric(), stringsAsFactors = FALSE)
    n_rows <- 1L
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
      exons     <- data.frame(gene_name = character(), start_cl = numeric(),
                              end_cl = numeric(), row = integer(),
                              fill_id = character(), stringsAsFactors = FALSE)
      arrows_df <- data.frame(gene_name = character(), row = integer(),
                              strand = character(), x_from = numeric(),
                              x_to = numeric(), stringsAsFactors = FALSE)
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
      
      n_rows     <- max(transcripts$row)
      row_lookup <- dplyr::distinct(dplyr::select(transcripts, gene_name, row),
                                    gene_name, .keep_all = TRUE)
      exons <- dplyr::left_join(exons, row_lookup, by = "gene_name") |>
        dplyr::filter(!is.na(row))
      
      # Arrow ticks
      arrow_spacing <- (xmax - xmin) / 15
      
      arrows_df <- do.call(dplyr::bind_rows, lapply(seq_len(nrow(transcripts)), function(k) {
        tx        <- transcripts[k, ]
        ex        <- exons[exons$gene_name == tx$gene_name, , drop = FALSE]
        sp        <- arrow_spacing
        gene_span <- tx$end_cl - tx$start_cl
        
        # Gene narrower than one spacing — single arrow at midpoint
        if (gene_span <= sp) {
          mid <- (tx$start_cl + tx$end_cl) / 2
          return(data.frame(
            gene_name = tx$gene_name, row = tx$row, strand = tx$strand,
            x_from = ifelse(tx$strand == "+", mid - sp * 0.35, mid + sp * 0.35),
            x_to   = ifelse(tx$strand == "+", mid + sp * 0.35, mid - sp * 0.35),
            stringsAsFactors = FALSE
          ))
        }
        
        pos <- seq(tx$start_cl + sp / 2, tx$end_cl - sp / 2, by = sp)
        if (length(pos) == 0) pos <- (tx$start_cl + tx$end_cl) / 2
        
        if (nrow(ex) > 0) {
          exon_buf   <- sp * 0.5
          in_exon    <- sapply(pos, function(p)
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
                return(data.frame(gene_name = character(), row = integer(),
                                  strand = character(), x_from = numeric(),
                                  x_to = numeric(), stringsAsFactors = FALSE))
              }
            } else {
              return(data.frame(gene_name = character(), row = integer(),
                                strand = character(), x_from = numeric(),
                                x_to = numeric(), stringsAsFactors = FALSE))
            }
          }
          pos <- pos_intron
        }
        
        data.frame(
          gene_name = tx$gene_name, row = tx$row, strand = tx$strand,
          x_from = ifelse(tx$strand == "+", pos - sp * 0.35, pos + sp * 0.35),
          x_to   = ifelse(tx$strand == "+", pos + sp * 0.35, pos - sp * 0.35),
          stringsAsFactors = FALSE
        )
      }))
      
      transcripts$is_target <- transcripts$gene_name == gene
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
  ccre_colors <- c("PLS" = "#FF0000", "pELS" = "#FFA500", "dELS" = "#FFCD00",
                   "CTCF-only" = "#00B0F0", "DNase-H3K4me3" = "#00B050")
  
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
  
  has_ccre    <- !is.null(ccre_df) && nrow(ccre_df) > 0
  ccre_extra  <- if (has_ccre) 0.55 else 0
  h_in        <- max(0.8, n_rows * row_height_in + 0.7 + ccre_extra)
  y_floor     <- if (has_ccre) ccre_ymin - 0.15 else -0.3
  
  # Single unified fill scale: gene colors + cCRE colors merged
  unified_fills <- c(gene_colors, ccre_colors)
  if (has_ccre) ccre_df$fill_id <- ccre_df$type
  
  # ── Build plot ─────────────────────────────────────────────────────────────
  p <- ggplot()
  
  if (nrow(transcripts) > 0)
    p <- p + geom_segment(data = transcripts,
                          aes(x = start_cl, xend = end_cl, y = row, yend = row,
                              color = gene_name),
                          linewidth = 0.5)
  
  if (nrow(exons) > 0)
    p <- p + geom_rect(data = exons,
                       aes(xmin = start_cl, xmax = end_cl,
                           ymin = row - 0.18, ymax = row + 0.18,
                           fill = fill_id),
                       color = NA)
  
  if (nrow(arrows_df) > 0)
    p <- p + geom_segment(data = arrows_df,
                          aes(x = x_from, xend = x_to, y = row, yend = row,
                              color = gene_name),
                          linewidth = 0.3,
                          arrow = arrow(length = unit(0.14, "cm"),
                                        type = "open", ends = "last"))
  
  # ── Labels: dynamic hjust so edge genes don't get clipped ─────────────────
  win_width <- xmax - xmin
  if (nrow(transcripts) > 0) {
    tx_labels <- transcripts |>
      dplyr::mutate(
        mid     = (start_cl + end_cl) / 2,
        # left-align if gene is in left 20% of window, right-align if in right 20%
        hjust   = dplyr::case_when(
          mid < xmin + win_width * 0.2 ~ 0,
          mid > xmax - win_width * 0.2 ~ 1,
          TRUE ~ 0.5
        ),
        label_x = dplyr::case_when(
          hjust == 0 ~ pmax(mid, xmin),
          hjust == 1 ~ pmin(mid, xmax),
          TRUE ~ mid
        )
      )
    
    target_tx <- dplyr::filter(tx_labels,  is_target)
    other_tx  <- dplyr::filter(tx_labels, !is_target)
    
    # Draw each hjust group separately (geom_text doesn't accept per-row hjust)
    for (hj in unique(other_tx$hjust)) {
      sub <- dplyr::filter(other_tx, hjust == hj)
      p <- p + geom_text(data = sub,
                         aes(x = label_x, y = row + 0.38, label = gene_name),
                         color = "#888888", size = label_size,
                         hjust = hj, vjust = 0, fontface = "italic")
    }
    for (hj in unique(target_tx$hjust)) {
      sub <- dplyr::filter(target_tx, hjust == hj)
      p <- p + geom_text(data = sub,
                         aes(x = label_x, y = row + 0.38, label = gene_name),
                         color = "#1a6e1a", size = label_size,
                         hjust = hj, vjust = 0, fontface = "bold.italic")
    }
  } else {
    p <- p + annotate("text", x = (xmin + xmax) / 2, y = 1, label = gene,
                      color = "#1a6e1a", size = label_size,
                      hjust = 0.5, vjust = 0, fontface = "bold.italic")
  }
  
  p <- p + scale_color_manual(values = gene_colors, guide = "none")
  
  # ── cCRE bars — fill merged into single scale to avoid double-scale conflict
  if (has_ccre) {
    p <- p +
      geom_hline(yintercept = 0, color = "grey80", linewidth = 0.3,
                 linetype = "dashed") +
      geom_rect(data = ccre_df,
                aes(xmin = start_cl, xmax = end_cl,
                    ymin = ccre_ymin, ymax = ccre_ymax, fill = fill_id),
                color = NA, inherit.aes = FALSE) +
      geom_text(data = data.frame(x = xmin, y = (ccre_ymin + ccre_ymax) / 2,
                                  label = "cCREs"),
                aes(x = x, y = y, label = label),
                hjust = 0, size = 2.8, color = "grey40", inherit.aes = FALSE) +
      scale_fill_manual(
        values = unified_fills,
        breaks = names(ccre_colors),   # only cCRE types appear in legend
        labels = c("PLS"           = "Promoter (PLS)",
                   "pELS"          = "Proximal enhancer (pELS)",
                   "dELS"          = "Distal enhancer (dELS)",
                   "CTCF-only"     = "CTCF insulator",
                   "DNase-H3K4me3" = "DNase+H3K4me3"),
        name  = "ENCODE cCRE",
        guide = guide_legend(override.aes = list(size = 3),
                             keyheight = unit(0.35, "cm"),
                             keywidth  = unit(0.35, "cm")))
  } else {
    p <- p + scale_fill_manual(values = unified_fills, guide = "none")
  }
  
  p <- p +
    scale_x_continuous(limits  = c(xmin, xmax),
                       expand  = expansion(mult = c(0, 0)),
                       labels  = scales::label_number(scale = 1e-6, suffix = " Mb",
                                                      accuracy = 0.01)) +
    scale_y_continuous(limits = c(y_floor, n_rows + 0.9), expand = c(0, 0)) +
    coord_cartesian(clip = "off") +
    labs(x = NULL, y = "Genes") +
    theme_minimal(base_size = 10) +
    theme(panel.grid      = element_blank(),
          axis.line.x     = element_line(color = "grey70"),
          axis.ticks.x    = element_line(color = "grey70"),
          axis.text.x     = element_text(size = 8),
          axis.title.y    = element_text(size = 9, angle = 90),
          axis.text.y     = element_blank(),
          axis.ticks.y    = element_blank(),
          legend.position = if (has_ccre) "right" else "none",
          legend.key.size = unit(0.35, "cm"),
          legend.text     = element_text(size = 7),
          legend.title    = element_text(size = 8),
          plot.margin     = margin(2, 40, 6, 40))
  
  list(plot = p, height_in = h_in)
}


###### Plot differential --------

plot_gene_differential <- function(
    gene,
    cell_type  = NULL,
    contrast   = "KOvCtrl_male",
    atac_obj   = obj_male,
    deg_dir    = "/Users/cnbr/Downloads/Seurat_scATAC-seq/DEG_results",
    save_pdf   = FALSE,
    out_dir    = "/Users/cnbr/Downloads/Seurat_scATAC-seq/coverage_plots_targets/interactive",
    extend_up  = 10000,
    extend_dn  = 10000
) {
  
  # ── 0. Resolve expression direction from DEG files ───────────────────────
  get_expr_direction <- function(gene, cell_type, contrast, deg_dir) {
    normalise <- function(x) tolower(gsub("[^A-Za-z0-9]", "", x))
    
    # List all DEG files, then filter — more robust than a single regex
    all_deg_files <- list.files(deg_dir, pattern = "DEG\\.csv$",
                                full.names = TRUE, ignore.case = TRUE)
    if (length(all_deg_files) == 0) return(NULL)
    
    # Filter by contrast (exact string in filename, not normalised)
    all_deg_files <- all_deg_files[grepl(contrast, basename(all_deg_files),
                                         fixed = TRUE)]
    
    # Filter by cell type if provided (normalised comparison)
    if (!is.null(cell_type)) {
      ct_norm    <- normalise(cell_type)
      file_norms <- normalise(basename(all_deg_files))
      all_deg_files <- all_deg_files[grepl(ct_norm, file_norms, fixed = TRUE)]
    }
    
    if (length(all_deg_files) == 0) return(NULL)
    
    for (f in all_deg_files) {
      deg <- tryCatch(read.csv(f, row.names = 1), error = function(e) NULL)
      if (is.null(deg)) next
      
      if (gene %in% rownames(deg)) {
        row_deg <- deg[gene, , drop = FALSE]
      } else if ("gene" %in% colnames(deg)) {
        row_deg <- deg[deg$gene == gene, , drop = FALSE]
      } else next
      
      if (nrow(row_deg) == 0) next
      
      lfc_col <- dplyr::case_when(
        "avg_log2FC" %in% colnames(row_deg) ~ "avg_log2FC",
        "avg_logFC"  %in% colnames(row_deg) ~ "avg_logFC",
        TRUE ~ NA_character_
      )
      pval_col <- dplyr::case_when(
        "p_val_fdr" %in% colnames(row_deg) ~ "p_val_fdr",
        "p_val_adj" %in% colnames(row_deg) ~ "p_val_adj",
        "p_val"     %in% colnames(row_deg) ~ "p_val",
        TRUE ~ NA_character_
      )
      if (is.na(lfc_col)) next
      
      lfc  <- row_deg[[lfc_col]]
      pval <- if (!is.na(pval_col)) row_deg[[pval_col]] else NA_real_
      
      direction <- dplyr::case_when(
        !is.na(pval) & pval < 0.05 & lfc >  0.1 ~ "Up in KO",
        !is.na(pval) & pval < 0.05 & lfc < -0.1 ~ "Down in KO",
        is.na(pval)                & lfc >  0.1  ~ "Up in KO (nom.)",
        is.na(pval)                & lfc < -0.1  ~ "Down in KO (nom.)",
        TRUE                                      ~ "Not DE"
      )
      return(list(
        direction = direction,
        lfc       = round(lfc, 3),
        pval      = if (!is.na(pval_col)) round(row_deg[[pval_col]], 4) else NA
      ))
    }
    return(NULL)
  }
  
  expr_info <- get_expr_direction(gene, cell_type, contrast, deg_dir)
  
  # ── 1. Subset object ─────────────────────────────────────────────────────
  obj_ct <- if (!is.null(cell_type)) {
    tryCatch(subset(atac_obj, subset = cell_type == !!cell_type),
             error = function(e) atac_obj)
  } else {
    atac_obj
  }
  if (ncol(obj_ct) < 20) stop("Too few cells for: ", cell_type)
  
  obj_ct$genotype <- factor(obj_ct$genotype, levels = c("Ctrl", "KO"))
  obj_ct          <- clean_annotation(obj_ct)
  
  gene_resolved <- resolve_gene_alias(gene, Annotation(obj_ct))
  if (gene_resolved != gene)
    message("Alias resolved: ", gene, " \u2192 ", gene_resolved)
  
  # ── 2. Compute gene-body-aware window ────────────────────────────────────
  gr_gene  <- gene_to_region(gene_resolved)
  gene_len <- if (!is.null(gr_gene)) width(gr_gene) else 0
  flank    <- max(extend_up, gene_len)
  
  region_str <- if (!is.null(gr_gene)) {
    paste0(
      as.character(seqnames(gr_gene)), "-",
      start(gr_gene) - flank, "-",
      end(gr_gene)   + flank
    )
  } else {
    gene_resolved
  }
  
  # ── 3. Coverage per genotype ──────────────────────────────────────────────
  obj_ctrl <- clean_annotation(subset(obj_ct, subset = genotype == "Ctrl"))
  obj_ko   <- clean_annotation(subset(obj_ct, subset = genotype == "KO"))
  obj_ctrl$genotype <- factor("Ctrl", levels = "Ctrl")
  obj_ko$genotype   <- factor("KO",   levels = "KO")
  
  df_ctrl <- get_cov_df(obj_ctrl, region_str, "Ctrl", 0, 0)
  df_ko   <- get_cov_df(obj_ko,   region_str, "KO",   0, 0)
  if (is.null(df_ctrl) && is.null(df_ko))
    stop("Could not extract coverage data for ", gene)
  
  n_ctrl <- ncol(obj_ctrl)
  n_ko   <- ncol(obj_ko)
  
  ref_x    <- if (!is.null(df_ko)) df_ko$x else df_ctrl$x
  x_range  <- range(ref_x)
  n_bins   <- length(ref_x)
  common_x <- seq(x_range[1], x_range[2], length.out = n_bins)
  
  y_ctrl <- if (!is.null(df_ctrl))
    approx(df_ctrl$x, df_ctrl$y / n_ctrl, xout = common_x, rule = 2)$y
  else rep(0, n_bins)
  
  y_ko <- if (!is.null(df_ko))
    approx(df_ko$x, df_ko$y / n_ko, xout = common_x, rule = 2)$y
  else rep(0, n_bins)
  
  gaussian_smooth <- function(y, sigma = 3) {
    hw  <- ceiling(3 * sigma)
    ker <- dnorm(-hw:hw, sd = sigma)
    ker <- ker / sum(ker)
    as.numeric(stats::filter(y, ker, circular = TRUE))
  }
  
  y_ctrl <- gaussian_smooth(as.numeric(y_ctrl), sigma = 3)
  y_ko   <- gaussian_smooth(as.numeric(y_ko),   sigma = 3)
  
  diff_df <- data.frame(x    = common_x,
                        diff = as.numeric(y_ko) - as.numeric(y_ctrl))
  
  split_segments <- function(x, y) {
    segs <- list(); cur_x <- c(); cur_y <- c()
    for (k in seq_len(length(x))) {
      if (k > 1) {
        y0 <- y[k - 1]; y1 <- y[k]
        if (!is.na(y0) && !is.na(y1) && sign(y0) != sign(y1) && (y1 - y0) != 0) {
          x_zero <- x[k-1] + (x[k] - x[k-1]) * (-y0 / (y1 - y0))
          cur_x  <- c(cur_x, x_zero); cur_y <- c(cur_y, 0)
          segs[[length(segs) + 1]] <- data.frame(x = cur_x, y = cur_y)
          cur_x <- c(x_zero); cur_y <- c(0)
        }
      }
      cur_x <- c(cur_x, x[k]); cur_y <- c(cur_y, y[k])
    }
    if (length(cur_x) > 0) segs[[length(segs) + 1]] <- data.frame(x = cur_x, y = cur_y)
    segs
  }
  
  diff_colors <- c("UP" = "#b2182b", "DOWN" = "#2166ac")
  segments    <- split_segments(diff_df$x, diff_df$diff)
  seg_pos_df  <- dplyr::bind_rows(lapply(segments, function(s) if (mean(s$y, na.rm=TRUE) >= 0) s else NULL))
  seg_neg_df  <- dplyr::bind_rows(lapply(segments, function(s) if (mean(s$y, na.rm=TRUE) <  0) s else NULL))
  
  yabs <- max(abs(diff_df$diff), na.rm = TRUE) * 1.05
  if (!is.finite(yabs) || yabs == 0) yabs <- 1
  
  # ── 4. Title with expression direction ───────────────────────────────────
  dir_label     <- if (!is.null(expr_info)) expr_info$direction else NULL
  dir_color_hex <- dplyr::case_when(
    !is.null(dir_label) & grepl("Up",   dir_label, ignore.case = TRUE) ~ "#b2182b",
    !is.null(dir_label) & grepl("Down", dir_label, ignore.case = TRUE) ~ "#2166ac",
    TRUE ~ "#555555"
  )
  
  expr_tag <- if (!is.null(expr_info)) {
    lfc_str  <- if (!is.na(expr_info$lfc))  paste0("log2FC=", expr_info$lfc)  else ""
    pval_str <- if (!is.na(expr_info$pval)) paste0(", p=", expr_info$pval)    else ""
    paste0("  |  <span style='color:", dir_color_hex, "'>**Expression: ",
           dir_label, "** (", lfc_str, pval_str, ")</span>")
  } else {
    ""
  }
  
  ct_label <- if (!is.null(cell_type)) cell_type else "All cell types"
  
  # ── 5. Coverage plot ──────────────────────────────────────────────────────
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
      title    = paste0("**", gene, "**", expr_tag),
      subtitle = paste0(
        "Cell type: ", ct_label,
        "  |  Full gene + \u00b1", round(flank / 1000, 3), "kb flanks",
        "  |  <span style='color:#b2182b'>**KO > Ctrl**</span>",
        "  <span style='color:#2166ac'>**KO < Ctrl**</span>"
      ),
      y = "Delta Coverage per cell (KO - Ctrl)",
      x = NULL
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
  
  # ── 6. Gene annotation track ─────────────────────────────────────────────
  gene_track <- make_gene_track(obj_ct, gene_resolved, flank, flank,
                                ccre_gr       = if (exists("ccre_gr")) ccre_gr else NULL,
                                row_height_in = 0.35,
                                label_size    = 3)
  p_annot <- gene_track$plot
  annot_h <- gene_track$height_in
  
  cov_h   <- 4.0
  title_h <- 0.8
  
  if (!is.null(p_annot)) {
    p_final <- p_cov / p_annot +
      plot_layout(heights = unit(c(cov_h, annot_h), "in"))
    total_h <- title_h + cov_h + annot_h + 0.5
  } else {
    p_final <- p_cov
    total_h <- title_h + cov_h + 0.5
  }
  
  # ── 7. Save or return ─────────────────────────────────────────────────────
  if (save_pdf) {
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    ct_safe   <- gsub("[^A-Za-z0-9]", "_", ct_label)
    gene_safe <- gsub("[^A-Za-z0-9]", "_", gene)
    outfile   <- file.path(out_dir, paste0(ct_safe, "__", gene_safe, ".pdf"))
    ggsave(outfile, p_final, width = 10, height = total_h, limitsize = FALSE)
    message("Saved: ", outfile)
  }
  print(p_final)
  invisible(p_final)
}

####### Query genes --------
query_gene_accessibility <- function(
    gene,
    cell_type     = NULL,
    contrast      = "KOvCtrl_male",
    da_dir        = "/Users/cnbr/Downloads/Seurat_scATAC-seq/DA_results",
    upstream      = 5000,
    downstream    = 500,
    pval_thresh   = 0.05,
    plot          = TRUE,
    save_pdf      = FALSE,
    out_dir       = "/Users/cnbr/Downloads/Seurat_scATAC-seq/coverage_plots_targets/interactive",
    include_links = TRUE
) {
  
  gr_gene <- gene_to_region(gene)
  if (is.null(gr_gene)) stop("Gene not found in TxDb: ", gene)
  
  tss_gr  <- promoters(gr_gene, upstream = upstream, downstream = downstream)
  body_gr <- gr_gene
  
  flank_gr <- GenomicRanges::GRanges(
    seqnames = seqnames(gr_gene),
    ranges   = IRanges::IRanges(
      start = max(1, start(gr_gene) - upstream),
      end   = end(gr_gene) + upstream
    )
  )
  
  linked_peaks <- character(0)
  if (include_links && exists("links_sig")) {
    linked_peaks <- links_sig$peak[links_sig$gene == gene]
    if (length(linked_peaks) > 0)
      message("Linked peaks from links_sig: ", length(linked_peaks))
  }
  
  all_files <- list.files(da_dir, pattern = paste0(contrast, "\\.csv$"), full.names = TRUE)
  all_files <- all_files[!grepl("DEG", basename(all_files))]
  if (length(all_files) == 0)
    stop("No files found for contrast '", contrast, "' in ", da_dir)
  
  if (!is.null(cell_type)) {
    normalise  <- function(x) tolower(gsub("[^A-Za-z0-9]", "", x))
    query_norm <- normalise(cell_type)
    file_norms <- normalise(basename(all_files))
    all_files  <- all_files[grepl(query_norm, file_norms, fixed = TRUE)]
    if (length(all_files) == 0)
      stop("No DA file found for cell type '", cell_type, "'")
  }
  
  sample_cols <- colnames(read.csv(all_files[1], row.names = 1, nrows = 1))
  lfc_col <- dplyr::case_when(
    "avg_log2FC" %in% sample_cols ~ "avg_log2FC",
    "avg_logFC"  %in% sample_cols ~ "avg_logFC",
    TRUE ~ NA_character_
  )
  if (is.na(lfc_col))
    stop("Cannot find log2FC column. Available: ", paste(sample_cols, collapse = ", "))
  
  message("Matched files:\n", paste(basename(all_files), collapse = "\n"))
  message("Using nominal p_val for significance (single-replicate design)")
  
  results <- lapply(all_files, function(f) {
    da <- read.csv(f, row.names = 1)
    if (nrow(da) == 0) return(NULL)
    
    peak_gr <- suppressWarnings(tryCatch(
      Signac::StringToGRanges(rownames(da), sep = c("-", "-")),
      error = function(e) Signac::StringToGRanges(rownames(da))
    ))
    peak_gr$peak <- rownames(da)
    
    promo_hits <- suppressWarnings(subsetByOverlaps(peak_gr, tss_gr))
    body_hits  <- suppressWarnings(subsetByOverlaps(peak_gr, body_gr))
    flank_hits <- suppressWarnings(subsetByOverlaps(peak_gr, flank_gr))
    
    all_hits <- unique(c(
      promo_hits$peak,
      body_hits$peak,
      flank_hits$peak,
      intersect(linked_peaks, rownames(da))
    ))
    if (length(all_hits) == 0) return(NULL)
    
    ct_name <- gsub(paste0("_", contrast, "\\.csv$"), "", basename(f))
    ct_name <- gsub("_", " ", ct_name)
    
    sub_da           <- da[all_hits, , drop = FALSE]
    sub_da$peak      <- all_hits
    sub_da$cell_type <- ct_name
    sub_da$region    <- dplyr::case_when(
      all_hits %in% promo_hits$peak & all_hits %in% body_hits$peak ~ "Promoter+body",
      all_hits %in% promo_hits$peak                                 ~ "Promoter",
      all_hits %in% body_hits$peak                                  ~ "Gene body",
      all_hits %in% flank_hits$peak                                 ~ "Genomic flank",
      TRUE                                                           ~ "Linked (distal)"
    )
    sub_da$direction <- dplyr::case_when(
      sub_da[[lfc_col]] >  0.1 & sub_da$p_val < pval_thresh ~ "More open in KO",
      sub_da[[lfc_col]] < -0.1 & sub_da$p_val < pval_thresh ~ "More open in Ctrl",
      TRUE ~ "NS"
    )
    
    out_cols <- intersect(
      c(lfc_col, "pct.1", "pct.2", "p_val", "p_val_fdr", "p_val_adj"),
      colnames(sub_da)
    )
    sub_da |>
      dplyr::select(cell_type, region, peak, dplyr::all_of(out_cols), direction) |>
      dplyr::arrange(p_val)
  })
  
  result <- as.data.frame(dplyr::bind_rows(results))
  if (nrow(result) == 0) {
    message("No peaks overlapping ", gene, " in any matched file.")
    return(invisible(NULL))
  }
  
  n_sig <- sum(result$direction != "NS")
  message("\n── ", gene, " | contrast: ", contrast)
  message("Peaks found: ", nrow(result),
          "  |  Significant (nominal p<", pval_thresh, "): ", n_sig)
  
  show <- if (n_sig > 0) result[result$direction != "NS", ] else head(result, 10)
  if (n_sig == 0)
    message("No significant peaks at nominal p<", pval_thresh,
            " — showing top 10 by p_val:")
  print(show, row.names = FALSE)
  
  if (plot && exists("plot_gene_differential", mode = "function")) {
    message("\nGenerating coverage plot...")
    tryCatch(
      plot_gene_differential(
        gene,
        cell_type = cell_type,
        save_pdf  = save_pdf,
        out_dir   = out_dir
      ),
      error = function(e) message("Coverage plot unavailable: ", e$message)
    )
  }
  
  invisible(result)
}




####### Use the code ----------
# ── Usage ─────────────────────────────────────────────────────────────────────
query_gene_accessibility("Dnmt1", cell_type = "PV interneurons", save_pdf = TRUE)

####### Batch plot DA_DEG overlap genes --------

batch_plot_da_deg_overlap <- function(
    csv_path  = "/Users/cnbr/Downloads/Seurat_scATAC-seq/DA_DEG_direction_per_peak_gene.csv",
    base_dir  = "/Users/cnbr/Downloads/Seurat_scATAC-seq/coverage_plots_targets",
    contrast  = "KOvCtrl_male",
    deg_dir   = "/Users/cnbr/Downloads/Seurat_scATAC-seq/DEG_results",
    save_pdf  = TRUE
) {
  
  # ── 1. Output folders ─────────────────────────────────────────────────────
  out_concordant  <- file.path(base_dir, "DA_DEG_overlap", "Concordant")
  out_discordant  <- file.path(base_dir, "DA_DEG_overlap", "Discordant")
  dir.create(out_concordant, showWarnings = FALSE, recursive = TRUE)
  dir.create(out_discordant, showWarnings = FALSE, recursive = TRUE)
  message("Concordant  → ", out_concordant)
  message("Discordant  → ", out_discordant)
  
  # ── 2. Load CSV and deduplicate to one row per gene × cell_type ───────────
  overlap <- read.csv(csv_path, stringsAsFactors = FALSE)
  
  # One job per gene × cell_type; take the concordance from the first matching row
  # (all rows for a given gene × cell_type should share the same concordance value)
  jobs <- overlap |>
    dplyr::distinct(gene, cell_type, .keep_all = TRUE) |>
    dplyr::select(gene, cell_type, concordance)
  
  message("Total unique gene × cell-type combinations: ", nrow(jobs))
  
  # ── 3. Helper: build expected filename (mirrors plot_gene_differential) ────
  expected_filename <- function(gene, cell_type) {
    ct_safe   <- gsub("[^A-Za-z0-9]", "_", cell_type)
    gene_safe <- gsub("[^A-Za-z0-9]", "_", gene)
    paste0(ct_safe, "__", gene_safe, ".pdf")
  }
  
  # ── 4. Loop ───────────────────────────────────────────────────────────────
  n_total   <- nrow(jobs)
  n_skipped <- 0L
  n_done    <- 0L
  n_failed  <- 0L
  
  for (i in seq_len(n_total)) {
    gene        <- jobs$gene[i]
    cell_type   <- jobs$cell_type[i]
    concordance <- jobs$concordance[i]
    
    out_dir <- if (grepl("Concordant", concordance, ignore.case = TRUE)) {
      out_concordant
    } else {
      out_discordant
    }
    
    outfile <- file.path(out_dir, expected_filename(gene, cell_type))
    
    # ── Skip if already plotted ──────────────────────────────────────────
    if (file.exists(outfile)) {
      message(sprintf("[%d/%d] SKIP (exists): %s | %s [%s]",
                      i, n_total, gene, cell_type, concordance))
      n_skipped <- n_skipped + 1L
      next
    }
    
    message(sprintf("[%d/%d] Plotting: %s | %s [%s]",
                    i, n_total, gene, cell_type, concordance))
    
    tryCatch({
      plot_gene_differential(
        gene      = gene,
        cell_type = cell_type,
        contrast  = contrast,
        deg_dir   = deg_dir,
        save_pdf  = save_pdf,
        out_dir   = out_dir
      )
      n_done <- n_done + 1L
    }, error = function(e) {
      message(sprintf("  ERROR: %s", conditionMessage(e)))
      n_failed <<- n_failed + 1L
    })
  }
  
  # ── 5. Summary ────────────────────────────────────────────────────────────
  message("\n── Batch complete ──────────────────────────────────")
  message(sprintf("  Total jobs : %d", n_total))
  message(sprintf("  Plotted    : %d", n_done))
  message(sprintf("  Skipped    : %d (already existed)", n_skipped))
  message(sprintf("  Failed     : %d", n_failed))
  message(sprintf("  Concordant : %s", out_concordant))
  message(sprintf("  Discordant : %s", out_discordant))
  
  invisible(list(total = n_total, done = n_done,
                 skipped = n_skipped, failed = n_failed))
}

# ── Run the batch ─────────────────────────────────────────────────────────────
batch_plot_da_deg_overlap()