# ══════════════════════════════════════════════════════════════════════════════
# SCENIC+ EXPORT FROM R
# ══════════════════════════════════════════════════════════════════════════════

library(Matrix)
library(Signac)
library(Seurat)

setwd("~/Downloads/Seurat_scATAC-seq")
dir.create("scenicplus_input", showWarnings = FALSE)

# ── 1. RNA raw counts ─────────────────────────────────────────────────────────
DefaultAssay(rna) <- "RNA"
rna_counts <- GetAssayData(rna, assay = "RNA", layer = "counts")

writeMM(rna_counts, "scenicplus_input/rna_counts.mtx")
write.table(rownames(rna_counts), "scenicplus_input/rna_gene_names.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(colnames(rna_counts), "scenicplus_input/rna_barcodes.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE)

# RNA metadata
rna_meta <- rna@meta.data[, c("cell_type", "genotype", "sex",
                              "nCount_RNA", "nFeature_RNA"), drop = FALSE]
rna_meta$cell_type <- as.character(rna_meta$cell_type)
rna_meta$genotype  <- as.character(rna_meta$genotype)
write.csv(rna_meta, "scenicplus_input/rna_metadata.csv")

# RNA UMAP
rna_umap <- as.data.frame(Embeddings(rna, reduction = "umap"))
colnames(rna_umap) <- c("UMAP_1", "UMAP_2")
write.csv(rna_umap, "scenicplus_input/rna_umap.csv")

message("RNA exported: ", nrow(rna_counts), " genes x ", ncol(rna_counts), " cells")

# ── 2. ATAC cell metadata ─────────────────────────────────────────────────────
atac_meta <- obj_male@meta.data[, intersect(
  c("cell_type", "genotype", "sex", "sample",
    "nCount_peaks", "nFeature_peaks",
    "blacklist_fraction", "nucleosome_signal", "TSS.enrichment"),
  colnames(obj_male@meta.data)
), drop = FALSE]
atac_meta$cell_type <- as.character(atac_meta$cell_type)
atac_meta$genotype  <- as.character(atac_meta$genotype)
atac_meta$barcode   <- rownames(atac_meta)

# SCENIC+ needs sample_id — use one sample since this is aggregated
atac_meta$sample_id <- "male_cortex"

write.csv(atac_meta, "scenicplus_input/atac_metadata.csv")

# ATAC UMAP
atac_umap <- as.data.frame(Embeddings(obj_male, reduction = "umap"))
colnames(atac_umap) <- c("UMAP_1", "UMAP_2")
write.csv(atac_umap, "scenicplus_input/atac_umap.csv")

message("ATAC metadata exported: ", nrow(atac_meta), " cells")

# ── 3. Verify fragment file ───────────────────────────────────────────────────
frag_path <- "~/Downloads/Seurat_scATAC-seq/fragments.tsv.gz"
tbi_path  <- paste0(frag_path, ".tbi")
message("\nFragment file exists: ", file.exists(frag_path))
message("Fragment index exists: ", file.exists(tbi_path))
message("\nAll files in scenicplus_input/:")
print(list.files("scenicplus_input"))
message("\n=== R export complete ===")