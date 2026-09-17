# 04_human_multiome_pilot.R
# Use the GSE220251 N1 sample as a pilot for human kidney RNA/ATAC analysis.
# This script tests overlap with the N1 ATAC peak catalogue; it does not test
# differential accessibility between cell types or donors.

suppressPackageStartupMessages({
  library(Seurat)
  library(GenomicRanges)
  library(IRanges)
  library(ggplot2)
})

dir.create("data/raw/GSE220251_multiome", recursive = TRUE, showWarnings = FALSE)
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
results_dir <- "Results/human_multiome_pilot"
figures_dir <- "Figures"
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

h5_file <- "data/raw/GSE220251_multiome/GSE220251_N1_filtered_feature_bc_matrix.h5"
h5_url <- paste0(
  "https://ftp.ncbi.nlm.nih.gov/geo/series/",
  "GSE220nnn/GSE220251/suppl/GSE220251_N1_filtered_feature_bc_matrix.h5"
)

if (!file.exists(h5_file)) {
  options(timeout = 1800)
  download.file(h5_url, h5_file, mode = "wb", method = "libcurl")
}

n1 <- Read10X_h5(h5_file, use.names = TRUE, unique.features = TRUE)
rna_counts <- n1[["Gene Expression"]]
atac_counts <- n1[["Peaks"]]
stopifnot(identical(colnames(rna_counts), colnames(atac_counts)))

# RNA QC and clustering ---------------------------------------------------
object <- CreateSeuratObject(
  rna_counts,
  project = "GSE220251_N1",
  min.cells = 0,
  min.features = 0
)
cells_before_qc <- ncol(object)

mitochondrial_genes <- grep("^MT-", rownames(object), value = TRUE)
object[["percent.mt"]] <- PercentageFeatureSet(object, features = mitochondrial_genes)
object <- subset(
  object,
  subset = nFeature_RNA >= 500 & nFeature_RNA <= 9000 &
    nCount_RNA >= 500 & nCount_RNA <= 60000 & percent.mt < 2
)
cells_after_qc <- ncol(object)

object <- NormalizeData(object, verbose = FALSE)
object <- FindVariableFeatures(object, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
object <- ScaleData(object, features = VariableFeatures(object), verbose = FALSE)
object <- RunPCA(object, features = VariableFeatures(object), npcs = 50, verbose = FALSE)
object <- FindNeighbors(object, dims = 1:20, verbose = FALSE)
object <- FindClusters(object, resolution = 0.5, algorithm = 1, random.seed = 123, verbose = FALSE)
object <- RunUMAP(object, dims = 1:20, seed.use = 123, verbose = FALSE)

# These labels were assigned after reviewing canonical kidney markers.
cluster_labels <- c(
  `0` = "Proximal tubule",
  `1` = "Proximal tubule",
  `2` = "B/plasma-lineage cells",
  `3` = "Thick ascending limb",
  `4` = "Collecting duct principal cell",
  `5` = "Distal convoluted tubule",
  `6` = "Thick ascending limb subtype",
  `7` = "Proximal tubule subtype",
  `8` = "Podocyte",
  `9` = "Type A intercalated cell",
  `10` = "Endothelial cell",
  `11` = "Macrophage/myeloid",
  `12` = "Fibroblast/interstitial",
  `13` = "Mixed PT/TAL-like tubular population",
  `14` = "T cell",
  `15` = "Type B/non-A intercalated cell"
)

clusters_found <- sort(unique(as.character(Idents(object))))
unlabelled_clusters <- setdiff(clusters_found, names(cluster_labels))
if (length(unlabelled_clusters)) {
  stop("No cell-type label is defined for cluster(s): ", paste(unlabelled_clusters, collapse = ", "))
}

object$CellType <- unname(cluster_labels[as.character(Idents(object))])
cluster_table <- data.frame(
  Cluster = clusters_found,
  CellType = unname(cluster_labels[clusters_found]),
  Number_of_cells = as.integer(table(factor(as.character(Idents(object)), levels = clusters_found))),
  stringsAsFactors = FALSE
)
write.csv(
  cluster_table,
  file.path(results_dir, "human_multiome_cluster_labels.csv"),
  row.names = FALSE
)

qc_summary <- data.frame(
  Sample = "GSE220251_N1",
  Cells_before_QC = cells_before_qc,
  Cells_after_QC = cells_after_qc,
  Cells_removed = cells_before_qc - cells_after_qc
)
write.csv(
  qc_summary,
  file.path(results_dir, "human_multiome_QC_summary.csv"),
  row.names = FALSE
)

p_umap <- DimPlot(object, reduction = "umap", group.by = "CellType", label = TRUE, repel = TRUE) +
  NoLegend()
ggsave(
  file.path(figures_dir, "FigureS1_Human_kidney_UMAP.png"),
  p_umap,
  width = 9,
  height = 6,
  dpi = 300
)

# Overlap with the N1 ATAC peak catalogue --------------------------------
human_loci <- read.csv(
  "Results/priority_loci_liftover/human_priority_loci_hg38.csv",
  stringsAsFactors = FALSE
)
if (nrow(human_loci) != 18) {
  warning("Expected 18 uniquely mapped loci; found ", nrow(human_loci), ".")
}

human_ranges <- GRanges(
  human_loci$Human_chr,
  IRanges(human_loci$Human_start, human_loci$Human_end)
)
names(human_ranges) <- human_loci$Probe_ID

peak_names <- rownames(atac_counts)
peak_parts <- regexec("^([^:]+):([0-9]+)-([0-9]+)$", peak_names)
peak_parts <- regmatches(peak_names, peak_parts)
valid_peaks <- lengths(peak_parts) == 4
if (!all(valid_peaks)) {
  stop(sum(!valid_peaks), " ATAC peak names do not follow chr:start-end format.")
}
peak_parts <- do.call(rbind, peak_parts)

atac_ranges <- GRanges(
  peak_parts[, 2],
  IRanges(as.integer(peak_parts[, 3]), as.integer(peak_parts[, 4]))
)
names(atac_ranges) <- peak_names

overlaps <- findOverlaps(human_ranges, atac_ranges, ignore.strand = TRUE)
overlap_counts <- tabulate(queryHits(overlaps), nbins = length(human_ranges))

nearest <- distanceToNearest(human_ranges, atac_ranges, ignore.strand = TRUE)
nearest_distance <- rep(NA_integer_, length(human_ranges))
nearest_peak <- rep(NA_character_, length(human_ranges))
nearest_distance[queryHits(nearest)] <- mcols(nearest)$distance
nearest_peak[queryHits(nearest)] <- names(atac_ranges)[subjectHits(nearest)]

overlapping_peaks <- rep(NA_character_, length(human_ranges))
for (i in seq_along(human_ranges)) {
  indices <- subjectHits(overlaps)[queryHits(overlaps) == i]
  if (length(indices)) {
    overlapping_peaks[i] <- paste(unique(names(atac_ranges)[indices]), collapse = "; ")
  }
}

atac_result <- human_loci
atac_result$ATAC_exact_overlap <- overlap_counts > 0
atac_result$Number_of_overlapping_peaks <- overlap_counts
atac_result$Overlapping_ATAC_peaks <- overlapping_peaks
atac_result$Nearest_ATAC_peak <- nearest_peak
atac_result$Distance_to_nearest_ATAC_peak_bp <- nearest_distance
write.csv(
  atac_result,
  file.path(results_dir, "human_multiome_ATAC_overlap.csv"),
  row.names = FALSE
)

saveRDS(object, "data/processed/GSE220251_N1_RNA_annotated.rds")
