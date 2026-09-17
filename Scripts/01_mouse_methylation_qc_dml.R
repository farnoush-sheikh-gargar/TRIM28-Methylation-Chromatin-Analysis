# 01_mouse_methylation_qc_dml.R
# Recover the published Light/Heavy DML set and run sample-level QC.
# Dataset: GSE229030, Trim28+/D9 Light versus Heavy mice.
#
# Produces:
#   Results/mouse_methylation_qc/
#   Data/processed/
#   Figures/

suppressPackageStartupMessages({
  library(GEOquery)
  library(data.table)
  library(limma)
  library(ggplot2)
  library(readxl)
})

dir.create("data/raw", recursive = TRUE, showWarnings = FALSE)
results_dir <- "Results/mouse_methylation_qc"
processed_dir <- "Data/processed"
figures_dir <- "Figures"
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

detection_p_cutoff <- 0.05
delta_beta_cutoff <- 0.05
nominal_p_cutoff <- 0.05

# Sample metadata ---------------------------------------------------------
gse <- getGEO("GSE229030", GSEMatrix = TRUE, getGPL = FALSE)
meta <- pData(gse[[1]])

is_light <- grepl("Trim28_light", meta$title, ignore.case = TRUE)
is_heavy <- grepl("Trim28_heavy", meta$title, ignore.case = TRUE)
keep <- is_light | is_heavy

sample_info <- data.frame(
  GEO_ID = rownames(meta)[keep],
  Sample = meta$title[keep],
  Group = ifelse(is_light[keep], "Light", "Heavy"),
  stringsAsFactors = FALSE
)

if (sum(sample_info$Group == "Light") != 9 || sum(sample_info$Group == "Heavy") != 15) {
  stop("The expected 9 Light and 15 Heavy samples were not found in the GEO metadata.")
}

meta_lh <- meta[sample_info$GEO_ID, , drop = FALSE]
green_files <- meta_lh$supplementary_file
array_id <- sub(".*_([0-9]+_R[0-9]+C[0-9]+)_Grn\\.idat\\.gz$", "\\1", green_files)

batch_columns <- grep("characteristics_ch1", colnames(meta_lh), value = TRUE)
batch_column <- batch_columns[vapply(batch_columns, function(column) {
  any(grepl("^batch:", meta_lh[[column]], ignore.case = TRUE))
}, logical(1))][1]

if (is.na(batch_column)) stop("Could not find the batch column in the GEO metadata.")
batch <- sub("^batch:\\s*", "", meta_lh[[batch_column]], ignore.case = TRUE)

mapping <- data.frame(
  GEO_ID = sample_info$GEO_ID,
  Sample = sample_info$Sample,
  Group = sample_info$Group,
  Batch = batch,
  Array_ID = array_id,
  Beta_column = paste0("Sample_", array_id),
  DetectionP_column = paste0("Detection_Pvalue_", array_id),
  stringsAsFactors = FALSE
)
write.csv(
  mapping,
  file.path(processed_dir, "GSE229030_sample_mapping.csv"),
  row.names = FALSE
)

# Processed methylation matrix -------------------------------------------
processed_file <- "data/raw/GSE229030_processed_matrix.txt.gz"
processed_url <- paste0(
  "https://ftp.ncbi.nlm.nih.gov/geo/series/",
  "GSE229nnn/GSE229030/suppl/GSE229030_processed_matrix.txt.gz"
)

if (!file.exists(processed_file)) {
  options(timeout = 1800)
  download.file(processed_url, processed_file, mode = "wb", method = "libcurl")
}

matrix_columns <- colnames(read.delim(processed_file, nrows = 0, check.names = FALSE))
required_columns <- c(mapping$Beta_column, mapping$DetectionP_column)
missing_columns <- setdiff(required_columns, matrix_columns)
if (length(missing_columns)) {
  stop("Columns missing from the processed matrix: ", paste(missing_columns, collapse = ", "))
}

columns_to_read <- c("ID_REF", required_columns)
methylation <- fread(
  processed_file,
  select = columns_to_read,
  na.strings = "NA",
  data.table = FALSE
)

beta <- as.matrix(methylation[, mapping$Beta_column, drop = FALSE])
detection_p <- as.matrix(methylation[, mapping$DetectionP_column, drop = FALSE])
rownames(beta) <- rownames(detection_p) <- methylation$ID_REF
colnames(beta) <- colnames(detection_p) <- mapping$GEO_ID

# QC and PCA --------------------------------------------------------------
sample_qc <- data.frame(
  GEO_ID = mapping$GEO_ID,
  Group = mapping$Group,
  Batch = mapping$Batch,
  Missing_Beta = colSums(is.na(beta)),
  Detection_Rate = colMeans(detection_p < detection_p_cutoff, na.rm = TRUE)
)
write.csv(
  sample_qc,
  file.path(results_dir, "mouse_methylation_sample_QC.csv"),
  row.names = FALSE
)

beta_complete <- beta[complete.cases(beta), , drop = FALSE]
if (nrow(beta_complete) != 250882) {
  stop("Expected 250,882 complete probes; found ", format(nrow(beta_complete), big.mark = ","), ".")
}

rm(methylation, beta, detection_p)
invisible(gc())

pca <- prcomp(
  t(beta_complete),
  center = TRUE,
  scale. = FALSE,
  rank. = 3
)
variance_explained <- 100 * pca$sdev^2 / sum(pca$sdev^2)
pca_data <- data.frame(
  GEO_ID = rownames(pca$x),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  PC3 = pca$x[, 3],
  Group = factor(mapping$Group[match(rownames(pca$x), mapping$GEO_ID)], levels = c("Light", "Heavy")),
  Batch = factor(mapping$Batch[match(rownames(pca$x), mapping$GEO_ID)])
)

p_pca <- ggplot(pca_data, aes(PC1, PC2, colour = Group, shape = Batch, label = GEO_ID)) +
  geom_point(size = 3.6) +
  geom_text(vjust = -0.7, size = 2.5, check_overlap = TRUE) +
  labs(
    title = "GSE229030: Trim28 Light versus Heavy",
    x = sprintf("PC1 (%.1f%%)", variance_explained[1]),
    y = sprintf("PC2 (%.1f%%)", variance_explained[2])
  ) +
  theme_classic()

ggsave(
  file.path(figures_dir, "Figure0_Mouse_methylation_PCA.png"),
  p_pca,
  width = 8,
  height = 6,
  dpi = 300
)
write.csv(
  pca_data,
  file.path(results_dir, "mouse_methylation_PCA_coordinates.csv"),
  row.names = FALSE
)
rm(pca)
invisible(gc())

# Batch-adjusted Light/Heavy comparison ----------------------------------
phenotype <- data.frame(
  Group = factor(mapping$Group, levels = c("Light", "Heavy")),
  Batch = factor(mapping$Batch),
  row.names = mapping$GEO_ID
)
stopifnot(identical(colnames(beta_complete), rownames(phenotype)))

design <- model.matrix(~ Batch + Group, data = phenotype)
fit <- eBayes(lmFit(beta_complete, design), robust = TRUE)
all_results <- topTable(fit, coef = "GroupHeavy", number = Inf, sort.by = "none")
all_results$Probe_ID <- rownames(all_results)
all_results$DeltaBeta_Heavy_minus_Light <- all_results$logFC
all_results$Direction <- ifelse(
  all_results$DeltaBeta_Heavy_minus_Light > 0,
  "Higher_methylation_in_Heavy",
  "Higher_methylation_in_Light"
)

# This nominal-p set is used as a sensitivity analysis. The published DMLs
# below remain the main set carried into the chromatin analysis.
is_sensitivity_dml <-
  abs(all_results$DeltaBeta_Heavy_minus_Light) >= delta_beta_cutoff &
  all_results$P.Value < nominal_p_cutoff
sensitivity_dml <- all_results[is_sensitivity_dml, ]

write.csv(
  all_results,
  file.path(results_dir, "mouse_methylation_all_probes.csv"),
  row.names = FALSE
)
write.csv(
  sensitivity_dml,
  file.path(results_dir, "mouse_methylation_sensitivity_DMLs.csv"),
  row.names = FALSE
)

# Published Fig. 5B DMLs --------------------------------------------------
source_workbook <- "data/raw/43018_2024_900_MOESM6_ESM.xlsx"
if (!file.exists(source_workbook)) {
  stop("Place 43018_2024_900_MOESM6_ESM.xlsx in data/raw and run the script again.")
}

figure5b <- read_excel(source_workbook, sheet = "Sheet1", skip = 1123, n_max = 250882)
if (!"DML" %in% names(figure5b)) stop("The DML column was not found in the source-data sheet.")

if (!"Probe_ID" %in% names(figure5b)) {
  overlap_counts <- vapply(figure5b, function(column) {
    sum(as.character(column) %in% all_results$Probe_ID, na.rm = TRUE)
  }, numeric(1))
  probe_column <- names(which.max(overlap_counts))
  if (!length(probe_column) || max(overlap_counts) == 0) {
    stop("Could not identify the probe-ID column in the source-data sheet.")
  }
  names(figure5b)[names(figure5b) == probe_column] <- "Probe_ID"
}

published_dml <- figure5b[
  !is.na(figure5b$DML) & figure5b$DML %in% c("Obese_UP", "Lean_UP"),
  ,
  drop = FALSE
]
published_dml$Paper_DML <- published_dml$DML

recovered <- merge(
  published_dml[, c("Probe_ID", "Paper_DML")],
  sensitivity_dml[, c("Probe_ID", "DeltaBeta_Heavy_minus_Light", "Direction")],
  by = "Probe_ID",
  all.x = TRUE,
  sort = FALSE
)
recovered <- recovered[match(published_dml$Probe_ID, recovered$Probe_ID), ]

if (nrow(recovered) != 1133 || anyNA(recovered$DeltaBeta_Heavy_minus_Light)) {
  stop("The 1,133 published DMLs were not fully recovered by the sensitivity analysis.")
}

write.csv(
  recovered,
  file.path(results_dir, "mouse_published_1133_DMLs.csv"),
  row.names = FALSE
)
saveRDS(
  beta_complete,
  file.path(processed_dir, "GSE229030_beta_complete.rds")
)

summary_lines <- c(
  paste("Light samples:", sum(mapping$Group == "Light")),
  paste("Heavy samples:", sum(mapping$Group == "Heavy")),
  paste("Complete probes:", nrow(beta_complete)),
  paste("Sensitivity DMLs:", nrow(sensitivity_dml)),
  paste("Published DMLs recovered:", nrow(recovered))
)
writeLines(
  summary_lines,
  file.path(results_dir, "mouse_methylation_summary.txt")
)
