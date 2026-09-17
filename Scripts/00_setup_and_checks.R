# 00_setup_and_checks.R
# Create the project folders and check the packages used by the analysis.
# Run this script from the repository root.

folders <- c(
  "data/raw",
  "data/processed",
  "data/reference/chain_files",
  "results",
  "figures"
)
invisible(lapply(folders, dir.create, recursive = TRUE, showWarnings = FALSE))

cran_packages <- c(
  "data.table", "dplyr", "ggplot2", "hdf5r", "readxl", "Seurat", "tidyr"
)
bioconductor_packages <- c(
  "GEOquery", "GenomicRanges", "IRanges", "knowYourCG", "limma",
  "rtracklayer", "sesameData"
)

required_packages <- c(cran_packages, bioconductor_packages)
available <- vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)

if (all(available)) {
  message("All required packages are available.")
} else {
  message("Missing packages: ", paste(required_packages[!available], collapse = ", "))
  message("See the project README for installation instructions.")
}
