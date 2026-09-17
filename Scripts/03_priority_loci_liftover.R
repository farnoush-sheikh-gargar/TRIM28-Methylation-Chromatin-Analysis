# 03_priority_loci_liftover.R
# Map the 31 priority mouse loci, and the tested-probe background, to hg38.

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(IRanges)
  library(rtracklayer)
  library(sesameData)
  library(ggplot2)
})

dir.create("data/reference/chain_files", recursive = TRUE, showWarnings = FALSE)
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
results_dir <- "Results/priority_loci_liftover"
figures_dir <- "Figures"
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

priority <- read.csv(
  "Results/mouse_chromatin_context/mouse_priority_loci_mm10.csv",
  stringsAsFactors = FALSE
)
manifest <- sesameData_getManifestGRanges(platform = "MM285", genome = "mm10")

make_probe_ranges <- function(probe_ids) {
  probe_ids <- intersect(probe_ids, names(manifest))
  ranges <- manifest[match(probe_ids, names(manifest))]
  names(ranges) <- probe_ids
  ranges
}

chain_gz <- "data/reference/chain_files/mm10ToHg38.over.chain.gz"
chain_file <- "data/reference/chain_files/mm10ToHg38.over.chain"

if (!file.exists(chain_gz)) {
  download.file(
    "https://hgdownload.soe.ucsc.edu/goldenPath/mm10/liftOver/mm10ToHg38.over.chain.gz",
    chain_gz,
    mode = "wb",
    method = "libcurl"
  )
}

if (!file.exists(chain_file)) {
  input <- gzfile(chain_gz, "rb")
  output <- file(chain_file, "wb")
  repeat {
    bytes <- readBin(input, "raw", n = 1024 * 1024)
    if (!length(bytes)) break
    writeBin(bytes, output)
  }
  close(input)
  close(output)
}

chain <- import.chain(chain_file)

lift_unique <- function(ranges) {
  lifted <- liftOver(ranges, chain)
  mapping_count <- elementNROWS(lifted)
  keep <- mapping_count == 1
  unique_ranges <- unlist(lifted[keep], use.names = FALSE)

  mapped <- data.frame(
    Probe_ID = names(ranges)[keep],
    Human_chr = as.character(seqnames(unique_ranges)),
    Human_start = start(unique_ranges),
    Human_end = end(unique_ranges),
    stringsAsFactors = FALSE
  )
  status <- data.frame(
    Probe_ID = names(ranges),
    Number_of_hg38_mappings = mapping_count,
    Mapping_status = ifelse(mapping_count == 0, "Unmapped", ifelse(mapping_count == 1, "Unique", "Multiple")),
    stringsAsFactors = FALSE
  )
  list(mapped = mapped, status = status)
}

# Priority loci -----------------------------------------------------------
priority_lift <- lift_unique(make_probe_ranges(priority$Probe_ID))
priority_hg38 <- merge(priority, priority_lift$mapped, by = "Probe_ID", all.y = TRUE, sort = FALSE)
priority_hg38 <- priority_hg38[match(priority_lift$mapped$Probe_ID, priority_hg38$Probe_ID), ]

write.csv(
  priority_hg38,
  file.path(results_dir, "human_priority_loci_hg38.csv"),
  row.names = FALSE
)
write.csv(
  priority_lift$status,
  file.path(results_dir, "priority_loci_liftover_status.csv"),
  row.names = FALSE
)

liftover_summary <- as.data.frame(table(priority_lift$status$Mapping_status), stringsAsFactors = FALSE)
names(liftover_summary) <- c("Mapping_status", "Number_of_loci")
write.csv(
  liftover_summary,
  file.path(results_dir, "priority_loci_liftover_summary.csv"),
  row.names = FALSE
)

liftover_plot <- ggplot(
  liftover_summary,
  aes(Mapping_status, Number_of_loci, fill = Mapping_status)
) +
  geom_col(width = 0.65, show.legend = FALSE) +
  geom_text(aes(label = Number_of_loci), vjust = -0.4) +
  labs(
    title = "Mouse priority loci mapped to hg38",
    x = NULL,
    y = "Number of loci"
  ) +
  theme_classic()

ggsave(
  file.path(figures_dir, "Figure3_Priority_loci_liftover_status.png"),
  liftover_plot,
  width = 6,
  height = 5,
  dpi = 300
)

if (nrow(priority_hg38) != 18) {
  warning("Expected 18 uniquely mapped priority loci; found ", nrow(priority_hg38), ".")
}

# Tested-probe background -------------------------------------------------
beta_complete <- readRDS(
  "Data/processed/GSE229030_beta_complete.rds"
)
background_lift <- lift_unique(make_probe_ranges(rownames(beta_complete)))
background_hg38 <- background_lift$mapped
standard_chromosomes <- c(paste0("chr", 1:22), "chrX", "chrY")
background_hg38 <- background_hg38[background_hg38$Human_chr %in% standard_chromosomes, ]

saveRDS(background_hg38, "data/processed/tested_probe_background_hg38.rds")
message("Unique priority mappings: ", nrow(priority_hg38), " / ", nrow(priority))
