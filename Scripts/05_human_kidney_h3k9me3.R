# 05_human_kidney_h3k9me3.R
# Exploratory comparison with H3K9me3 peaks from one adult kidney donor.
# The donor-specific result is used as supporting context, not population-level evidence.

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(IRanges)
  library(rtracklayer)
  library(ggplot2)
})

dir.create("data/raw/human_H3K9me3", recursive = TRUE, showWarnings = FALSE)
dir.create("data/reference/chain_files", recursive = TRUE, showWarnings = FALSE)
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
results_dir <- "Results/human_kidney_h3k9me3"
figures_dir <- "Figures"
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

distance_thresholds <- c(0, 1000, 5000, 10000)

# H3K9me3 peaks: hg19 to hg38 --------------------------------------------
peak_file <- "data/raw/human_H3K9me3/BI.Adult_Kidney.H3K9me3.27.broadPeak.gz"
peak_url <- paste0(
  "https://egg2.wustl.edu/roadmap/data/byFileType/peaks/unconsolidated/",
  "broadPeak/ucsc_compatible/BI.Adult_Kidney.H3K9me3.27.broadPeak.gz"
)
if (!file.exists(peak_file)) {
  download.file(peak_url, peak_file, mode = "wb", method = "libcurl")
}

peaks <- fread(peak_file, header = FALSE, sep = "\t")
if (ncol(peaks) < 9) stop("The H3K9me3 broadPeak file has fewer than nine columns.")
setnames(
  peaks,
  1:9,
  c("chrom", "chromStart", "chromEnd", "name", "score", "strand", "signalValue", "pValue", "qValue")
)

# BED coordinates are 0-based; GRanges coordinates are 1-based.
hg19_ranges <- GRanges(peaks$chrom, IRanges(peaks$chromStart + 1L, peaks$chromEnd))

chain_gz <- "data/reference/chain_files/hg19ToHg38.over.chain.gz"
chain_file <- "data/reference/chain_files/hg19ToHg38.over.chain"
if (!file.exists(chain_gz)) {
  download.file(
    "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz",
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
lifted <- liftOver(hg19_ranges, chain)
mapped_peak_count <- sum(elementNROWS(lifted) > 0)
hg38_ranges <- unlist(lifted[elementNROWS(lifted) > 0], use.names = FALSE)
standard_chromosomes <- c(paste0("chr", 1:22), "chrX", "chrY")
hg38_ranges <- hg38_ranges[as.character(seqnames(hg38_ranges)) %in% standard_chromosomes]
hg38_ranges <- reduce(hg38_ranges)
saveRDS(hg38_ranges, "data/processed/human_kidney_H3K9me3_donor27_hg38.rds")

peak_summary <- data.frame(
  Peaks_in_hg19_file = length(hg19_ranges),
  Peaks_with_at_least_one_hg38_mapping = mapped_peak_count,
  Reduced_hg38_peak_regions = length(hg38_ranges)
)
write.csv(
  peak_summary,
  file.path(results_dir, "human_H3K9me3_liftover_summary.csv"),
  row.names = FALSE
)

# Priority loci -----------------------------------------------------------
human_loci <- read.csv(
  "Results/priority_loci_liftover/human_priority_loci_hg38.csv",
  stringsAsFactors = FALSE
)
human_ranges <- GRanges(
  human_loci$Human_chr,
  IRanges(human_loci$Human_start, human_loci$Human_end)
)
names(human_ranges) <- human_loci$Probe_ID

exact_hits <- findOverlaps(human_ranges, hg38_ranges, ignore.strand = TRUE)
exact_overlap <- tabulate(queryHits(exact_hits), nbins = length(human_ranges)) > 0
nearest <- distanceToNearest(human_ranges, hg38_ranges, ignore.strand = TRUE)
nearest_distance <- rep(NA_integer_, length(human_ranges))
nearest_distance[queryHits(nearest)] <- mcols(nearest)$distance

human_loci$H3K9me3_exact_overlap <- exact_overlap
human_loci$Distance_to_nearest_H3K9me3_bp <- nearest_distance
write.csv(
  human_loci,
  file.path(results_dir, "human_priority_loci_H3K9me3.csv"),
  row.names = FALSE
)

# Mapped-probe background -------------------------------------------------
background <- readRDS("data/processed/tested_probe_background_hg38.rds")
background_ranges <- GRanges(
  background$Human_chr,
  IRanges(background$Human_start, background$Human_end)
)
names(background_ranges) <- background$Probe_ID
background_ranges <- background_ranges[!names(background_ranges) %in% human_loci$Probe_ID]

background_hits <- findOverlaps(background_ranges, hg38_ranges, ignore.strand = TRUE)
background_exact <- tabulate(queryHits(background_hits), nbins = length(background_ranges)) > 0
background_nearest <- distanceToNearest(background_ranges, hg38_ranges, ignore.strand = TRUE)
background_distance <- rep(NA_integer_, length(background_ranges))
background_distance[queryHits(background_nearest)] <- mcols(background_nearest)$distance

fisher_results <- do.call(rbind, lapply(distance_thresholds, function(threshold) {
  target_flag <- if (threshold == 0) {
    exact_overlap
  } else {
    !is.na(nearest_distance) & nearest_distance <= threshold
  }
  background_flag <- if (threshold == 0) {
    background_exact
  } else {
    !is.na(background_distance) & background_distance <= threshold
  }

  contingency_table <- matrix(c(
    sum(target_flag), sum(!target_flag),
    sum(background_flag), sum(!background_flag)
  ), nrow = 2, byrow = TRUE)
  test <- fisher.test(contingency_table, alternative = "greater")

  data.frame(
    Threshold_bp = threshold,
    Target_In = sum(target_flag),
    Target_Total = length(target_flag),
    Target_Percent = 100 * mean(target_flag),
    Background_In = sum(background_flag),
    Background_Total = length(background_flag),
    Background_Percent = 100 * mean(background_flag),
    Odds_Ratio = unname(test$estimate),
    P_Value = test$p.value
  )
}))

fisher_results$FDR <- p.adjust(fisher_results$P_Value, method = "BH")
write.csv(
  fisher_results,
  file.path(results_dir, "human_H3K9me3_target_vs_background.csv"),
  row.names = FALSE
)

plot_data <- rbind(
  data.frame(
    Threshold_bp = fisher_results$Threshold_bp,
    Group = "Priority loci",
    Percent = fisher_results$Target_Percent
  ),
  data.frame(
    Threshold_bp = fisher_results$Threshold_bp,
    Group = "Mapped probe background",
    Percent = fisher_results$Background_Percent
  )
)

h3k9me3_plot <- ggplot(
  plot_data,
  aes(Threshold_bp, Percent, colour = Group, group = Group)
) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.8) +
  scale_x_continuous(
    breaks = distance_thresholds,
    labels = scales::label_comma()
  ) +
  labs(
    title = "Human kidney H3K9me3 proximity",
    x = "Distance threshold (bp)",
    y = "Loci within threshold (%)",
    colour = NULL
  ) +
  theme_classic()

ggsave(
  file.path(figures_dir, "FigureS2_Human_kidney_H3K9me3_enrichment.png"),
  h3k9me3_plot,
  width = 7,
  height = 5,
  dpi = 300
)
