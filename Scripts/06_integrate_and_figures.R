# 06_integrate_and_figures.R
# Combine the human ATAC-peak and kidney H3K9me3 results.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

results_dir <- "Results/integrate_and_figures"
figures_dir <- "Figures"
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

h3k9me3_distance_bp <- 1000

atac <- read.csv(
  "Results/human_multiome_pilot/human_multiome_ATAC_overlap.csv",
  stringsAsFactors = FALSE
)
h3k9me3 <- read.csv(
  "Results/human_kidney_h3k9me3/human_priority_loci_H3K9me3.csv",
  stringsAsFactors = FALSE
)

required_atac_columns <- c(
  "Probe_ID", "ATAC_exact_overlap", "Overlapping_ATAC_peaks",
  "Nearest_ATAC_peak", "Distance_to_nearest_ATAC_peak_bp"
)
required_h3_columns <- c("Probe_ID", "H3K9me3_exact_overlap", "Distance_to_nearest_H3K9me3_bp")
if (!all(required_atac_columns %in% names(atac))) stop("Required ATAC columns are missing.")
if (!all(required_h3_columns %in% names(h3k9me3))) stop("Required H3K9me3 columns are missing.")

integrated <- merge(
  h3k9me3,
  atac[, required_atac_columns],
  by = "Probe_ID",
  all.x = TRUE,
  sort = FALSE
)
integrated <- integrated[match(h3k9me3$Probe_ID, integrated$Probe_ID), ]

integrated$H3K9me3_within_1kb <-
  !is.na(integrated$Distance_to_nearest_H3K9me3_bp) &
  integrated$Distance_to_nearest_H3K9me3_bp <= h3k9me3_distance_bp
integrated$ATAC_exact_overlap[is.na(integrated$ATAC_exact_overlap)] <- FALSE

integrated$Chromatin_context <- with(integrated, ifelse(
  ATAC_exact_overlap & H3K9me3_within_1kb,
  "ATAC_overlap_and_near_H3K9me3",
  ifelse(
    ATAC_exact_overlap,
    "ATAC_overlap_only",
    ifelse(H3K9me3_within_1kb, "Near_H3K9me3_only", "Neither")
  )
))
write.csv(
  integrated,
  file.path(results_dir, "human_chromatin_context_integrated.csv"),
  row.names = FALSE
)

context_counts <- as.data.frame(table(integrated$Chromatin_context), stringsAsFactors = FALSE)
names(context_counts) <- c("Chromatin_context", "Number_of_loci")
write.csv(
  context_counts,
  file.path(results_dir, "human_chromatin_context_summary.csv"),
  row.names = FALSE
)

figure_data <- integrated %>%
  select(Probe_ID, ATAC_exact_overlap, H3K9me3_within_1kb) %>%
  rename(
    `ATAC peak overlap` = ATAC_exact_overlap,
    `Within 1 kb of H3K9me3` = H3K9me3_within_1kb
  ) %>%
  pivot_longer(-Probe_ID, names_to = "Chromatin_feature", values_to = "Present")
figure_data$Probe_ID <- factor(figure_data$Probe_ID, levels = rev(integrated$Probe_ID))

p_context <- ggplot(figure_data, aes(Chromatin_feature, Probe_ID, fill = Present)) +
  geom_tile(colour = "white", linewidth = 0.7) +
  scale_fill_manual(
    values = c(`TRUE` = "black", `FALSE` = "grey90"),
    labels = c(`FALSE` = "No", `TRUE` = "Yes"),
    name = "Feature present"
  ) +
  labs(
    title = "Chromatin context of human-mapped priority loci",
    x = NULL,
    y = "Human-mapped priority locus"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 20, hjust = 1),
    plot.title = element_text(face = "bold")
  )

ggsave(
  file.path(figures_dir, "Figure2_Human_chromatin_context.png"),
  p_context,
  width = 7,
  height = 7,
  dpi = 300
)
ggsave(
  file.path(figures_dir, "Figure2_Human_chromatin_context.pdf"),
  p_context,
  width = 7,
  height = 7
)
