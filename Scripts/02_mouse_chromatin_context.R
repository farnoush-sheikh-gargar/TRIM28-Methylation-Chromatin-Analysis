# 02_mouse_chromatin_context.R
# Test whether the published DMLs overlap repressive chromatin annotations.
#
# Inputs:
#   Results/mouse_methylation_qc/mouse_published_1133_DMLs.csv
#   Data/processed/GSE229030_beta_complete.rds
#
# Produces enrichment tables, a 31-locus priority set, and Figure 1.

suppressPackageStartupMessages({
  library(sesameData)
  library(knowYourCG)
  library(GenomicRanges)
  library(ggplot2)
})

results_dir <- "Results/mouse_chromatin_context"
figures_dir <- "Figures"
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

dml <- read.csv(
  "Results/mouse_methylation_qc/mouse_published_1133_DMLs.csv",
  stringsAsFactors = FALSE
)
beta_complete <- readRDS(
  "Data/processed/GSE229030_beta_complete.rds"
)
tested_probes <- rownames(beta_complete)

# MM285 coordinates -------------------------------------------------------
manifest <- sesameData_getManifestGRanges(platform = "MM285", genome = "mm10")
mapped_probes <- intersect(dml$Probe_ID, names(manifest))
dml_gr <- manifest[match(mapped_probes, names(manifest))]

coordinates <- data.frame(
  Probe_ID = names(dml_gr),
  Chromosome = as.character(seqnames(dml_gr)),
  Position = start(dml_gr),
  stringsAsFactors = FALSE
)
dml_annotated <- merge(dml, coordinates, by = "Probe_ID", all.x = TRUE, sort = FALSE)
dml_annotated <- dml_annotated[match(dml$Probe_ID, dml_annotated$Probe_ID), ]
write.csv(
  dml_annotated,
  file.path(results_dir, "mouse_DML_mm10_annotation.csv"),
  row.names = FALSE
)

# Chromatin annotations ---------------------------------------------------
histone_marks <- getDBs("MM285.HMconsensus")
tf_binding <- getDBs("MM285.TFBSconsensus")
chromhmm <- getDBs("MM285.chromHMM")

required_marks <- c("H3K9me3", "H3K27me3")
if (!all(required_marks %in% names(histone_marks))) {
  stop("H3K9me3 or H3K27me3 is missing from the MM285 histone-mark database.")
}
if (!"TRIM28" %in% names(tf_binding)) {
  stop("TRIM28 is missing from the MM285 TFBS database.")
}
if (!"Het" %in% names(chromhmm)) {
  stop("Het is missing from the MM285 chromHMM database.")
}

chromatin_sets <- c(
  histone_marks[required_marks],
  tf_binding[intersect(c("TRIM28", "EZH2", "SUZ12", "RNF2", "BMI1", "CBX2", "CBX3", "CBX7", "CBX8"), names(tf_binding))],
  chromhmm[intersect(c("Het", "ReprPC", "ReprPCWk"), names(chromhmm))]
)
chromatin_sets <- lapply(chromatin_sets, intersect, tested_probes)

all_dml <- intersect(dml$Probe_ID, tested_probes)
heavy_up <- intersect(dml$Probe_ID[dml$Paper_DML == "Obese_UP"], tested_probes)
light_up <- intersect(dml$Probe_ID[dml$Paper_DML == "Lean_UP"], tested_probes)

enrichment_test <- function(query, query_name, annotation, annotation_name, universe) {
  query <- intersect(query, universe)
  annotation <- intersect(annotation, universe)
  background <- setdiff(universe, query)

  contingency_table <- matrix(c(
    sum(query %in% annotation), sum(!query %in% annotation),
    sum(background %in% annotation), sum(!background %in% annotation)
  ), nrow = 2, byrow = TRUE)

  test <- fisher.test(contingency_table, alternative = "greater")
  data.frame(
    Query_Group = query_name,
    Annotation = annotation_name,
    Query_Total = length(query),
    Query_In_Annotation = sum(query %in% annotation),
    Query_Percent = 100 * mean(query %in% annotation),
    Background_Percent = 100 * length(annotation) / length(universe),
    Odds_Ratio = unname(test$estimate),
    P_Value = test$p.value,
    stringsAsFactors = FALSE
  )
}

query_sets <- list(All_DML = all_dml, Heavy_UP = heavy_up, Light_UP = light_up)
enrichment_results <- do.call(rbind, lapply(names(query_sets), function(query_name) {
  do.call(rbind, lapply(names(chromatin_sets), function(annotation_name) {
    enrichment_test(
      query_sets[[query_name]], query_name,
      chromatin_sets[[annotation_name]], annotation_name,
      tested_probes
    )
  }))
}))
enrichment_results$FDR <- ave(
  enrichment_results$P_Value,
  enrichment_results$Query_Group,
  FUN = function(p) p.adjust(p, method = "BH")
)
write.csv(
  enrichment_results,
  file.path(results_dir, "mouse_chromatin_enrichment.csv"),
  row.names = FALSE
)

# Direction-specific comparison -----------------------------------------
h3k9me3 <- chromatin_sets[["H3K9me3"]]
h3k27me3 <- chromatin_sets[["H3K27me3"]]
trim28 <- chromatin_sets[["TRIM28"]]

comparison_sets <- list(
  H3K9me3 = h3k9me3,
  H3K27me3 = h3k27me3,
  TRIM28_H3K9me3 = intersect(trim28, h3k9me3),
  H3K9me3_H3K27me3 = intersect(h3k9me3, h3k27me3),
  Het = chromatin_sets[["Het"]]
)

compare_directions <- function(annotation) {
  heavy_in <- sum(heavy_up %in% annotation)
  light_in <- sum(light_up %in% annotation)
  table <- matrix(c(
    heavy_in, length(heavy_up) - heavy_in,
    light_in, length(light_up) - light_in
  ), nrow = 2, byrow = TRUE)
  test <- fisher.test(table, alternative = "two.sided")

  data.frame(
    Heavy_Total = length(heavy_up),
    Heavy_In = heavy_in,
    Heavy_Percent = 100 * heavy_in / length(heavy_up),
    Light_Total = length(light_up),
    Light_In = light_in,
    Light_Percent = 100 * light_in / length(light_up),
    Odds_Ratio_Heavy_vs_Light = unname(test$estimate),
    P_Value = test$p.value
  )
}

direction_results <- do.call(rbind, lapply(names(comparison_sets), function(name) {
  cbind(Annotation = name, compare_directions(comparison_sets[[name]]))
}))
direction_results$FDR <- p.adjust(direction_results$P_Value, method = "BH")
write.csv(
  direction_results,
  file.path(results_dir, "mouse_heavy_vs_light_chromatin.csv"),
  row.names = FALSE
)

# Heavy-up DMLs in both TRIM28 and H3K9me3 define the priority set.
priority_ids <- Reduce(intersect, list(heavy_up, trim28, h3k9me3))
priority_loci <- dml_annotated[dml_annotated$Probe_ID %in% priority_ids, , drop = FALSE]
if (nrow(priority_loci) != 31) {
  stop("Expected 31 priority loci; found ", nrow(priority_loci), ".")
}

write.csv(
  priority_loci,
  file.path(results_dir, "mouse_priority_loci_mm10.csv"),
  row.names = FALSE
)
writeLines(
  priority_loci$Probe_ID,
  file.path(results_dir, "mouse_priority_probe_IDs.txt")
)

# Main mouse figure -------------------------------------------------------
plot_order <- c("H3K9me3", "H3K27me3", "TRIM28_H3K9me3", "H3K9me3_H3K27me3", "Het")
plot_data <- direction_results[match(plot_order, direction_results$Annotation), ]
plot_data <- rbind(
  data.frame(Annotation = plot_data$Annotation, Group = "Heavy-up DMLs", Percent = plot_data$Heavy_Percent),
  data.frame(Annotation = plot_data$Annotation, Group = "Light-up DMLs", Percent = plot_data$Light_Percent)
)
plot_data$Annotation <- factor(plot_data$Annotation, levels = plot_order)

p_mouse <- ggplot(plot_data, aes(Annotation, Percent, fill = Group)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  labs(
    title = "Direction-specific overlap with repressive chromatin",
    x = NULL,
    y = "DMLs overlapping annotation (%)"
  ) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

ggsave(
  file.path(figures_dir, "Figure1_Mouse_chromatin_comparison.png"),
  p_mouse,
  width = 8,
  height = 5,
  dpi = 300
)
