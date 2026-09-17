# 99_session_info.R
# Record the R session and the versions of packages used by the project.

results_dir <- "Results/session_info"
figures_dir <- "Figures"
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

project_packages <- c(
  "data.table", "dplyr", "GEOquery", "GenomicRanges", "ggplot2", "hdf5r",
  "IRanges", "knowYourCG", "limma", "readxl", "rtracklayer", "sesameData",
  "Seurat", "tidyr"
)

package_versions <- vapply(project_packages, function(package) {
  if (requireNamespace(package, quietly = TRUE)) {
    as.character(packageVersion(package))
  } else {
    "not installed"
  }
}, character(1))

output <- c(
  "PROJECT PACKAGE VERSIONS",
  paste(names(package_versions), package_versions, sep = ": "),
  "",
  "SESSION INFO",
  capture.output(sessionInfo())
)
writeLines(
  output,
  file.path(results_dir, "session_info.txt")
)

package_table <- data.frame(
  Package = names(package_versions),
  Version = unname(package_versions),
  Status = ifelse(package_versions == "not installed", "Not installed", "Installed"),
  stringsAsFactors = FALSE
)
package_table$Package <- factor(
  package_table$Package,
  levels = rev(package_table$Package)
)

session_plot <- ggplot2::ggplot(
  package_table,
  ggplot2::aes(Version, Package, colour = Status)
) +
  ggplot2::geom_point(size = 2.8) +
  ggplot2::labs(
    title = paste("Project package versions —", R.version.string),
    x = "Installed version",
    y = NULL,
    colour = NULL
  ) +
  ggplot2::theme_classic() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
    plot.title = ggplot2::element_text(size = 11)
  )

ggplot2::ggsave(
  file.path(figures_dir, "FigureS3_R_session_package_versions.png"),
  session_plot,
  width = 9,
  height = 6,
  dpi = 300
)
