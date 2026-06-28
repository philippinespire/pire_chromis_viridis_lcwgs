#!/usr/bin/env Rscript

# ==============================================================================
# SCRIPT NAME:  plot_dystruct.R
# PURPOSE:      Visualizes DYSTRUCT ancestry proportions into a PDF barplot.
#
# USAGE:        Rscript plot_dystruct.R <input.theta> <sample_ids.txt|samplesheet.csv> <output.pdf>
#
# INPUTS:       1. <input.theta>      - Whitespace-separated DYSTRUCT ancestry
#                                      matrix with one row per sample and one
#                                      column per ancestral population (K).
#               2. <sample_ids.txt|samplesheet.csv>
#                                    - Either one sample label per line, in the
#                                      same row order as the theta file, or a
#                                      CSV with a header row containing a
#                                      sample column (preferred name: sample) and
#                                      optional era column for grouping.
#               3. <output.pdf>       - Target filename for the output PDF.
#
# REQUIREMENTS: Base R installation.
# DATE:         June 2026
# AUTHOR:       Malin Pinsky with assistance from GitHub Copilot with the GPT-5.3-Codex model.
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("\nUsage: Rscript plot_dystruct.R <input.theta> <sample_ids.txt|samplesheet.csv> <output.pdf>\n", call. = FALSE)
}

theta_file  <- args[1]
sample_file <- args[2]
output_pdf  <- args[3]

if (!file.exists(theta_file)) stop(paste("Error: File not found -", theta_file))
if (!file.exists(sample_file)) stop(paste("Error: File not found -", sample_file))

theta <- read.table(theta_file, header = FALSE)

sample_data <- NULL
sample_ids <- NULL
sample_era <- NULL

is_csv_samplesheet <- grepl("\\.csv$", tolower(sample_file))

if (is_csv_samplesheet) {
  sample_data <- read.csv(sample_file, stringsAsFactors = FALSE, check.names = FALSE)
  if (!"sample" %in% colnames(sample_data)) {
    stop("Error: CSV samplesheet must contain a 'sample' column.")
  }
  sample_ids <- as.character(sample_data$sample)
  if ("era" %in% colnames(sample_data)) {
    sample_era <- as.character(sample_data$era)
  }
} else {
  sample_ids <- readLines(sample_file, warn = FALSE)
}

if (nrow(theta) != length(sample_ids)) {
  stop(paste(
    "Error: Row mismatch! Theta file has", nrow(theta), "rows, but sample file has", length(sample_ids), "entries."
  ))
}

if (ncol(theta) < 2) {
  stop("Error: Theta file must contain at least 2 population columns.")
}

K <- ncol(theta)

dir.create(dirname(output_pdf), recursive = TRUE, showWarnings = FALSE)

# Remove blank labels if present, but keep row alignment strict.
if (any(sample_ids == "")) {
  stop("Error: Sample id file contains blank lines.")
}

if (any(is.na(theta))) {
  stop("Error: Theta file contains missing values; expected a complete numeric matrix.")
}

plot_data <- data.frame(
  sample = sample_ids,
  stringsAsFactors = FALSE
)

if (!is.null(sample_era)) {
  if (length(sample_era) != length(sample_ids)) {
    stop("Error: samplesheet 'era' column length does not match theta rows.")
  }
  plot_data$era <- sample_era
  plot_data <- plot_data[order(plot_data$era), , drop = FALSE]
  theta <- theta[match(plot_data$sample, sample_ids), , drop = FALSE]
}

# Use the input row order and sample IDs as labels unless we sorted by era.
sample_labels <- plot_data$sample

# Dynamic palette selection.
if (K <= 12) {
  my_colors <- palette.colors(n = K, palette = "Set 3")
} else {
  my_colors <- grDevices::rainbow(K)
}

pdf(file = output_pdf, width = 10, height = 6)
par(mar = c(10, 4, 3, 1))

bp <- barplot(
  t(as.matrix(theta)),
  col = my_colors,
  border = NA,
  space = 0,
  names.arg = sample_labels,
  las = 2,
  cex.names = 0.5,
  xlab = "",
  ylab = "Ancestry Proportions",
  main = paste("DYSTRUCT Ancestry Proportions (K =", K, ")")
)

# If era information is present, add dividers and group labels.
if (!is.null(sample_era)) {
  group_counts <- table(plot_data$era)
  group_sizes <- as.numeric(group_counts)
  group_ends <- cumsum(group_sizes)
  group_starts <- c(1, head(group_ends, -1) + 1)

  if (length(group_sizes) > 1) {
    for (i in seq_len(length(group_sizes) - 1)) {
      divider_x <- (bp[group_ends[i]] + bp[group_ends[i] + 1]) / 2
      abline(v = divider_x, col = "white", lwd = 2, lty = 2)
    }
  }

  group_centers <- mapply(function(start_idx, end_idx) mean(bp[start_idx:end_idx]), group_starts, group_ends)
  axis(side = 1,
       at = group_centers,
       labels = names(group_counts),
       tick = FALSE,
       line = 7,
       font = 2,
       cex.axis = 0.9)
}

# Add a light reference grid for readability.
abline(h = seq(0, 1, by = 0.25), col = "gray90", lty = 1, lwd = 0.7)

# Put grid behind bars by redrawing the bars after the grid would be too expensive;
# instead, make the grid subtle and keep it unobtrusive.
box()

dev.off()

cat("Success: Detected K =", K, "and saved plot to:", output_pdf, "\n")
