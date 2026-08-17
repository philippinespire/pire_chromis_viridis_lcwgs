#!/usr/bin/env Rscript

# DATE:         August 2026
# AUTHOR:       Malin Pinsky with assistance from Google Gemini
# ==============================================================================
# Script: plot_tp.R
# Description: Calculates mean nucleotide diversity (pi / tP per site) and 95% 
#              bootstrap confidence intervals across historical and modern 
#              timepoints for a given region. Outputs a summary TSV table and 
#              a comparison plot.
# 
# Usage: 
#   Rscript plot_tp.R --historic_file <hist.pestPG> --modern_file <mod.pestPG> [options]
#
# Required Arguments:
#   --historic_file <path> : Path to the historic ANGSD .pestPG theta file
#   --modern_file <path>   : Path to the modern ANGSD .pestPG theta file
#
# Optional Arguments:
#   --region <string>      : Region or population label (Default: "Region")
#   --out_png <path>       : Output path for PNG plot (Default: "pi_historic_vs_modern.png")
#   --out_tsv <path>       : Output path for summary table (Default: "pi_summary.tsv")
#   --n_boot <int>         : Number of bootstrap iterations for CIs (Default: 10000)
#
# Example:
#   Rscript plot_tp.R \
#     --historic_file=Bali_historic.pestPG \
#     --modern_file=Bali_modern.pestPG \
#     --region=Bali \
#     --out_png=pi_historic_vs_modern_Bali.png \
#     --out_tsv=pi_summary_Bali.tsv
#
# Dependencies: boot
# ==============================================================================

# Dependencies
suppressPackageStartupMessages({
  library(boot)
})

args <- commandArgs(trailingOnly = TRUE)

opt <- list(
  historic_file = NULL,
  modern_file   = NULL,
  region        = "Region",
  out_png       = "pi_historic_vs_modern.png",
  out_tsv       = "pi_summary.tsv",
  n_boot        = 10000
)

# Parse command line options
if (length(args) > 0) {
  i <- 1
  while (i <= length(args)) {
    arg <- args[i]
    if (grepl("^--", arg)) {
      parts <- if (grepl("=", arg)) strsplit(sub("^--", "", arg), "=")[[1]] else c(sub("^--", "", arg), args[i + 1])
      key <- parts[1]
      val <- parts[2]
      i <- if (grepl("=", arg)) i + 1 else i + 2
      if (key %in% names(opt)) {
        opt[[key]] <- if (is.numeric(opt[[key]])) as.numeric(val) else val
      }
    } else { i <- i + 1 }
  }
}

if (is.null(opt$historic_file) || is.null(opt$modern_file)) {
  stop("Error: Must provide both --historic_file and --modern_file.", call. = FALSE)
}

read_tp <- function(path, era_label) {
  dat <- read.table(
    file = path, header = TRUE, sep = "\t", check.names = FALSE,
    stringsAsFactors = FALSE, comment.char = "", quote = ""
  )

  if (!("tP" %in% names(dat))) {
    stop(sprintf("Column 'tP' not found in file: %s", path))
  }

  tp_vals     <- suppressWarnings(as.numeric(dat$tP))
  nsites_vals <- suppressWarnings(as.numeric(dat$nSites))
  keep        <- is.finite(tp_vals) & is.finite(nsites_vals) & nsites_vals > 0
  
  data.frame(
    era    = era_label,
    tP     = tp_vals[keep],
    nSites = nsites_vals[keep],
    stringsAsFactors = FALSE
  )
}

bootstrap_tp <- function(tp_vals, nsites_vals, n_boot = 10000, seed = 42) {
  set.seed(seed)
  n <- length(tp_vals)
  obs_tp <- sum(tp_vals) / sum(nsites_vals)
  
  boot_stats <- vapply(seq_len(n_boot), function(i) {
    idx <- sample.int(n, n, replace = TRUE)
    sum(tp_vals[idx]) / sum(nsites_vals[idx])
  }, numeric(1))
  
  data.frame(
    n        = n,
    mean_tP  = obs_tp,
    ci_lower = as.numeric(quantile(boot_stats, 0.025)),
    ci_upper = as.numeric(quantile(boot_stats, 0.975))
  )
}

historic_df <- read_tp(opt$historic_file, "Historical")
modern_df   <- read_tp(opt$modern_file, "Modern")
all_df      <- rbind(historic_df, modern_df)

era_order <- c("Historical", "Modern")
summary_list <- lapply(era_order, function(era_name) {
  rows <- all_df$era == era_name
  out  <- bootstrap_tp(all_df$tP[rows], all_df$nSites[rows], n_boot = opt$n_boot)
  out$era <- era_name
  out
})

summary_df     <- do.call(rbind, summary_list)
summary_df$era <- factor(summary_df$era, levels = era_order)

# Save Summary TSV
out_summary <- summary_df[, c("era", "n", "mean_tP", "ci_lower", "ci_upper")]
write.table(out_summary, file = opt$out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)

# Generate Plot
cols  <- c("Historical" = "#1f78b4", "Modern" = "#e31a1c")
x_pos <- seq_along(era_order)
y_lim <- range(summary_df$ci_lower, summary_df$ci_upper)
y_pad <- if (diff(y_lim) == 0) 0.01 else 0.08 * diff(y_lim)

png(filename = opt$out_png, width = 1800, height = 1200, res = 200, type = "cairo")
par(mar = c(6, 6, 4, 2) + 0.1, las = 1)

plot(
  x = x_pos,
  y = summary_df$mean_tP,
  pch = 19, cex = 2,
  col = cols[as.character(summary_df$era)],
  xaxt = "n", xlab = "Era", ylab = "Mean pi (tP per site)",
  ylim = c(y_lim[1] - y_pad, y_lim[2] + y_pad)
)

arrows(
  x0 = x_pos, y0 = summary_df$ci_lower,
  x1 = x_pos, y1 = summary_df$ci_upper,
  angle = 90, code = 3, length = 0.08, lwd = 2,
  col = cols[as.character(summary_df$era)]
)

axis(1, at = x_pos, labels = era_order)
grid(nx = NA, ny = NULL, col = "gray85", lty = "dotted")
title(main = paste("Nucleotide Diversity (pi) -", opt$region))

invisible(dev.off())
message(sprintf("Plot saved to %s", opt$out_png))