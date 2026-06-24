#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

historic_file <- if (length(args) >= 1) {
  args[1]
} else {
  "/archive/carpenterlab/pire/mpinsky/pire_chromis_viridis_lcwgs/nf-pipelines/nf-angsd-diversity/results/angsd_pop_theta/CviAPal_historic.pestPG"
}

modern_file <- if (length(args) >= 2) {
  args[2]
} else {
  "/archive/carpenterlab/pire/mpinsky/pire_chromis_viridis_lcwgs/nf-pipelines/nf-angsd-diversity/results/angsd_pop_theta/CviCPal_modern.pestPG"
}

output_file <- if (length(args) >= 3) {
  args[3]
} else {
  "output/tp_historic_vs_modern_mean_ci.png"
}

read_tp <- function(path, era_label) {
  dat <- read.table(
    file = path,
    header = TRUE,
    sep = "\t",
    check.names = FALSE,
    stringsAsFactors = FALSE,
    comment.char = "",
    quote = ""
  )

  if (!("tP" %in% names(dat))) {
    stop(sprintf("Column 'tP' not found in file: %s", path))
  }

  tp_vals     <- suppressWarnings(as.numeric(dat$tP))
  nsites_vals <- suppressWarnings(as.numeric(dat$nSites))
  keep <- is.finite(tp_vals) & is.finite(nsites_vals) & nsites_vals > 0
  tp_vals     <- tp_vals[keep]
  nsites_vals <- nsites_vals[keep]

  if (length(tp_vals) < 2) {
    stop(sprintf("Need at least 2 finite tP/nSites rows in file: %s", path))
  }

  data.frame(
    era    = era_label,
    tP     = tp_vals,
    nSites = nsites_vals,
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
    ci_lower = quantile(boot_stats, 0.025),
    ci_upper = quantile(boot_stats, 0.975)
  )
}

historic_df <- read_tp(historic_file, "Historical")
modern_df <- read_tp(modern_file, "Modern")
all_df <- rbind(historic_df, modern_df)

era_order <- c("Historical", "Modern")
message("Bootstrapping CIs (10 000 resamples per era) ...")
summary_list <- lapply(era_order, function(era_name) {
  rows <- all_df$era == era_name
  out  <- bootstrap_tp(all_df$tP[rows], all_df$nSites[rows])
  out$era <- era_name
  out
})
summary_df <- do.call(rbind, summary_list)
summary_df$era <- factor(summary_df$era, levels = era_order)

cols <- c("Historical" = "#1f78b4", "Modern" = "#e31a1c")
x_pos <- seq_along(era_order)
y_lim <- range(summary_df$ci_lower, summary_df$ci_upper)
y_pad <- 0.08 * diff(y_lim)
if (!is.finite(y_pad) || y_pad == 0) {
  y_pad <- 1
}

png(filename = output_file, width = 1800, height = 1200, res = 200, type = "cairo")
par(mar = c(6, 6, 4, 2) + 0.1, las = 1)

plot(
  x = x_pos,
  y = summary_df$mean_tP,
  pch = 19,
  cex = 2,
  col = cols[as.character(summary_df$era)],
  xaxt = "n",
  xlab = "Era",
  ylab = "Mean tP with 95% CI",
  ylim = c(y_lim[1] - y_pad, y_lim[2] + y_pad)
)

arrows(
  x0 = x_pos,
  y0 = summary_df$ci_lower,
  x1 = x_pos,
  y1 = summary_df$ci_upper,
  angle = 90,
  code = 3,
  length = 0.08,
  lwd = 2,
  col = cols[as.character(summary_df$era)]
)

axis(1, at = x_pos, labels = era_order)
grid(nx = NA, ny = NULL, col = "gray85", lty = "dotted")
title(main = "Historical vs Modern tP per site")

invisible(dev.off())

summary_df <- summary_df[, c("era", "n", "mean_tP", "ci_lower", "ci_upper")]
write.table(
  summary_df,
  file = stdout(),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message(sprintf("Plot written to: %s", output_file))
