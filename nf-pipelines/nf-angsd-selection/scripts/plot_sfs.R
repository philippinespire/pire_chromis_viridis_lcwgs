#!/usr/bin/env Rscript

library(ggplot2)
library(dplyr)
library(tidyr)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript plot_sfs.R <output_png> <2dsfs_file1> [2dsfs_file2 ...]")
}

out_png   <- args[1]
sfs_files <- args[-1]

# Function to parse 2D SFS vector and extract 1D Historical and Modern spectra
parse_2dsfs <- function(file_path) {
  if (!file.exists(file_path) || file.size(file_path) == 0) return(NULL)
  
  sfs_2d <- scan(file_path, quiet = TRUE)
  L <- length(sfs_2d)
  if (L <= 1) return(NULL)
  
  # Determine matrix dimensions dim_A (Historical) and dim_C (Modern)
  sqrt_L <- round(sqrt(L))
  if (sqrt_L * sqrt_L == L) {
    dim_A <- sqrt_L
    dim_C <- sqrt_L
  } else {
    odds <- seq(3, floor(sqrt(L)), by = 2)
    odds <- odds[L %% odds == 0]
    odds <- odds[(L / odds) %% 2 == 1]
    if (length(odds) == 0) return(NULL)
    dim_A <- tail(odds, 1)
    dim_C <- L / dim_A
  }
  
  # Reshape to matrix (ANGSD writes in C row-major order: row = Hist, col = Modern)
  mat <- matrix(sfs_2d, nrow = dim_A, ncol = dim_C, byrow = TRUE)
  hist_sfs <- rowSums(mat)
  mod_sfs  <- colSums(mat)
  
  if (length(hist_sfs) <= 1 || length(mod_sfs) <= 1) return(NULL)
  
  # Exclude monomorphic bin (index 1) and calculate relative proportions
  hist_poly <- hist_sfs[-1]
  mod_poly  <- mod_sfs[-1]
  
  hist_prop <- hist_poly / sum(hist_poly)
  mod_prop  <- mod_poly / sum(mod_poly)

  # Density scaling
  hist_density <- hist_prop * length(hist_prop)
  mod_density  <- mod_prop  * length(mod_prop)
  
  # Relative allele frequency bins (0 to 1)
  hist_freq <- seq_along(hist_prop) / length(hist_prop)
  mod_freq  <- seq_along(mod_prop)  / length(mod_prop)
  
  region_name <- gsub("_hist_vs_mod\\.2dsfs$", "", basename(file_path))
  region_name <- gsub("_", " ", region_name)
  
  df_hist <- data.frame(
    Frequency  = hist_freq,
    Density    = hist_density,
    Timepoint  = "Historical",
    Title      = region_name,
    stringsAsFactors = FALSE
  )
  
  df_mod <- data.frame(
    Frequency  = mod_freq,
    Density    = mod_density,
    Timepoint  = "Modern",
    Title      = region_name,
    stringsAsFactors = FALSE
  )
  
  rbind(df_hist, df_mod)
}

# Process all input 2D SFS files
data_list <- lapply(sfs_files, parse_2dsfs)
df_combined <- bind_rows(data_list)

if (is.null(df_combined) || nrow(df_combined) == 0) {
  stop("No valid 2D SFS data found.")
}

num_panels <- length(unique(df_combined$Title))
num_cols   <- min(4, num_panels)
num_rows   <- ceiling(num_panels / num_cols)

p <- ggplot(df_combined, aes(x = Frequency, y = Density, color = Timepoint, linetype = Timepoint)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 0.5) +
  scale_color_manual(values = c("Historical" = "#1f77b4", "Modern" = "#d62728")) +
  scale_linetype_manual(values = c("Historical" = "solid", "Modern" = "solid")) +
  facet_wrap(~ Title, scales = "free", ncol = num_cols) +
  labs(
    title = "Site Frequency Spectrum Comparison: Neutral SNPs",
    x = "Derived Allele Frequency (excl. fixed sites)",
    y = "Probability Density (Area = 1.0)",
    color = "Timepoint",
    linetype = "Timepoint"
  ) +
  theme_bw(base_size = 10) +
  theme(
    legend.position = "top",
    legend.background = element_rect(color = "grey80", fill = "white"),
    legend.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 13),
    strip.background = element_rect(fill = "grey90"),
    strip.text = element_text(face = "bold", size = 8)
  )

plot_width  <- max(6, num_cols * 3.5)
plot_height <- max(4, num_rows * 2.8)

ggsave(out_png, p, width = plot_width, height = plot_height, dpi = 300)
message("✓ Saved plot: ", out_png)