#!/bin/bash

# ==============================================================================
# SCRIPT NAME: plot_softclip_vs_depth.sbatch
# DESCRIPTION: Extract and merge soft-clipped base statistics with mean 
#              read depth across individual samples, generating a scatter plot.
#
# Malin Pinsky, July 2026, with help from Gemini
#
# INPUT REQUIREMENTS:
#   1. Soft-clip Stats File:
#      - Location: output/softclip_analysis-iridian/per_individual_stats.txt
#      - Format: Tab-delimited (Sample, SoftClippedBases, TotalBases, Fraction)
#      - Sample name format example: CviAPal001.merged.L98.realn
#
#   2. Individual Depth Files
#      - Location: nf-pipelines/nf-trim-generode/results-iridian/depth/*.dpstats.txt
#      - Format: Single mean depth number on the first line of each file
#      - Sample format example: CviAPal001.Q25.L103.bam.dpstats.txt
#
# SAMPLE MATCHING LOGIC:
#   - Both datasets are matched on the prefix before the first dot (`.`)
#     e.g., "CviAPal001.merged.L98.realn" and "CviAPal001.Q25.L103.bam.dpstats.txt"
#     both parse to "CviAPal001".
#
# ENVIRONMENT DEPENDENCIES:
#   - Module: container_env R
#   - Execution Wrapper: crun
#   - R Libraries: ggplot2, dplyr (must be installed inside the container)
#
# HOW TO RUN:
#   Submit to Slurm scheduler:
#     sbatch plot_softclip_vs_depth.sbatch
#
# OUTPUTS:
#   - Plot Output:  output/softclip_analysis-iridiansoftclip_vs_depth_plot.png (300 DPI PNG)
#   - Log Outputs:  logs/softclip_vs_depth_<JOBID>.log
#                   logs/softclip_vs_depth_<JOBID>.err
# ==============================================================================

#SBATCH --job-name=softclip_vs_depth
#SBATCH --output=logs/softclip_vs_depth_%j.log
#SBATCH --error=logs/softclip_vs_depth_%j.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G

# Load the required containerized R environment
module load container_env R

# ------------------------------------------------------------------------------
# 1. Write the R script to a temporary file
# ------------------------------------------------------------------------------
cat << 'EOF' > plot_softclip_depth.R
library(ggplot2)
library(dplyr)

# Paths
softclip_path <- "output/softclip_analysis-iridian/per_individual_stats.txt"
depth_dir     <- "nf-pipelines/nf-trim-generode/results-iridian/depth"
output_plot   <- "output/softclip_analysis-iridian/softclip_vs_depth_plot.png"

# Load soft-clip data
cat("Reading soft-clip analysis stats...\n")
softclip_df <- read.table(softclip_path, header = TRUE, stringsAsFactors = FALSE)

# Load individual depth files
cat("Reading depth files...\n")
depth_files <- list.files(depth_dir, pattern = "\\.dpstats\\.txt$", full.names = TRUE)

depth_list <- lapply(depth_files, function(file_path) {
  fname <- basename(file_path)
  # Read single average depth number
  val <- as.numeric(readLines(file_path, n = 1, warn = FALSE))
  
  data.frame(
    DepthFile = fname,
    MeanDepth = val,
    stringsAsFactors = FALSE
  )
})

depth_df <- do.call(rbind, depth_list)

# Extract a common Sample ID key for merging
# Adjust regex if your depth filenames use a different naming scheme
softclip_df$SampleID <- gsub("\\..*", "", softclip_df$Sample)
depth_df$SampleID    <- gsub("\\..*", "", depth_df$DepthFile)

# Merge datasets
merged_df <- inner_join(softclip_df, depth_df, by = "SampleID")

cat(sprintf("Successfully matched %d samples.\n", nrow(merged_df)))

# Generate plot
p <- ggplot(merged_df, aes(x = MeanDepth, y = Fraction)) +
  geom_point(color = "#2b5c8f", size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, color = "#d95f02", linetype = "dashed") +
  theme_light(base_size = 14) +
  labs(
    title = "Soft-clipped Base Fraction vs. Mean Read Depth",
    subtitle = paste("n =", nrow(merged_df), "individuals"),
    x = "Mean Read Depth (X)",
    y = "Fraction of Soft-Clipped Bases"
  ) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

# Save plot output
ggsave(output_plot, plot = p, width = 8, height = 6, dpi = 300)
cat(paste("Plot successfully saved to:", output_plot, "\n"))

EOF

# ------------------------------------------------------------------------------
# 2. Execute the R script via crun
# ------------------------------------------------------------------------------
crun Rscript plot_softclip_depth.R

# Clean up the script file afterward
rm plot_softclip_depth.R