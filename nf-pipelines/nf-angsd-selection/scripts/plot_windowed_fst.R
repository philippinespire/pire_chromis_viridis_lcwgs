#!/usr/bin/env Rscript

# Script to plot the windowed FST output by angsd
# Malin Pinsky, July 2026, with the assistance of Google Gemini.
# ==============================================================================
# Script: plot_windowed_fst.R
# Description: Generates a Manhattan plot from ANGSD/realSFS sliding window 
#              Fst output. Stitches chromosomes together on a continuous x-axis.
# 
# Usage: 
#   Rscript plot_windowed_fst.R <input_windows.txt> <output_plot.png>
#
# Arguments:
#   <input_windows.txt>  : The sliding window output file from realSFS.
#                          Must have a header (e.g., chr midPos Nsites Fst).
#                          The columns must be in the order: chr, midPos, Nsites, Fst.
#   <output_plot.png>    : The desired filename for the output plot. The extension
#                          dictates the file type (e.g., .png, .pdf, .jpg).
#
# Example:
#   Rscript plot_windowed_fst.R fst_windows_50k.txt fst_manhattan.png
#
# Dependencies: gplot2, dplyr
# ==============================================================================

# 1. Handle command-line arguments using base R
args <- commandArgs(trailingOnly = TRUE)

# Check if exactly two arguments are provided
if (length(args) != 2) {
  stop("Error: Incorrect number of arguments.\nUsage: Rscript plot_windowed_fst.R <input_windows.txt> <output_plot.png>", call. = FALSE)
}

input_file <- args[1]
output_plot <- args[2]

# Load required libraries silently
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))

# 2. Read the data
cat("Reading data from:", input_file, "\n")
fst_data <- read.table(input_file, header = TRUE, stringsAsFactors = FALSE, row.names = NULL)

# Rename the columns
colnames(fst_data) <- c("region","chr", "midPos", "Nsites", "Fst")

# 3. Prepare data for a continuous Manhattan plot
cat("Processing data for Manhattan plot...\n")
plot_data <- fst_data %>%
  filter(!is.na(Fst)) %>%
  mutate(midPos = as.numeric(midPos),
         Fst = as.numeric(Fst)) %>%
  group_by(chr) %>%
  summarise(chr_len = max(midPos), .groups = "drop") %>%
  mutate(tot = cumsum(as.numeric(chr_len)) - chr_len) %>%
  select(-chr_len) %>%
  left_join(fst_data, ., by = "chr") %>%
  arrange(chr, midPos) %>%
  mutate(BPcum = midPos + tot)

# 4. Calculate the center of each chromosome for x-axis labeling
axis_df <- plot_data %>%
  group_by(chr) %>%
  summarize(center = (max(BPcum) + min(BPcum)) / 2, .groups = "drop")

# 5. Create the Manhattan plot
cat("Generating plot...\n")
p <- ggplot(plot_data, aes(x = BPcum, y = Fst)) +
  geom_point(aes(color = as.factor(chr)), alpha = 0.8, size = 1.2) +
  scale_color_manual(values = rep(c("#276FBF", "#183059"), length(unique(plot_data$chr)))) +
  scale_x_continuous(label = axis_df$chr, breaks = axis_df$center) +
  scale_y_continuous(expand = expansion(mult = c(0.01, 0.05))) +
  labs(title = "Fst Sliding Window Manhattan Plot",
       x = "Chromosome",
       y = expression(F[ST])) +
  theme_classic() +
  theme(
    legend.position = "none",
    panel.border = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    axis.title = element_text(size = 12, face = "bold")
  )

# 6. Save the plot
ggsave(output_plot, plot = p, width = 12, height = 5, dpi = 300)

cat("Plot successfully saved to:", output_plot, "\n")