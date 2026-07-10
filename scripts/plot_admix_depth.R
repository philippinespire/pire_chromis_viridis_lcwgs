#!/usr/bin/env Rscript

# Script to plot admixture proportions vs. read depths
# Reads read depth from the *.dpstats.txt file output by the nf-pipelines BAM_QC process (samtools)
# Read admixture proportions from a PCAngsd .Q file and sample names from the angsd bamlist.txt file.
#
# Requirements:
#  ggplot2
#
# Run as:
# Rscript plot_admix_depth.R <path/to/bamlist.txt> <path/to/.Q> <path/to/dpstats_dir> <column_index_for_admixture> <output.png>
#
# Malin Pinsky, 2026
# Developed with the assistance of Google Gemini (July 2026).



# Load required library for plotting
if (!require(ggplot2, quietly = TRUE)) {
  stop("The 'ggplot2' package is required but not installed.")
}

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  cat("Usage: Rscript plot_admix_depth.R <bamlist.txt> <.Q_file> <dpstats_dir> <admix_column> <output_plot.png>\n")
  cat("Example: Rscript plot_admix_depth.R bamlist.txt Cvi.pcangsd.cvi-only.admix.2.Q ./dpstats 1 output.png\n")
  quit(status = 1)
}

bamlist_file <- args[1]
q_file       <- args[2]
dpstats_dir  <- args[3]
k_col        <- as.integer(args[4])
output_plot  <- args[5]

# 1. Read bamlist and extract individual names
# Extracts everything before the first period (e.g., "CviAPal001.merged..." -> "CviAPal001")
cat("Reading bamlist\n")
bamlist_paths <- readLines(bamlist_file, warn = FALSE)
inds_bam <- sub("\\..*", "", basename(bamlist_paths))

# 2. Read Q file and attach individual names
cat("Reading Q file\n")
q_data <- read.table(q_file, header = FALSE)
if (k_col < 1 || k_col > ncol(q_data)) {
  stop(paste("Error: admix_column must be between 1 and", ncol(q_data)))
}

# Create a dataframe for admixture proportions
admix_df <- data.frame(
  Ind = inds_bam,
  AdmixProp = q_data[, k_col],
  stringsAsFactors = FALSE
)

# 3. Read depth stats from the directory
cat("Reading dpstats\n")
dp_files <- list.files(dpstats_dir, pattern = "\\.dpstats\\.txt$", full.names = TRUE)
if (length(dp_files) == 0) {
  stop("No .dpstats.txt files found in the specified directory.")
}

dp_data <- data.frame(Ind = character(), Depth = numeric(), stringsAsFactors = FALSE)

for (f in dp_files) {
  # Read the single depth value from the file
  val <- as.numeric(readLines(f, warn = FALSE)[1])
  # Extract individual ID from the filename (e.g., "CviAPal001.Q25..." -> "CviAPal001")
  ind <- sub("\\..*", "", basename(f))
  dp_data <- rbind(dp_data, data.frame(Ind = ind, Depth = val, stringsAsFactors = FALSE))
}

# 4. Merge the datasets by individual ID
cat("Merging datasets\n")
merged_data <- merge(admix_df, dp_data, by = "Ind")

if (nrow(merged_data) == 0) {
  stop("Merge failed: No matching individual IDs between the bamlist/.Q file and the .dpstats.txt files.")
}

# Extract the group identifier (4th character: A or C from CviA or CviC)
merged_data$Group <- substr(merged_data$Ind, 4, 4)
merged_data$Group <- factor(merged_data$Group, levels = c("A", "C"))

# 5. Generate the Plot
cat("Plotting\n")
p <- ggplot(merged_data, aes(y = AdmixProp, x = Depth, color = Group)) +
  geom_point(alpha = 0.7, size = 3) +
  scale_color_manual(values = c("A" = "#1f77b4", "C" = "#d62728"), drop = FALSE) + # color by group
  theme_light() +
  labs(
    title = paste("Admixture Proportion vs. Sequencing Depth"),
    subtitle = paste("Plotting Cluster", k_col),
    y = paste("Admixture Proportion (Cluster", k_col, ")"),
    x = "Sequencing Depth (dpstats)",
    color = "Era"
  )

# 6. Save the plot
ggsave(output_plot, plot = p, width = 7, height = 5, dpi = 300)
cat("Plot successfully saved to:", output_plot, "\n")
