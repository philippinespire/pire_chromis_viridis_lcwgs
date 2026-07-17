#!/usr/bin/env Rscript

# ==============================================================================
# SCRIPT NAME:  plot_admixture.R
# PURPOSE:      Visualizes PCAngsd Admixture proportions (K=3) into a high-quality 
#               PDF barplot, sorted and grouped by historical vs. modern metadata.
#
# USAGE:        Rscript plot_admixture.R <input.Q> <metadata.csv> <output.pdf>
#
# INPUTS:       1. <input.Q>    - Space-separated PCAngsd admixture matrix.
#               2. <metadata.csv> - CSV file containing sample metadata. Must 
#                                 include a "era" column (e.g., historic/modern).
#                                 Row order must match the .Q file exactly.
#               3. <output.pdf>  - Target filename for the output PDF.
#
# REQUIREMENTS: Base R installation.
# DATE:         June 2026
# AUTHOR:       Malin Pinsky with assistance from GitHub Copilot with the GPT-5.3-Codex model.
# ==============================================================================


# 1. Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Validate that all 3 required arguments are present
if (length(args) < 3) {
  stop("\nUsage: Rscript plot_admixture.R <input.Q> <metadata.csv> <output.pdf>\n", call. = FALSE)
}

q_file     <- args[1]
meta_file  <- args[2]
output_pdf <- args[3]

# 2. Read data
if (!file.exists(q_file)) stop(paste("Error: File not found -", q_file))
if (!file.exists(meta_file)) stop(paste("Error: File not found -", meta_file))

Q    <- read.table(q_file)
meta <- read.csv(meta_file)

# Automatically infer K from the number of columns in the .Q file
K <- ncol(Q)

# Check for matching row counts
if (nrow(Q) != nrow(meta)) {
  stop(paste("Error: Row mismatch! .Q file has", nrow(Q), "rows, but CSV has", nrow(meta), "rows."))
}

# 3. Combine and sort data by the era column
# This groups all 'historic' individuals together and all 'modern' individuals together
combined <- cbind(meta, Q)

if (!("era" %in% colnames(meta))) {
  stop(paste("Error: Column 'era' not found in metadata file.", sep=""))
}

sorted_data <- combined[order(combined$era), ]

# Separate back out into the sorted Q matrix
# Assumes Q columns are the last K columns of the sorted data frame
Q_sorted <- sorted_data[, (ncol(sorted_data)-K+1):ncol(sorted_data)]

# Choose sample labels for the x-axis; prefer common sample ID column names.
sample_label_candidates <- c("sample", "sample_id", "sampleid", "individual", "ind", "id", "name")
meta_colnames_lower <- tolower(colnames(sorted_data))
matched_label_col <- sample_label_candidates[sample_label_candidates %in% meta_colnames_lower]

if (length(matched_label_col) > 0) {
  sample_label_col <- colnames(sorted_data)[which(meta_colnames_lower == matched_label_col[1])[1]]
} else {
  sample_label_col <- colnames(meta)[1]
}

sample_labels <- as.character(sorted_data[[sample_label_col]])

# Count how many individuals are in each category for labeling later 
group_counts <- table(sorted_data$era)

# 4. Open PDF graphics device 
pdf(file = output_pdf, width = 10, height = 6)
	
# 5. Set margins (increase bottom margin for sample and era labels)
par(mar = c(10, 4, 3, 1))

# 6. Generate color palette dynamically based on the inferred K
if (K <= 12) {
  my_colors <- palette.colors(n = K, palette = "Set 3")
} else {
  my_colors <- rainbow(K)
}

# 7. Generate the barplot
bp <- barplot(t(Q_sorted), 
              col    = my_colors, 
              border = NA, 
              space  = 0,
              names.arg = sample_labels,
              las = 2,
              cex.names = 0.5,
              xlab   = "",          # Turned off default label to manually place group names
              ylab   = "Ancestry Proportions",
              main   = paste("PCAngsd Admixture Proportions (K =", K, ")"))

# 8. Add dividers and group labels below the x-axis
# Draw vertical lines between era groups and place centered group labels
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

# 8. Close and save the PDF
dev.off()

cat("Success: Detected K =", K, "and saved plot to:", output_pdf, "\n")