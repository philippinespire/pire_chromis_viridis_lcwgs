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
# Assumes Q columns are the last 3 columns of the sorted data frame
Q_sorted <- sorted_data[, (ncol(sorted_data)-2):ncol(sorted_data)]

# Count how many individuals are in each category for labeling later 
group_counts <- table(sorted_data$era)

# 4. Open PDF graphics device 
pdf(file = output_pdf, width = 10, height = 6)
	
# 5. Set margins (increased bottom margin to 5 for population labels) 
par(mar = c(5, 4, 3, 1))

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
              xlab   = "",          # Turned off default label to manually place group names
              ylab   = "Ancestry Proportions",
              main   = paste("PCAngsd Admixture Proportions (K =", K, ")"))

# 8. Add dividers and group labels below the x-axis
# Calculate the split point between groups based on individual bar widths
split_point <- group_counts[1]

# Draw a vertical line separating the categories
abline(v = split_point, col = "white", lwd = 2, lty = 2)

# Place group labels centered underneath each section
# group_counts[1] is the size of group 1; group_counts[2] is the size of group 2
axis(side = 1, 
     at = c(group_counts[1] / 2, group_counts[1] + (group_counts[2] / 2)), 
     labels = names(group_counts), 
     tick = FALSE, 
     line = 1, 
     font = 2,      # Bold text
     cex.axis = 1.2)

# 8. Close and save the PDF
dev.off()

cat("Success: Detected K =", K, "and saved plot to:", output_pdf, "\n")