#!/usr/bin/env Rscript

# ==============================================================================
# SCRIPT NAME:  pcangsd.R
# PURPOSE:      Calculates variance explained and generates 2-panel PCA plots 
#               (PC1 vs PC2 and PC2 vs PC3) from a PCAngsd covariance matrix, 
#               color-coded by population and shape-coded by sampling era.
#
# USAGE:        Rscript pcangsd.R --cov <input.cov> \
#                                 --samplesheet <samplesheet.csv> \
#                                 --out <output.pdf> \
#                                 [--species <species_name>]
#
# INPUT ARGUMENTS:
#   --cov         Path to space-separated PCAngsd output covariance matrix (.cov)
#   --samplesheet CSV file containing metadata. Must include 'sample', 'pop', and 
#                 'era' (historic/modern) columns matching matrix row order.
#   --out         Target PDF file path for saved PCA plots
#   --species     [Optional] Name of target species for main plot title 
#                 (Default: "Species")
#
# OUTPUTS:       PDF containing a two-panel scatter plot (PC1/PC2 & PC2/PC3).
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)

# simple args parser: --key value pairs
alist <- list()
for (i in seq(1, length(args), by = 2)) alist[[sub("^--","", args[i])]] <- args[i+1]

cov_path       <- alist$cov
samplesheet    <- alist$samplesheet
species        <- if (!is.null(alist$species)) alist$species else "Species"
out_pdf        <- alist$out

# read samplesheet; expect columns: sample, pop, bam, era (historic|modern)
meta <- read.csv(samplesheet, header = TRUE)

C <- as.matrix(read.table(cov_path))
e <- eigen(C)
e_value<-e$values

pc1_variance<-e_value[1]/sum(e_value)*100
pc2_variance<-e_value[2]/sum(e_value)*100
pc3_variance<-e_value[3]/sum(e_value)*100

# Assign colours by population

colours_colorblind <- c(
  "#332288", "#88CCEE", "#44AA99", "#117733",
  "#DDCC77", "#CC6677", "#AA4499", "#882255",
  "#999933", "#661100"
)

uniq_pops <- unique(meta$pop)
n_pops <- length(uniq_pops)
pop_palette <- setNames(colours_colorblind[1:n_pops], uniq_pops)
col_vec <- pop_palette[meta$pop]

# Assign shapes by era (historic vs modern)
shape_map <- c(historic = 24, modern = 21)
pch_vec <- shape_map[meta$era]


pdf(out_pdf, width = 10, height = 5)
layout(matrix(1:2, 1, 2, byrow = TRUE), respect = TRUE)
par(oma = c(0, 0, 3, 0))
par(pty = "s")
par(mar = c(4, 5, 1, 1))

plot(e$vectors[, 1:2],
     bg = col_vec,
     pch = pch_vec,
     cex = 1.2,
     xlab = paste("PC1 ", round(pc1_variance, 2), "%"),
     ylab = paste("PC2 ", round(pc2_variance, 2), "%"))

legend("topleft",
       legend = uniq_pops,
       pt.bg = pop_palette[uniq_pops],
       pch = 21,
       cex = 0.8)

plot(e$vectors[, 2:3],
     bg = col_vec,
     pch = pch_vec,
     cex = 1.2,
     xlab = paste("PC2 ", round(pc2_variance, 2), "%"),
     ylab = paste("PC3 ", round(pc3_variance, 2), "%"))

mtext(paste0("PCAngsd - ", species), outer = TRUE, cex = 1, font = 2, line = 1)
dev.off()