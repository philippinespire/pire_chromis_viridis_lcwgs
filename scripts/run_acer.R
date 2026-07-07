#!/usr/bin/env Rscript

# Script to iteratively use Chi-sq and CMH tests from the ACER package
# to calculate genomic sites under selection and Ne.
# Sample size is calculated from the number of individuals in the ANGSD input files.
# Needs acer_helpers.r
# Original by Marianne Dehasque 2026
# Updated to handle one or more population with ANGSD input by Malin Pinsky, July 2026, 
#   with the assistance of Google Gemini.

# USE INSTRUCTIONS
# 
# On Wahab, probably need to run:
#   module load container_env R
# How to run for a single population:
# crun Rscript run_acer.R \
#  --hist_mafs=/path/to/hist1.mafs.gz \
#  --mod_mafs=/path/to/mod1.mafs.gz \
#  --region_names=Pop1 \
#  --out_dir=/path/to/output_directory \
#  --helpers=scripts/acer_helpers.R \
#  --ne_generations=114 \ # Number of generations for Ne calculation
#  --test_gen_start=0 \ # Start generation for testing
#  --test_gen_end=113 \ # End generation for testing
#  --fdr_cutoff=0.05 \ # FDR cutoff for multiple testing correction
#  --max_rounds=20 \ # Maximum number of iterations
#  --n_boot=1000 \ # Number of bootstrap replicates
#  --min_ind=4
#
# How to run for multiple populations:
# crun Rscript run_acer.R \
#  --hist_mafs=/path/hist1.mafs.gz,/path/hist2.mafs.gz \
#  --mod_mafs=/path/mod1.mafs.gz,/path/mod2.mafs.gz \
#  --region_names=Pop1,Pop2 \
#  --out_dir=/path/to/output_directory \
#  --helpers=scripts/acer_helpers.R \
#  --ne_generations=114 \ # Number of generations for Ne calculation
#  --test_gen_start=0 \ # Start generation for testing
#  --test_gen_end=113 \ # End generation for testing
#  --fdr_cutoff=0.05 \ # FDR cutoff for multiple testing correction
#  --max_rounds=20 \ # Maximum number of iterations
#  --n_boot=1000 \ # Number of bootstrap replicates
#  --min_ind=4

suppressPackageStartupMessages({
  library(ACER)
  library(boot)
  library(ggplot2)
})

# -----------------------------
# Base R Command-Line Arguments
# -----------------------------
args <- commandArgs(trailingOnly = TRUE)

# Set defaults
opt <- list(
  hist_mafs = NULL,     # Comma-separated list of historical maf files
  mod_mafs = NULL,      # Comma-separated list of modern maf files
  region_names = NULL,  # Comma-separated list of region names
  out_dir = "./results/selection",
  helpers = "acer_helpers.R",
  ne_generations = 114,
  test_gen_start = 0,
  test_gen_end = 113,
  fdr_cutoff = 0.05,
  max_rounds = 20,
  n_boot = 1000,
  min_ind = 4
)

# Parse --key=value arguments
if (length(args) > 0) {
  for (arg in args) {
    if (grepl("^--", arg) && grepl("=", arg)) {
      parts <- strsplit(sub("^--", "", arg), "=")[[1]]
      key <- parts[1]
      val <- parts[2]
      
      if (key %in% names(opt)) {
        # Cast to numeric if the default is numeric
        if (is.numeric(opt[[key]])) {
          opt[[key]] <- as.numeric(val)
        } else {
          opt[[key]] <- val
        }
      } else {
        warning(paste("Unknown argument ignored:", key))
      }
    }
  }
}

# Validate mandatory arguments
if (is.null(opt$hist_mafs) || is.null(opt$mod_mafs) || is.null(opt$region_names)) {
  stop("Error: --hist_mafs, --mod_mafs, and --region_names must all be provided.", call.=FALSE)
}

hist_files <- strsplit(opt$hist_mafs, ",")[[1]]
mod_files  <- strsplit(opt$mod_mafs, ",")[[1]]
reg_names  <- strsplit(opt$region_names, ",")[[1]]

if (length(hist_files) != length(mod_files) || length(hist_files) != length(reg_names)) {
  stop("Error: --hist_mafs, --mod_mafs, and --region_names must have the same number of comma-separated items.", call.=FALSE)
}

# -----------------------------
# Setup & Source Helpers
# -----------------------------
if (!file.exists(opt$helpers)) {
  stop(paste("Helper script not found at:", opt$helpers))
}
source(opt$helpers)

dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)
test_generations <- c(opt$test_gen_start, opt$test_gen_end)

# -----------------------------
# Input and Pre-filtering
# -----------------------------
merged_regions_list <- list()
regions <- list()

for (i in seq_along(reg_names)) {
  r_name <- reg_names[i]
  h_file <- hist_files[i]
  m_file <- mod_files[i]
  
  message(sprintf("Processing Region: %s...", r_name))
  hist_df <- read.table(gzfile(h_file), header = TRUE, stringsAsFactors = FALSE)
  mod_df  <- read.table(gzfile(m_file), header = TRUE, stringsAsFactors = FALSE)
  
  af_col_hist <- ifelse("knownEM" %in% names(hist_df), "knownEM", "freq")
  af_col_mod  <- ifelse("knownEM" %in% names(mod_df), "knownEM", "freq")
  
  merged <- merge(hist_df, mod_df, by = c("chromo", "position"), suffixes = c("_H", "_M"))
  
  # Ensure alleles match; drop totally differing alleles
  valid_sites <- (merged$major_H == merged$major_M & merged$minor_H == merged$minor_M) | 
                 (merged$major_H == merged$minor_M & merged$minor_H == merged$major_M)
  merged <- merged[valid_sites, , drop = FALSE]
  
  # Adjust modern frequencies for flipped sites
  is_flipped <- merged$major_H == merged$minor_M & merged$minor_H == merged$major_M
  merged$AF_mod_adj <- merged[[paste0(af_col_mod, "_M")]]
  merged$AF_mod_adj[is_flipped] <- 1 - merged$AF_mod_adj[is_flipped]
  
  # Format specific to this region
  reg_df <- data.frame(
    CHR = merged$chromo,
    BP  = merged$position
  )
  reg_df[[paste0(r_name, "_A_AF")]] <- merged[[paste0(af_col_hist, "_H")]]
  reg_df[[paste0(r_name, "_C_AF")]] <- merged$AF_mod_adj
  reg_df[[paste0(r_name, "_A_N")]]  <- merged$nInd_H
  reg_df[[paste0(r_name, "_C_N")]]  <- merged$nInd_M
  
  merged_regions_list[[i]] <- reg_df
  regions[[r_name]] <- c(A = paste0(r_name, "_A"), C = paste0(r_name, "_C"))
}

# Combine all regions (inner join to test only shared sites across all provided regions)
message("Combining all regions into a single dataframe...")
df_unfiltered <- merged_regions_list[[1]]
if (length(merged_regions_list) > 1) {
  for (i in 2:length(merged_regions_list)) {
    df_unfiltered <- merge(df_unfiltered, merged_regions_list[[i]], by = c("CHR", "BP"), all = FALSE)
  }
}

# Apply minimum individuals filter across ALL columns ending in _N
all_n_cols <- grep("_N$", names(df_unfiltered), value = TRUE)
df <- df_unfiltered[apply(df_unfiltered[, all_n_cols, drop = FALSE] > opt$min_ind, 1, all), , drop = FALSE]

message(sprintf("Total SNPs ready for testing (shared & passing filters): %d", nrow(df)))

# -----------------------------
# Iterative pipeline
# -----------------------------
message("Running iterative ACER selection scan...")
iterative_result <- run_iterative_selection(
  df = df,
  regions = regions,
  ne_generations = opt$ne_generations,
  test_generations = test_generations,
  fdr_cutoff = opt$fdr_cutoff,
  max_rounds = opt$max_rounds
)

iteration_summary <- build_iteration_summary(iterative_result, regions)
final_outputs <- build_final_outputs(iterative_result$final, regions)

final_test_results <- final_outputs$test_results
final_cmh_results_full <- final_outputs$cmh_results_full

# -----------------------------
# Write outputs
# -----------------------------
message("Saving results...")
write.table(iteration_summary,   file.path(opt$out_dir, "iteration_summary.tsv"),   sep = "\t", quote = FALSE, row.names = FALSE)
write.table(final_test_results,  file.path(opt$out_dir, "final_test_results.tsv"),  sep = "\t", quote = FALSE, row.names = FALSE)
write.table(final_cmh_results_full, file.path(opt$out_dir, "final_cmh_results_full.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# Bootstrap Ne CIs from the final neutral SNP set
message("Bootstrapping Ne from final neutral SNP set...")
final_selected_idx <- iterative_result$final$selected_idx
df_neutral_final <- if (length(final_selected_idx) > 0) df[-final_selected_idx, , drop = FALSE] else df
neutral_data_final <- build_region_data(df_neutral_final, regions)
ne_boot <- bootstrap_ne_by_region(
  neutral_data_final$mafs, 
  neutral_data_final$covs,
  regions, 
  generations = opt$ne_generations,
  n_boot = opt$n_boot
)

write.table(ne_boot, file.path(opt$out_dir, "ne_bootstrap.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

message("\n--- ACER Summary ---")
message("Converged after ", iterative_result$n_rounds, " rounds.")
message("Total SNPs tested: ", nrow(iterative_result$final$df))
message("Total SNPs under selection: ", iterative_result$n_total_removed)
message("Neutral SNPs remaining: ", iterative_result$n_remaining)
message("--------------------\n")

# -----------------------------
# Chi-squared Manhattan plots
# -----------------------------
chisq_dat <- final_test_results
if (nrow(chisq_dat) > 0 && all(c("CHR", "BP") %in% names(chisq_dat))) {
  chisq_dat$CHR <- as.character(chisq_dat$CHR)
  chisq_dat$BP <- as.numeric(chisq_dat$BP)
  chisq_dat <- chisq_dat[!is.na(chisq_dat$BP), , drop = FALSE]
  
  if (nrow(chisq_dat) > 0) {
    chr_levels <- unique(chisq_dat$CHR)
    chisq_dat$CHR <- factor(chisq_dat$CHR, levels = chr_levels)
    
    chr_max <- aggregate(BP ~ CHR, chisq_dat, max)
    chr_max <- chr_max[match(levels(chisq_dat$CHR), chr_max$CHR), ]
    chr_max$cumstart <- c(0, cumsum(head(chr_max$BP, -1)))
    
    chisq_dat$cumBP <- chisq_dat$BP + chr_max$cumstart[match(chisq_dat$CHR, chr_max$CHR)]
    
    for (r_name in names(regions)) {
      pval_col <- paste0(r_name, "_chisq_pval")
      fdr_col  <- paste0(r_name, "_chisq_fdr")
      
      if (pval_col %in% names(chisq_dat)) {
        # Filter out NA and 0 p-values for log transformation
        plt_dat <- chisq_dat[!is.na(chisq_dat[[pval_col]]) & chisq_dat[[pval_col]] > 0, ]
        
        if (nrow(plt_dat) > 0) {
          plt_dat$logp <- -log10(plt_dat[[pval_col]])
          plt_dat$is_fdr_sig <- !is.na(plt_dat[[fdr_col]]) & plt_dat[[fdr_col]] < opt$fdr_cutoff
          
          p <- ggplot(plt_dat, aes(x = cumBP, y = logp)) +
            geom_point(color = "grey50", alpha = 0.6, size = 0.6) +
            geom_point(data = plt_dat[plt_dat$is_fdr_sig, , drop = FALSE],
                       color = "blue", alpha = 0.9, size = 0.8) +
            scale_x_continuous(labels = NULL, breaks = NULL) +
            labs(
              x = "Genomic position",
              y = "-log10(Chi-sq p-value)",
              title = paste("Chi-squared Manhattan Plot (converged) -", r_name),
              subtitle = paste0("Blue: FDR < ", opt$fdr_cutoff)
            ) +
            theme_bw() +
            theme(legend.position = "none")
          
          out_name <- paste0("chisq_manhattan_", r_name, "_final.png")
          ggsave(file.path(opt$out_dir, out_name), p, width = 12, height = 5, dpi = 300)
          message(paste("Chi-sq Manhattan plot for", r_name, "saved successfully."))
        }
      }
    }
  }
}

# -----------------------------
# CMH Manhattan plot 
# -----------------------------
if (length(regions) > 1) {
  dat <- final_cmh_results_full
  dat <- dat[!is.na(dat$cmh_pval) & dat$cmh_pval > 0, , drop = FALSE]
  
  if (nrow(dat) > 0 && all(c("CHR", "BP") %in% names(dat))) {
    dat$CHR <- as.character(dat$CHR)
    dat$BP <- as.numeric(dat$BP)
    dat <- dat[!is.na(dat$BP), , drop = FALSE]
  
    if (nrow(dat) > 0) {
      chr_levels <- unique(dat$CHR)
      dat$CHR <- factor(dat$CHR, levels = chr_levels)
  
      chr_max <- aggregate(BP ~ CHR, dat, max)
      chr_max <- chr_max[match(levels(dat$CHR), chr_max$CHR), ]
      chr_max$cumstart <- c(0, cumsum(head(chr_max$BP, -1)))
  
      dat$cumBP <- dat$BP + chr_max$cumstart[match(dat$CHR, chr_max$CHR)]
      dat$logp <- -log10(dat$cmh_pval)
      dat$is_fdr_sig <- !is.na(dat$cmh_fdr) & dat$cmh_fdr < opt$fdr_cutoff
  
      p <- ggplot(dat, aes(x = cumBP, y = logp)) +
        geom_point(color = "grey50", alpha = 0.6, size = 0.6) +
        geom_point(data = dat[dat$is_fdr_sig, , drop = FALSE],
                   color = "red", alpha = 0.9, size = 0.8) +
        scale_x_continuous(labels = NULL, breaks = NULL) +
        labs(
          x = "Genomic position",
          y = "-log10(CMH p-value)",
          title = "CMH Manhattan Plot (converged)",
          subtitle = paste0("Red: FDR < ", opt$fdr_cutoff)
        ) +
        theme_bw() +
        theme(legend.position = "none")
  
      ggsave(file.path(opt$out_dir, "cmh_manhattan_final.png"), p, width = 12, height = 5, dpi = 300)
      message("Manhattan plot saved successfully.")
    }
  }
} else {
  message("Only one region provided; skipping CMH Manhattan plot.")
}