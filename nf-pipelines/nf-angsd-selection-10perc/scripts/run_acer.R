#!/usr/bin/env Rscript

# ==============================================================================
# SCRIPT NAME:  run_acer.R
# PURPOSE:      Executes an iterative genome-wide selection scan and temporal Ne 
#               estimation using chi-squared and CMH tests from the ACER package 
#               across historical and modern population timepoints.
#
# USAGE:        Rscript run_acer.R --hist_mafs <hist1.mafs.gz,hist2.mafs.gz> \
#                                  --mod_mafs <mod1.mafs.gz,mod2.mafs.gz> \
#                                  --region_names <reg1,reg2> [options]
#
# REQUIRED ARGUMENTS:
#   --hist_mafs     Comma-separated list of historical ANGSD .mafs.gz file paths
#   --mod_mafs      Comma-separated list of modern ANGSD .mafs.gz file paths
#   --region_names  Comma-separated list of region/population names corresponding
#                   to the order of the MAF files
#
# OPTIONAL ARGUMENTS:
#   --out_dir       Output results directory (Default: ./results/selection)
#   --helpers       Path to acer_helpers.R script (Default: scripts/acer_helpers.R)
#   --generations   Elapsed generations between sampling points (Default: 114)
#   --fdr_cutoff    False discovery rate threshold for candidate loci (Default: 0.05)
#   --max_rounds    Maximum iterations for Ne/selection scan convergence (Default: 20)
#   --n_boot        Bootstrap iterations for Ne confidence intervals (Default: 1000)
#   --min_ind       Minimum individual count threshold per site (Default: 4)
#
# OUTPUTS:       tsv files for iteration summary, final test results, CMH 
#               statistics, and Ne bootstrap confidence intervals.
# ==============================================================================

suppressPackageStartupMessages({
  library(ACER)
  library(boot)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)

opt <- list(
  hist_mafs      = NULL,
  mod_mafs       = NULL,
  region_names   = NULL,
  out_dir        = "./results/selection",
  helpers        = "scripts/acer_helpers.R",
  generations    = 114,
  fdr_cutoff     = 0.05,
  max_rounds     = 20,
  n_boot         = 1000,
  min_ind        = 4
)

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

if (is.null(opt$hist_mafs) || is.null(opt$mod_mafs) || is.null(opt$region_names)) {
  stop("Error: Must provide --hist_mafs, --mod_mafs, and --region_names.", call.=FALSE)
}

if (!file.exists(opt$helpers)) stop(paste("Helper script not found at:", opt$helpers))
source(opt$helpers)

dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)

# Automatically derive ne_generations and test_generations
ne_generations   <- opt$generations
test_generations <- c(0, opt$generations - 1)

hist_files <- strsplit(opt$hist_mafs, ",")[[1]]
mod_files  <- strsplit(opt$mod_mafs, ",")[[1]]
reg_names  <- strsplit(opt$region_names, ",")[[1]]

merged_regions_list <- list()
regions <- list()

for (i in seq_along(reg_names)) {
  r_name <- reg_names[i]
  hist_df <- read.table(gzfile(hist_files[i]), header = TRUE, stringsAsFactors = FALSE)
  mod_df  <- read.table(gzfile(mod_files[i]), header = TRUE, stringsAsFactors = FALSE)
  
  af_col_hist <- ifelse("knownEM" %in% names(hist_df), "knownEM", "freq")
  af_col_mod  <- ifelse("knownEM" %in% names(mod_df), "knownEM", "freq")
  
  merged <- merge(hist_df, mod_df, by = c("chromo", "position"), suffixes = c("_H", "_M"))
  
  # Ensure site allele compatibility
  valid_sites <- (merged$major_H == merged$major_M & merged$minor_H == merged$minor_M) | 
                 (merged$major_H == merged$minor_M & merged$minor_H == merged$major_M)
  merged <- merged[valid_sites, , drop = FALSE]
  
  is_flipped <- merged$major_H == merged$minor_M & merged$minor_H == merged$major_M
  merged$AF_mod_adj <- merged[[paste0(af_col_mod, "_M")]]
  merged$AF_mod_adj[is_flipped] <- 1 - merged$AF_mod_adj[is_flipped]
  
  reg_df <- data.frame(CHR = merged$chromo, BP = merged$position)
  reg_df[[paste0(r_name, "_A_AF")]] <- merged[[paste0(af_col_hist, "_H")]]
  reg_df[[paste0(r_name, "_C_AF")]] <- merged$AF_mod_adj
  reg_df[[paste0(r_name, "_A_N")]]  <- merged$nInd_H
  reg_df[[paste0(r_name, "_C_N")]]  <- merged$nInd_M
  
  merged_regions_list[[i]] <- reg_df
  regions[[r_name]] <- c(A = paste0(r_name, "_A"), C = paste0(r_name, "_C"))
}

# Merge and filter across regions
df_unfiltered <- merged_regions_list[[1]]
if (length(merged_regions_list) > 1) {
  for (i in 2:length(merged_regions_list)) {
    df_unfiltered <- merge(df_unfiltered, merged_regions_list[[i]], by = c("CHR", "BP"), all = FALSE)
  }
}

af_cols <- grep("_AF$", names(df_unfiltered), value = TRUE)
df_unfiltered <- df_unfiltered[complete.cases(df_unfiltered[, af_cols]), , drop = FALSE]

all_n_cols <- grep("_N$", names(df_unfiltered), value = TRUE)
df <- df_unfiltered[apply(df_unfiltered[, all_n_cols, drop = FALSE] >= opt$min_ind, 1, all), , drop = FALSE]

# Ensure polymorphic sites only (drop sites fixed at 0 or 1 across all samples)
poly_mask <- apply(df[, af_cols, drop = FALSE], 1, function(x) !all(x == 0) && !all(x == 1))
df <- df[poly_mask, , drop = FALSE]

message(sprintf("Total polymorphic SNPs passing filters: %d", nrow(df)))

# Iterative selection scan
iterative_result <- run_iterative_selection(
  df = df, 
  regions = regions, 
  ne_generations = ne_generations,
  test_generations = test_generations, 
  fdr_cutoff = opt$fdr_cutoff, 
  max_rounds = opt$max_rounds
)

iteration_summary <- build_iteration_summary(iterative_result, regions)
final_outputs <- build_final_outputs(iterative_result$final, regions)

write.table(iteration_summary, file.path(opt$out_dir, "iteration_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(final_outputs$test_results, file.path(opt$out_dir, "final_test_results.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(final_outputs$cmh_results_full, file.path(opt$out_dir, "final_cmh_results_full.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# Bootstrap Ne
final_selected_idx <- iterative_result$final$selected_idx
df_neutral_final <- if (length(final_selected_idx) > 0) df[-final_selected_idx, , drop = FALSE] else df
neutral_data_final <- build_region_data(df_neutral_final, regions)
ne_boot <- bootstrap_ne_by_region(neutral_data_final$mafs, neutral_data_final$covs, regions, generations = ne_generations, n_boot = opt$n_boot)
write.table(ne_boot, file.path(opt$out_dir, "ne_bootstrap.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# -----------------------------
# Chi-squared Manhattan plots
# -----------------------------
chisq_dat <- final_outputs$test_results
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
  dat <- final_outputs$cmh_results_full
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
      message("CMH Manhattan plot saved successfully.")
    }
  }
} else {
  message("Only one region provided; skipping CMH Manhattan plot.")
}