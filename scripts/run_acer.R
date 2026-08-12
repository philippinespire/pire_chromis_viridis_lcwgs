#!/usr/bin/env Rscript

# Script to iteratively use Chi-sq and CMH tests from the ACER package
# to calculate genomic sites under selection and Ne.
# Sample size is calculated from the number of individuals in the ANGSD input files or Beagle file.
# Needs acer_helpers.R
# Original by Marianne Dehasque 2026
# Updated to handle one or more population with ANGSD input by Malin Pinsky, July 2026, 
#   with the assistance of Google Gemini.

# USE INSTRUCTIONS
# 
# On Wahab, probably need to run:
#   module load container_env R
#
# -----------------------------------------------------------------------------
# MODE 1: MAF Files Input
# -----------------------------------------------------------------------------
# crun Rscript run_acer.R \
#  --hist_mafs=/path/to/hist1.mafs.gz \
#  --mod_mafs=/path/to/mod1.mafs.gz \
#  --region_names=Pop1 \
#  --out_dir=/path/to/output_directory \
#  --helpers=scripts/acer_helpers.R \
#  --ne_generations=114 \
#  --test_gen_start=0 \
#  --test_gen_end=113 \
#  --fdr_cutoff=0.05 \
#  --max_rounds=20 \
#  --n_boot=1000 \
#  --min_ind=4
#
# -----------------------------------------------------------------------------
# MODE 2: Beagle + Sample Metadata CSV Input
# -----------------------------------------------------------------------------
# crun Rscript run_acer.R \
#  --beagle=/path/to/all_samples.beagle.gz \
#  --sample_info=/path/to/samples.csv \
#  --ind_col=sample \         # Column name for sample ID (default: sample)
#  --group_col=era \       # Column name for era/group (default: era)
#  --region_col=region \     # Column name for region (optional)
#  --hist_group=historic \   # Group label for historic samples (default: historic)
#  --mod_group=modern \      # Group label for modern samples (default: modern)
#  --out_dir=/path/to/output_directory \
#  --helpers=scripts/acer_helpers.R \
#  --ne_generations=114 \
#  --test_gen_start=0 \
#  --test_gen_end=113 \
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

opt <- list(
  # Mode 1: MAF input arguments
  hist_mafs = NULL,
  mod_mafs = NULL,
  region_names = NULL,
  
  # Mode 2: Beagle input arguments
  beagle = NULL,
  sample_info = NULL,
  ind_col = "sample",
  group_col = "group",
  region_col = "region",
  hist_group = "historic",
  mod_group = "modern",

  # Pipeline options
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

# Parse both --key=value and --key value styles
if (length(args) > 0) {
  i <- 1
  while (i <= length(args)) {
    arg <- args[i]
    if (grepl("^--", arg)) {
      if (grepl("=", arg)) {
        parts <- strsplit(sub("^--", "", arg), "=")[[1]]
        key <- parts[1]
        val <- paste(parts[-1], collapse = "=")
        i <- i + 1
      } else {
        key <- sub("^--", "", arg)
        val <- args[i + 1]
        i <- i + 2
      }
      
      if (key %in% names(opt)) {
        if (is.numeric(opt[[key]])) {
          opt[[key]] <- as.numeric(val)
        } else {
          opt[[key]] <- val
        }
      } else {
        warning(paste("Unknown argument ignored:", key))
      }
    } else {
      i <- i + 1
    }
  }
}

# Determine input mode
has_mafs <- !is.null(opt$hist_mafs) && !is.null(opt$mod_mafs) && !is.null(opt$region_names)
has_beagle <- !is.null(opt$beagle) && !is.null(opt$sample_info)

if (!has_mafs && !has_beagle) {
  stop("Error: Must provide either (--hist_mafs, --mod_mafs, and --region_names) OR (--beagle and --sample_info).", call.=FALSE)
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

merged_regions_list <- list()
regions <- list()

# -------------------------------------------------------------
# INPUT PROCESSING MODE 1: Legacy MAF files
# -------------------------------------------------------------
if (has_mafs) {
  message("Input mode: MAF files")
  hist_files <- strsplit(opt$hist_mafs, ",")[[1]]
  mod_files  <- strsplit(opt$mod_mafs, ",")[[1]]
  reg_names  <- strsplit(opt$region_names, ",")[[1]]

  if (length(hist_files) != length(mod_files) || length(hist_files) != length(reg_names)) {
    stop("Error: --hist_mafs, --mod_mafs, and --region_names must have the same number of items.", call.=FALSE)
  }

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
    
    valid_sites <- (merged$major_H == merged$major_M & merged$minor_H == merged$minor_M) | 
                   (merged$major_H == merged$minor_M & merged$minor_H == merged$major_M)
    merged <- merged[valid_sites, , drop = FALSE]
    
    is_flipped <- merged$major_H == merged$minor_M & merged$minor_H == merged$major_M
    merged$AF_mod_adj <- merged[[paste0(af_col_mod, "_M")]]
    merged$AF_mod_adj[is_flipped] <- 1 - merged$AF_mod_adj[is_flipped]
    
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

# -------------------------------------------------------------
# INPUT PROCESSING MODE 2: Beagle + Sample Metadata CSV
# -------------------------------------------------------------
} else {
  message("Input mode: Beagle + CSV metadata")
  
  if (!file.exists(opt$sample_info)) stop(paste("Sample info file not found:", opt$sample_info))
  meta <- read.csv(opt$sample_info, stringsAsFactors = FALSE)
  
  if (!opt$ind_col %in% names(meta)) stop(paste("Column", opt$ind_col, "not found in sample_info CSV."))
  if (!opt$group_col %in% names(meta)) stop(paste("Column", opt$group_col, "not found in sample_info CSV."))
  
  if (opt$region_col %in% names(meta)) {
    reg_names <- unique(meta[[opt$region_col]])
  } else if (!is.null(opt$region_names)) {
    reg_names <- strsplit(opt$region_names, ",")[[1]]
  } else {
    reg_names <- "Pop1"
    meta$region_temp <- "Pop1"
    opt$region_col <- "region_temp"
  }
  
  message(paste("Reading Beagle file:", opt$beagle))
  beagle_con <- if (grepl("\\.gz$", opt$beagle)) gzfile(opt$beagle) else file(opt$beagle)
  beagle_df <- read.table(beagle_con, header = TRUE, stringsAsFactors = FALSE)
  
  # Parse CHR and BP
  marker_split <- do.call(rbind, strsplit(beagle_df$marker, "_(?=[^_]+$)", perl = TRUE))
  if (ncol(marker_split) == 2) {
    chr_vec <- marker_split[, 1]
    bp_vec  <- as.numeric(marker_split[, 2])
  } else {
    chr_vec <- beagle_df$marker
    bp_vec  <- seq_len(nrow(beagle_df))
  }
  
  # Extract individual sample IDs or total count from Beagle headers
  gl_cols <- colnames(beagle_df)[4:ncol(beagle_df)]
  if (length(gl_cols) %% 3 != 0) stop("Error: Beagle genotype likelihood columns are not a multiple of 3.")
  n_beagle_ind <- length(gl_cols) / 3
  beagle_samples <- unique(sub("[_\\.][012]$", "", gl_cols))
  
  # Check if sample IDs match or if fallback to 1-to-1 order is required
  use_positional_matching <- FALSE
  if (!any(meta[[opt$ind_col]] %in% beagle_samples)) {
    message("Notice: Beagle header contains generic IDs (e.g. Ind0, Ind1). Falling back to positional order matching.")
    if (nrow(meta) != n_beagle_ind) {
      stop(sprintf("Error: CSV contains %d rows, but Beagle file contains %d individuals.", nrow(meta), n_beagle_ind))
    }
    use_positional_matching <- TRUE
    meta$ind_index <- seq_len(nrow(meta))
  }
  
  # Calculation function using sample 1-based indices
  calc_af_from_indices <- function(df_gl, ind_indices) {
    if (length(ind_indices) == 0) return(list(AF = rep(NA, nrow(df_gl)), N = rep(0, nrow(df_gl))))
    
    tot_dosage <- numeric(nrow(df_gl))
    n_ind      <- numeric(nrow(df_gl))
    
    for (k in ind_indices) {
      col_start <- 3 + (k - 1) * 3 + 1
      p0 <- df_gl[[col_start]]
      p1 <- df_gl[[col_start + 1]]
      p2 <- df_gl[[col_start + 2]]
      
      is_valid <- (abs(p0 - 1/3) > 1e-4) | (abs(p1 - 1/3) > 1e-4) | (abs(p2 - 1/3) > 1e-4)
      dosage <- p1 + 2 * p2
      
      tot_dosage <- tot_dosage + ifelse(is_valid, dosage, 0)
      n_ind      <- n_ind + as.integer(is_valid)
    }
    
    af <- ifelse(n_ind > 0, tot_dosage / (2 * n_ind), NA)
    list(AF = af, N = n_ind)
  }

  for (i in seq_along(reg_names)) {
    r_name <- reg_names[i]
    message(sprintf("Processing Region: %s from Beagle data...", r_name))
    
    region_meta <- if (opt$region_col %in% names(meta)) meta[meta[[opt$region_col]] == r_name, ] else meta
    
    if (use_positional_matching) {
      hist_indices <- region_meta[tolower(region_meta[[opt$group_col]]) == tolower(opt$hist_group), "ind_index"]
      mod_indices  <- region_meta[tolower(region_meta[[opt$group_col]]) == tolower(opt$mod_group), "ind_index"]
    } else {
      hist_samples <- region_meta[tolower(region_meta[[opt$group_col]]) == tolower(opt$hist_group), opt$ind_col]
      mod_samples  <- region_meta[tolower(region_meta[[opt$group_col]]) == tolower(opt$mod_group), opt$ind_col]
      hist_indices <- match(hist_samples, beagle_samples)
      mod_indices  <- match(mod_samples, beagle_samples)
    }
    
    hist_stats <- calc_af_from_indices(beagle_df, hist_indices)
    mod_stats  <- calc_af_from_indices(beagle_df, mod_indices)
    
    reg_df <- data.frame(
      CHR = chr_vec,
      BP  = bp_vec
    )
    reg_df[[paste0(r_name, "_A_AF")]] <- hist_stats$AF
    reg_df[[paste0(r_name, "_C_AF")]] <- mod_stats$AF
    reg_df[[paste0(r_name, "_A_N")]]  <- hist_stats$N
    reg_df[[paste0(r_name, "_C_N")]]  <- mod_stats$N
    
    merged_regions_list[[i]] <- reg_df
    regions[[r_name]] <- c(A = paste0(r_name, "_A"), C = paste0(r_name, "_C"))
  }
}

# -----------------------------
# Combine and Filter Regions
# -----------------------------
message("Combining all regions into a single dataframe...")
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
      message("CMH Manhattan plot saved successfully.")
    }
  }
} else {
  message("Only one region provided; skipping CMH Manhattan plot.")
}