#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(scales)
})

args <- commandArgs(trailingOnly = TRUE)
stats_file   <- args[1] # master_sequencing_stats.tsv
contig_file  <- args[2] # contig_depth_summary.tsv
amber_dir    <- args[3] # directory containing *.amber_MQ25.txt
mapd_dir     <- args[4] # directory containing mapdamage outputs
out_dir      <- args[5] # publish directory

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Define color palette for era grouping across all plots
era_colors <- c("historical" = "#D55E00", "modern" = "#0072B2")

# --- 1. Master Sequencing Stats & Individual Depth Histogram ---
if (file.exists(stats_file) && file.info(stats_file)$size > 0) {
  stats <- read_tsv(stats_file, show_col_types = FALSE)
  
  # Write formatted table with requested columns
  write_csv(stats, file.path(out_dir, "individual_sequencing_summary_table.csv"))
  
  # Plot 1: Histogram of Average Depth by Individual
  p1 <- stats %>%
    filter(!is.na(total_cov) & total_cov > 0) %>%
    ggplot(aes(x = total_cov, fill = era)) +
    geom_histogram(bins = 30, alpha = 0.7, position = "identity", color = "black") +
    scale_fill_manual(values = era_colors) +
    theme_bw() +
    labs(title = "Average Depth Distribution by Individual", x = "Mean Depth (X) [Log10 Scale]", y = "Count")
  
  ggsave(file.path(out_dir, "depth_by_individual_histogram.png"), p1, width = 7, height = 5)
}

# --- 2. Histogram of Average Depth by Contig ---
if (file.exists(contig_file) && file.info(contig_file)$size > 0) {
  c_depth <- read_tsv(contig_file, show_col_types = FALSE)
  write_tsv(c_depth, file.path(out_dir, "contig_depth_summary.tsv"))

  c_depth_filtered <- c_depth %>% filter(!is.na(depth) & depth > 0)
  
  p2 <- ggplot(c_depth_filtered, aes(x = depth, fill = era)) +
    geom_histogram(bins = 50, alpha = 0.7, position = "identity", color = "black") +
    scale_x_log10(
      breaks = c(0.01, 0.1, 0.5, 1, 2, 5, 10, 50, 100, 500, 1000, 2500),
      labels = scales::label_comma(drop0trailing = TRUE)
    ) +
    scale_fill_manual(values = era_colors) +
    theme_bw() +
    labs(
      title = "Average Depth Distribution across Contigs",
      x = "Mean Depth (X) [Log10 Scale]",
      y = "Contig Count",
      fill = "Era"
    ) +
    theme(axis.text.x = element_text(angle = -45, hjust = 0))
  
  ggsave(file.path(out_dir, "depth_by_contig_histogram.png"), p2, width = 7, height = 5)
}

# --- 3 & 4. AMBER Summary Plots ---
amber_files <- list.files(amber_dir, pattern = "\\.txt$", full.names = TRUE)

if (length(amber_files) > 0) {
  
  parse_single_amber <- function(filepath) {
    lines <- readLines(filepath)
    sep_idx <- grep("^---", lines)
    
    # Extract sample name from header or filename
    sample_line <- lines[grep("^sample:", lines)]
    sample_id <- if (length(sample_line) > 0) trimws(sub("^sample:\\s*", "", sample_line[1])) else sub("\\.amber_MQ25\\.txt$", "", basename(filepath))
    
    get_block <- function(header_pattern) {
      h_idx <- grep(header_pattern, lines)
      if (length(h_idx) == 0) return(NULL)
      next_sep <- sep_idx[sep_idx > h_idx[1]][1]
      if (is.na(next_sep)) next_sep <- length(lines) + 1
      
      block_text <- lines[h_idx[1]:(next_sep - 1)]
      block_text <- block_text[nchar(trimws(block_text)) > 0]
      df <- read.table(text = paste(block_text, collapse = "\n"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
      df$sample <- sample_id
      return(df)
    }
    
    list(
      mismatch = get_block("^READ_LENGTH\tMISMATCH_RATE"),
      counts   = get_block("^READ_LENGTH\tREAD_COUNTS"),
      dist     = get_block("^DISTANCE_FROM_READ_END")
    )
  }
  
  parsed_amber <- lapply(amber_files, parse_single_amber)
  counts_df   <- bind_rows(lapply(parsed_amber, `[[`, "counts"))
  mismatch_df <- bind_rows(lapply(parsed_amber, `[[`, "mismatch"))
  dist_df     <- bind_rows(lapply(parsed_amber, `[[`, "dist"))

  # Plot 3: AMBER Read Length Distribution (% Reads)
  if (!is.null(counts_df) && nrow(counts_df) > 0) {
    counts_plot_data <- counts_df %>%
      group_by(sample) %>%
      mutate(percent_reads = (READ_COUNTS / sum(READ_COUNTS)) * 100) %>%
      ungroup() %>%
      left_join(stats %>% select(sample, era) %>% unique(), by = "sample")

    p3 <- ggplot(counts_plot_data, aes(x = READ_LENGTH, y = percent_reads, group = sample, color = era)) +
      geom_line(alpha = 0.6, linewidth = 0.7) +
      scale_color_manual(values = era_colors) +
      theme_bw() +
      labs(title = "AMBER: Read Length vs % Reads", x = "Read Length (bp)", y = "% Reads", color = "Era")
    
    ggsave(file.path(out_dir, "amber_read_length_distribution.png"), p3, width = 8, height = 5)
  }
  
  # Plot 4: Read-End Mismatch Trends (3-Panel)
  if (!is.null(dist_df) && nrow(dist_df) > 0) {
    dist_long <- dist_df %>%
      pivot_longer(
        cols = c("C.to.T", "CpG.to.TG", "other"),
        names_to = "mutation_type",
        values_to = "mismatch_freq"
      ) %>%
      mutate(
        mutation_type = case_when(
          mutation_type %in% c("C.to.T", "C_to_T") ~ "C to T",
          mutation_type %in% c("CpG.to.TG", "CpG_to_TG") ~ "CpG to TG",
        TRUE ~ "Other"
        ),
        mutation_type = factor(mutation_type, levels = c("C to T", "CpG to TG", "Other"))
      ) %>%
      left_join(stats %>% select(sample, era) %>% unique(), by = "sample")

    p4 <- ggplot(dist_long, aes(x = DISTANCE_FROM_READ_END, y = mismatch_freq, group = sample, color = era)) +
      geom_line(alpha = 0.5) +
      facet_wrap(~ mutation_type, scales = "fixed", nrow = 1) +
      scale_color_manual(values = era_colors) +
      theme_bw() +
      labs(title = "AMBER: Mismatch Frequency vs Distance from Read End", x = "Distance from Read End (bp)", y = "Frequency", color = "Era")
    
    ggsave(file.path(out_dir, "amber_mismatch_frequencies.png"), p4, width = 11, height = 4.5)
  }
}

# --- 5. mapDamage Overlay Plots ---
mapd_files <- list.files(mapd_dir, pattern = "misincorporation\\.txt$", full.names = TRUE)

if (length(mapd_files) > 0) {
  parse_single_mapdamage <- function(f) {
    tryCatch({
      fname <- basename(f)
      sample_name <- sub("_.*", "", fname)
      
      # Read tab-delimited table, ignoring lines starting with '#'
      df <- read.table(f, header = TRUE, comment.char = "#", sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
      df$sample <- sample_name
      return(df)
    }, error = function(e) NULL)
  }
  
  mapd_list <- lapply(mapd_files, parse_single_mapdamage)
  mapd_df <- bind_rows(mapd_list)
  
  if (!is.null(mapd_df) && nrow(mapd_df) > 0 && all(c("Pos", "End") %in% colnames(mapd_df))) {
    if (exists("stats") && "era" %in% colnames(stats)) {
      mapd_df <- left_join(mapd_df, stats %>% select(sample, era) %>% unique(), by = "sample")
    } else {
      mapd_df$era <- "unknown"
    }
    
    # Identify substitution columns (e.g., C>T, G>A)
    sub_cols <- intersect(c("C>T", "G>A", "C_T", "G_A", "C.T", "G.A"), colnames(mapd_df))
    
    if (length(sub_cols) > 0) {
    mapd_long <- mapd_df %>%
        pivot_longer(cols = all_of(sub_cols), names_to = "mutation", values_to = "count") %>%
        mutate(
          total_bases = ifelse(End == "5p", C + count, G + count),
          frequency   = ifelse(total_bases > 0, count / total_bases, 0),
          panel       = case_when(
            mutation %in% c("C>T", "C_T", "C.T") ~ "C to T (5' Deamination)",
            mutation %in% c("G>A", "G_A", "G.A") ~ "G to A (3' Deamination)",
            TRUE ~ "Other Substitutions"
          )
        )
    
      p5 <- ggplot(mapd_long, aes(x = Pos, y = frequency, group = interaction(sample, End), color = era)) +
        geom_line(alpha = 0.5) +
      facet_wrap(~panel, scales = "free_y") +
        scale_color_manual(values = era_colors) +
      theme_bw() +
        labs(title = "mapDamage: Mutation Frequencies vs Distance from Ends", x = "Distance from End (bp)", y = "Frequency", color = "Era")
    
      ggsave(file.path(out_dir, "mapdamage_overlay.png"), p5, width = 12, height = 5)
    }
  }
}

# --- Plot 6: 12-Subpanel Summary Boxplots by Era ---
  metric_cols <- c("total_seqs", "total_mapped", "total_uniq", "total_MQ25",
                   "endogenous", "complexity", "grr", "total_cov",
                   "read_min", "read_max", "read_median", "read_mean")
  
  if (all(metric_cols %in% colnames(stats))) {
    stats_long <- stats %>%
      pivot_longer(
        cols = all_of(metric_cols),
        names_to = "metric",
        values_to = "value"
      ) %>%
      mutate(
        metric = factor(metric, levels = metric_cols),
        era    = factor(era, levels = c("historical", "modern"))
      )

    p6 <- ggplot(stats_long, aes(x = era, y = value, fill = era)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.5) +
      geom_jitter(color = "grey40", width = 0.15, size = 1.2, alpha = 0.7) +
      facet_wrap(~ metric, scales = "free_y", ncol = 4) +
      scale_fill_manual(values = era_colors) +
      theme_bw(base_size = 11) +
      labs(
        title = "Individual Sequencing Summary Metrics by Era",
        x = "Era",
        y = "Value",
        fill = "Era"
      ) +
      theme(
        legend.position = "bottom",
        axis.text.x = element_text(angle = -30, hjust = 0)
      )

    ggsave(file.path(out_dir, "individual_summary_metrics_boxplots.png"), p6, width = 12, height = 9)
}