#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)
stats_file   <- args[1] # master_sequencing_stats.tsv
contig_file  <- args[2] # contig_depth_summary.tsv
amber_dir    <- args[3] # directory containing *.amber_MQ25.txt
mapd_dir     <- args[4] # directory containing mapdamage outputs
out_dir      <- args[5] # publish directory

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# --- 1. Master Sequencing Stats & Individual Depth Histogram ---
if (file.exists(stats_file) && file.info(stats_file)$size > 0) {
  stats <- read_tsv(stats_file, show_col_types = FALSE)
  
  # Write formatted table with requested columns
  write_csv(stats, file.path(out_dir, "individual_sequencing_summary_table.csv"))
  
  # Plot 1: Histogram of Average Depth by Individual
  p1 <- ggplot(stats, aes(x = total_cov, fill = era)) +
    geom_histogram(bins = 30, alpha = 0.7, position = "identity", color = "black") +
    scale_fill_manual(values = c("historical" = "#D55E00", "modern" = "#0072B2")) +
    theme_bw() +
    labs(title = "Average Depth Distribution by Individual", x = "Mean Depth (X)", y = "Count")
  
  ggsave(file.path(out_dir, "plot1_depth_by_individual_histogram.pdf"), p1, width = 7, height = 5)
}

# --- 2. Histogram of Average Depth by Contig ---
if (file.exists(contig_file) && file.info(contig_file)$size > 0) {
  c_depth <- read_tsv(contig_file, col_names = c("contig", "mean_depth"), show_col_types = FALSE)
  
  p2 <- ggplot(c_depth, aes(x = mean_depth)) +
    geom_histogram(bins = 50, fill = "#009E73", color = "black", alpha = 0.7) +
    theme_bw() +
    labs(title = "Average Depth Distribution across Contigs", x = "Mean Depth (X)", y = "Contig Count")
  
  ggsave(file.path(out_dir, "plot2_depth_by_contig_histogram.pdf"), p2, width = 7, height = 5)
}

# --- 4 & 5. AMBER Summary Plots ---
amber_files <- list.files(amber_dir, pattern = "\\.amber_MQ25\\.txt$", full.names = TRUE)
if (length(amber_files) > 0) {
  amber_data_list <- lapply(amber_files, function(f) {
    tryCatch({
      sample_name <- gsub("\\.amber_MQ25\\.txt$", "", basename(f))
      df <- read_tsv(f, show_col_types = FALSE)
      df$sample <- sample_name
      return(df)
    }, error = function(e) NULL)
  })
  
  amber_df <- bind_rows(amber_data_list)
  
  if (exists("stats") && "era" %in% colnames(stats)) {
    amber_df <- left_join(amber_df, stats %>% select(sample, era), by = "sample")
  } else {
    amber_df$era <- "unknown"
  }
  
  # Plot 4: AMBER Read Length vs % Reads
  if (all(c("read_length", "pct_reads") %in% colnames(amber_df))) {
    p4 <- ggplot(amber_df, aes(x = read_length, y = pct_reads, group = sample, color = era)) +
      geom_line(alpha = 0.5) +
      scale_color_manual(values = c("historical" = "#D55E00", "modern" = "#0072B2")) +
      theme_bw() +
      labs(title = "AMBER: Read Length vs % Reads", x = "Read Length (bp)", y = "% Reads")
    
    ggsave(file.path(out_dir, "plot4_amber_read_length_distribution.pdf"), p4, width = 8, height = 5)
  }
  
  # Plot 5: Mismatch Frequency vs Distance from Ends (Panelled via facet_wrap)
  if (all(c("position", "mismatch_freq", "mutation_type") %in% colnames(amber_df))) {
    amber_mismatch <- amber_df %>%
      mutate(category = case_when(
        mutation_type %in% c("CpG>TpG", "CpG_to_TpG") ~ "CpG to TpG",
        mutation_type %in% c("C>T", "C_to_T") ~ "C to T",
        TRUE ~ "Other"
      ))
    
    p5 <- ggplot(amber_mismatch, aes(x = position, y = mismatch_freq, group = sample, color = era)) +
      geom_line(alpha = 0.4) +
      facet_wrap(~category, scales = "free_y") +
      scale_color_manual(values = c("historical" = "#D55E00", "modern" = "#0072B2")) +
      theme_bw() +
      labs(title = "AMBER: Mismatch Frequency vs Distance from End", x = "Distance from End (bp)", y = "Frequency")
    
    ggsave(file.path(out_dir, "plot5_amber_mismatch_frequencies.pdf"), p5, width = 10, height = 4)
  }
}

# --- 6. mapDamage Overlay Plots ---
mapd_files <- list.files(mapd_dir, pattern = "5pptable\\.txt|3pptable\\.txt", full.names = TRUE)
if (length(mapd_files) > 0) {
  mapd_list <- lapply(mapd_files, function(f) {
    tryCatch({
      fname <- basename(f)
      sample_name <- sub("_.*", "", fname)
      df <- read_tsv(f, show_col_types = FALSE)
      df$sample <- sample_name
      return(df)
    }, error = function(e) NULL)
  })
  
  mapd_df <- bind_rows(mapd_list)
  if (!is.null(mapd_df) && nrow(mapd_df) > 0) {
    if (exists("stats") && "era" %in% colnames(stats)) {
      mapd_df <- left_join(mapd_df, stats %>% select(sample, era), by = "sample")
    } else {
      mapd_df$era <- "unknown"
    }
    
    # Categorize mutation/damage panels
    mapd_long <- mapd_df %>%
      pivot_longer(cols = -c(sample, era, Pos, End), names_to = "mutation", values_to = "frequency") %>%
      mutate(panel = case_when(
        mutation %in% c("C>T", "C_T") ~ "C to T",
        mutation %in% c("G>A", "G_A") ~ "G to A",
        mutation %in% c("SoftClip", "Soft_Clip") ~ "Soft-clipped bases",
        mutation %in% c("Deletions", "Del") ~ "Deletions relative to reference",
        mutation %in% c("Insertions", "Ins") ~ "Insertions relative to reference",
        TRUE ~ "All other substitutions"
      ))
    
    p6 <- ggplot(mapd_long, aes(x = Pos, y = frequency, group = interaction(sample, End), color = era)) +
      geom_line(alpha = 0.4) +
      facet_wrap(~panel, scales = "free_y") +
      scale_color_manual(values = c("historical" = "#D55E00", "modern" = "#0072B2")) +
      theme_bw() +
      labs(title = "mapDamage: Mutation Frequencies vs Distance from Ends", x = "Distance (bp)", y = "Frequency")
    
    ggsave(file.path(out_dir, "plot6_mapdamage_overlay.pdf"), p6, width = 12, height = 7)
  }
}