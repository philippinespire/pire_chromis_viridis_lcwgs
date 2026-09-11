# Helper functions for iteratively using Chi-sq and CMH tests from the ACER package
# to calculate genomic sites under selection and Ne
# Used by run_acer.r
# Original by Marianne Dehasque 2026
# Updated to handle one population by Malin Pinsky, July 2026, with the assistance of Google Gemini.


build_region_data <- function(df, regions) {
  mafs <- list()
  covs <- list()

  for (r in names(regions)) {
    pops <- regions[[r]]
    af_cols <- paste0(pops, "_AF")
    n_cols <- paste0(pops, "_N")

    mafs[[r]] <- df[, af_cols, drop = FALSE]
    covs[[r]] <- df[, n_cols, drop = FALSE]

    colnames(mafs[[r]]) <- c("A", "C")
    colnames(covs[[r]]) <- c("A", "C")
  }

  list(mafs = mafs, covs = covs)
}

# Calculates temporal Ne using Jorde & Ryman (2007) unbiased Fs'
jrNe2 <- function(maf1, maf2, n1, n2, generations) {
  Fsnum <- (maf1 - maf2)^2 + (1 - maf1 - (1 - maf2))^2

  z <- (maf1 + maf2) / 2
  z2 <- ((1 - maf1) + (1 - maf2)) / 2
  Fsdenom <- z * (1 - z) + z2 * (1 - z2)
  Fs <- sum(Fsnum) / sum(Fsdenom)

  sl <- 2 / (1 / n1 + 1 / n2)

  S <- length(maf1) * 2 / sum(2 / sl)
  S2 <- length(maf2) * 2 / sum(2 / n2)
  Fsprime <- (Fs * (1 - 1 / (4 * S)) - 1 / S) / ((1 + Fs / 4) * (1 - 1 / (2 * S2)))

  generations / (2 * Fsprime)
}

calc_ne_by_region <- function(mafs, covs, regions, generations) {
  setNames(
    sapply(names(regions), function(r) {
      if (nrow(mafs[[r]]) == 0) {
        return(NA_real_)
      }

      jrNe2(
        maf1 = mafs[[r]]$A,
        maf2 = mafs[[r]]$C,
        n1 = covs[[r]]$A,
        n2 = covs[[r]]$C,
        generations = generations
      )
    }),
    names(regions)
  )
}

jrNe2boot <- function(data, indices, generations) {
  maf1 <- data$freq1[indices]
  maf2 <- data$freq2[indices]
  n1   <- data$n1[indices]
  n2   <- data$n2[indices]

  Fsnum  <- (maf1 - maf2)^2 + (1 - maf1 - (1 - maf2))^2
  z      <- (maf1 + maf2) / 2
  z2     <- ((1 - maf1) + (1 - maf2)) / 2
  Fsdenom <- z * (1 - z) + z2 * (1 - z2)
  Fs     <- sum(Fsnum) / sum(Fsdenom)

  sl     <- 2 / (1 / n1 + 1 / n2)
  S      <- length(maf1) * 2 / sum(2 / sl)
  S2     <- length(maf2) * 2 / sum(2 / n2)
  Fsprime <- (Fs * (1 - 1 / (4 * S)) - 1 / S) / ((1 + Fs / 4) * (1 - 1 / (2 * S2)))
  Ne     <- generations / (2 * Fsprime)
  if (Ne < 0) Ne <- Inf
  Ne
}

bootstrap_ne_by_region <- function(mafs, covs, regions, generations, n_boot = 1000) {
  results <- lapply(names(regions), function(r) {
    dat <- data.frame(
      freq1 = mafs[[r]]$A,
      freq2 = mafs[[r]]$C,
      n1    = covs[[r]]$A,
      n2    = covs[[r]]$C
    )

    if (nrow(dat) == 0) {
      return(data.frame(region = r, Ne = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_,
                        stringsAsFactors = FALSE))
    }

    b <- boot::boot(data = dat, statistic = jrNe2boot, R = n_boot, generations = generations)
    ci <- tryCatch(
      boot::boot.ci(b, type = "perc")$percent[4:5],
      error = function(e) c(NA_real_, NA_real_)
    )

    data.frame(region = r, Ne = b$t0, ci_lower = ci[1], ci_upper = ci[2],
               stringsAsFactors = FALSE)
  })

  do.call(rbind, results)
}

run_chisq_by_region <- function(mafs, covs, regions, ne_by_region, gen) {
  pval_results <- list()

  for (r in names(regions)) {
    pval_vec <- adapted.chisq.test(
      freq = mafs[[r]],
      coverage = covs[[r]],
      Ne = c(ne_by_region[[r]], ne_by_region[[r]]),
      gen = gen
    )

    pval_df <- data.frame(pval = as.numeric(pval_vec))
    pval_df$fdr <- p.adjust(pval_df$pval, method = "fdr")
    pval_results[[r]] <- pval_df
  }

  pval_results
}

build_timepoint_matrix <- function(region_list, timepoint) {
  out <- do.call(cbind, lapply(region_list, function(x) x[, timepoint, drop = FALSE]))
  colnames(out) <- names(region_list)
  as.matrix(out)
}

build_cmh_inputs <- function(mafs, covs) {
  all_sel_A <- build_timepoint_matrix(mafs, "A")
  all_sel_C <- build_timepoint_matrix(mafs, "C")
  all_cov_A <- build_timepoint_matrix(covs, "A")
  all_cov_C <- build_timepoint_matrix(covs, "C")

  list(
    all_sel = cbind(all_sel_A, all_sel_C),
    all_cov = cbind(all_cov_A, all_cov_C)
  )
}

run_cmh <- function(mafs, covs, ne_by_region, regions, gen) {
  cmh_inputs <- build_cmh_inputs(mafs, covs)

  cmh_pval <- adapted.cmh.test(
    freq = cmh_inputs$all_sel,
    coverage = cmh_inputs$all_cov,
    Ne = as.numeric(ne_by_region[names(regions)]),
    gen = gen,
    repl = seq_along(regions),
    order = 1
  )

  cmh_pval_df <- data.frame(cmh_pval = as.numeric(cmh_pval))
  cmh_pval_df$cmh_fdr <- p.adjust(cmh_pval_df$cmh_pval, method = "fdr")
  cmh_pval_df
}

collect_selected_idx <- function(chisq_results, cmh_pval_df, fdr_cutoff) {
  chisq_selected_idx <- unique(unlist(
    lapply(chisq_results, function(x) which(!is.na(x$fdr) & x$fdr < fdr_cutoff)),
    use.names = FALSE
  ))

  if (is.null(chisq_selected_idx)) {
    chisq_selected_idx <- integer(0)
  }

  cmh_selected_idx <- which(!is.na(cmh_pval_df$cmh_fdr) & cmh_pval_df$cmh_fdr < fdr_cutoff)

  list(
    chisq_selected_idx = chisq_selected_idx,
    cmh_selected_idx = cmh_selected_idx,
    selected_idx = sort(unique(c(chisq_selected_idx, cmh_selected_idx)))
  )
}

build_chisq_table <- function(chisq_results, regions) {
  chisq_cols <- lapply(names(regions), function(r) {
    data.frame(
      pval = chisq_results[[r]]$pval,
      fdr = chisq_results[[r]]$fdr,
      stringsAsFactors = FALSE
    )
  })

  names(chisq_cols) <- names(regions)

  col_list <- list()
  for (r in names(regions)) {
    col_list[[paste0(r, "_chisq_pval")]] <- chisq_cols[[r]]$pval
    col_list[[paste0(r, "_chisq_fdr")]] <- chisq_cols[[r]]$fdr
  }

  as.data.frame(col_list, stringsAsFactors = FALSE)
}

build_af_summary <- function(mafs, regions) {
  col_list <- list()
  for (r in names(regions)) {
    col_list[[paste0(r, "_A")]] <- mafs[[r]]$A
    col_list[[paste0(r, "_C")]] <- mafs[[r]]$C
    col_list[[paste0(r, "_delta")]] <- mafs[[r]]$C - mafs[[r]]$A
  }
  as.data.frame(col_list, stringsAsFactors = FALSE)
}

run_single_iteration <- function(df, regions, ne_generations, test_generations, fdr_cutoff, ne_by_region = NULL) {
  region_data <- build_region_data(df, regions)
  mafs <- region_data$mafs
  covs <- region_data$covs

  if (is.null(ne_by_region)) {
    ne_by_region <- calc_ne_by_region(mafs, covs, regions, generations = ne_generations)
  }
  
  chisq_results <- run_chisq_by_region(mafs, covs, regions, ne_by_region, gen = test_generations)
  
  # Only run CMH if there are multiple regions. Otherwise, return NA dataframe.
  if (length(regions) > 1) {
    cmh_pval_df <- run_cmh(mafs, covs, ne_by_region, regions, gen = test_generations)
  } else {
    cmh_pval_df <- data.frame(
      cmh_pval = rep(NA_real_, nrow(df)),
      cmh_fdr  = rep(NA_real_, nrow(df))
    )
  }
  
  selected <- collect_selected_idx(chisq_results, cmh_pval_df, fdr_cutoff = fdr_cutoff)

  list(
    df = df,
    mafs = mafs,
    covs = covs,
    ne_by_region = ne_by_region,
    chisq_results = chisq_results,
    cmh_pval_df = cmh_pval_df,
    chisq_selected_idx = selected$chisq_selected_idx,
    cmh_selected_idx = selected$cmh_selected_idx,
    selected_idx = selected$selected_idx
  )
}

run_iterative_selection <- function(df, regions, ne_generations, test_generations, fdr_cutoff, max_rounds) {
  rounds <- list()
  selected_idx <- integer(0)  # accumulated selected SNP indices

  for (i in seq_len(max_rounds)) {
    # Estimate Ne from neutral SNPs only (original df minus all accumulated selected)
    df_neutral <- if (length(selected_idx) > 0) df[-selected_idx, , drop = FALSE] else df
    neutral_data <- build_region_data(df_neutral, regions)
    ne_by_region <- calc_ne_by_region(neutral_data$mafs, neutral_data$covs, regions,
                                      generations = ne_generations)

    # Test the FULL original dataset with the updated Ne
    round_res <- run_single_iteration(
      df = df,
      regions = regions,
      ne_generations = ne_generations,
      test_generations = test_generations,
      fdr_cutoff = fdr_cutoff,
      ne_by_region = ne_by_region
    )
    round_res$n_snps_neutral <- nrow(df_neutral)
    rounds[[i]] <- round_res

    # Converge when the selected set is identical to the previous round
    new_selected_idx <- sort(unique(round_res$selected_idx))
    if (identical(new_selected_idx, selected_idx)) {
      break
    }
    selected_idx <- new_selected_idx
  }

  list(
    rounds = rounds,
    final = rounds[[length(rounds)]],
    n_rounds = length(rounds),
    n_total_removed = length(selected_idx),
    n_remaining = nrow(df) - length(selected_idx)
  )
}

build_iteration_summary <- function(iterative_result, regions) {
  round_rows <- lapply(seq_along(iterative_result$rounds), function(i) {
    rr <- iterative_result$rounds[[i]]
    out <- data.frame(
      round = i,
      n_snps_tested = nrow(rr$df),
      n_snps_ne_input = rr$n_snps_neutral,
      n_snps_chisq_fdr_lt_0_05 = length(rr$chisq_selected_idx),
      n_snps_cmh_fdr_lt_0_05 = length(rr$cmh_selected_idx),
      n_snps_selected_union = length(rr$selected_idx),
      stringsAsFactors = FALSE
    )

    for (r in names(regions)) {
      out[[paste0("Ne_", r)]] <- as.numeric(rr$ne_by_region[[r]])
    }

    out
  })

  do.call(rbind, round_rows)
}

build_final_outputs <- function(final_round, regions) {
  id_cols <- c("CHR", "SNP", "BP")
  id_cols <- id_cols[id_cols %in% names(final_round$df)]
  ids <- final_round$df[, id_cols, drop = FALSE]

  chisq_table <- build_chisq_table(final_round$chisq_results, regions)
  test_results <- cbind(ids, chisq_table, final_round$cmh_pval_df)

  af_summary <- build_af_summary(final_round$mafs, regions)
  cmh_results_full <- cbind(ids, af_summary, final_round$cmh_pval_df)

  list(
    test_results = test_results,
    cmh_results_full = cmh_results_full
  )
}