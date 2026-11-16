# Step 1: CpG site multi-criteria filtering.
# Moved from Aim1_Code/methylation_pipeline_v4_shared.R.

#' Create indicator variables for literature-based CpG filtering
#'
#' @param all_cpgs Character vector of all CpG names in the dataset
#' @param literature_studies Named list where each element is a character
#' vector of CpG names from a published study (e.g., Gruzieva, Wang, Parikh,
#' Broseus)
#' @return List containing: overall_indicator (binary), by_study (matrix of
#' indicators), n_studies (count per CpG), literature_cpgs (all unique
#' literature CpGs)
#' @export
create_literature_indicator <- function(all_cpgs, literature_studies) {
  # Combine all CpGs from literature into a single unique set
  all_literature_cpgs <- unique(unlist(literature_studies))

  # Create binary indicator for each study (1 = CpG present in study)
  study_indicators <- sapply(literature_studies, function(study_cpgs) {
    as.integer(all_cpgs %in% study_cpgs)
  })

  # Create overall indicator (1 = CpG appears in at least one study)
  in_any_study <- as.integer(all_cpgs %in% all_literature_cpgs)

  # Count number of studies each CpG appears in (replication score)
  n_studies <- rowSums(study_indicators)

  list(
    overall_indicator = in_any_study,
    by_study = study_indicators,
    n_studies = n_studies,
    literature_cpgs = all_literature_cpgs
  )
}

#' Filter CpG sites using three criteria: variance, association, and literature
#' support
#'
#' This function implements a multi-criteria filtering approach to reduce the
#' dimensionality of methylation data before machine learning. CpGs must meet
#' at least min_criteria of the three filtering criteria to be retained.
#'
#' @param beta_matrix Numeric matrix of beta values (CpGs in rows, samples in
#' columns)
#' @param pm25_data Numeric vector (continuous) or factor (dichotomous) of
#' PM2.5 exposure
#' @param literature_studies Optional named list of CpG vectors from published
#' studies
#' @param var_percentile Percentile threshold for variance filter (default 80
#' = top 20%)
#' @param pval_threshold P-value threshold for association filter (default
#' 0.05)
#' @param min_criteria Minimum number of criteria a CpG must meet (default 2)
#' @param n_cores Number of cores for parallel computation (default: half of
#' available cores, as set in the original analysis script)
#' @return List containing filtered beta matrix, selected CpGs, and filtering
#' statistics
#' @export
filter_cpg_sites <- function(
  beta_matrix, pm25_data,
  literature_studies = NULL,
  var_percentile = 80,
  pval_threshold = 0.05,
  min_criteria = 2,
  n_cores = max(1, floor(parallel::detectCores() / 2))
) {
  message("Starting multi-criteria CpG filtering...")
  message(paste("Input CpGs:", nrow(beta_matrix)))

  all_cpgs <- rownames(beta_matrix)
  n_cpgs <- nrow(beta_matrix)

  # Initialize matrix to track which criteria each CpG meets (rows=CpGs,
  # cols=criteria)
  criteria_matrix <- matrix(0,
    nrow = n_cpgs, ncol = 3,
    dimnames = list(
      all_cpgs,
      c("Variance", "Association", "Literature")
    )
  )

  # ============================================
  # CRITERION 1: Variance-based filtering
  # Rationale: CpGs with low variance across samples are unlikely to be
  # informative for predicting exposure differences
  # ============================================
  message("\n--- Criterion 1: Variance-based filtering ---")
  cpg_variance <- rowVars(beta_matrix, na.rm = TRUE)
  var_threshold <- quantile(cpg_variance, var_percentile / 100, na.rm = TRUE)
  criteria_matrix[, "Variance"] <- as.integer(cpg_variance > var_threshold)
  n_variance <- sum(criteria_matrix[, "Variance"])
  message(paste(
    "CpGs passing variance filter (top", 100 - var_percentile, "%): ",
    n_variance, " (", round(100 * n_variance / n_cpgs, 1), "%)",
    sep = ""
  ))

  # ============================================
  # CRITERION 2: Association-based filtering
  # Rationale: Retain CpGs that show preliminary association with PM2.5
  # Uses correlation test for continuous outcomes, t-test for dichotomous
  # ============================================
  message("\n--- Criterion 2: Association-based filtering ---")

  # Create progress bar for p-value calculations
  pb <- progress_bar$new(
    format = paste0(
      "  Computing p-values [:bar] :percent | ",
      "ETA: :eta | :current/:total CpGs"
    ),
    total = n_cpgs,
    clear = FALSE,
    width = 80
  )

  # Calculate association p-values with progress tracking
  # Using parallel processing via foreach for large datasets
  if (n_cpgs > 1000 && n_cores > 1) {
    # Parallel computation for large datasets
    message(paste("  Using parallel computation with", n_cores, "cores..."))

    # Split CpGs into chunks for parallel processing
    chunk_size <- ceiling(n_cpgs / n_cores)
    chunks <- split(1:n_cpgs, ceiling(seq_along(1:n_cpgs) / chunk_size))

    pvals <-
      foreach(chunk = chunks, .combine = c, .packages = c("stats")) %dopar% {
        chunk_pvals <- sapply(chunk, function(i) {
          cpg_values <- beta_matrix[i, ]
          if (is.numeric(pm25_data)) {
            # Equivalent to lm() test for p-value unless covariates
            tryCatch(cor.test(cpg_values, pm25_data)$p.value,
              error = function(e) 1
            ) # are present in the model. Then, use lm()
          } else {
            tryCatch(t.test(cpg_values ~ pm25_data)$p.value,
              error = function(e) 1
            )
          }
        })
        chunk_pvals
      }
    names(pvals) <- all_cpgs
    pb$terminate()
  } else {
    # Sequential computation with progress bar for smaller datasets
    pvals <- numeric(n_cpgs)
    names(pvals) <- all_cpgs

    for (i in 1:n_cpgs) {
      cpg_values <- beta_matrix[i, ]
      if (is.numeric(pm25_data)) {
        pvals[i] <- tryCatch(cor.test(cpg_values, pm25_data)$p.value,
          error = function(e) 1
        )
      } else {
        pvals[i] <- tryCatch(t.test(cpg_values ~ pm25_data)$p.value,
          error = function(e) 1
        )
      }
      pb$tick()
    }
  }

  criteria_matrix[, "Association"] <- as.integer(pvals < pval_threshold)
  n_association <- sum(criteria_matrix[, "Association"])
  message(paste("CpGs passing association filter (p < ", pval_threshold, "): ",
    n_association, " (", round(100 * n_association / n_cpgs, 1), "%)",
    sep = ""
  ))

  # ============================================
  # CRITERION 3: Literature-based filtering
  # Rationale: Prioritize CpGs previously associated with air pollution
  # in published epigenome-wide association studies
  # ============================================
  message("\n--- Criterion 3: Literature-based filtering ---")
  if (!is.null(literature_studies)) {
    lit_indicators <- create_literature_indicator(all_cpgs, literature_studies)
    criteria_matrix[, "Literature"] <- lit_indicators$overall_indicator
    n_literature <- sum(criteria_matrix[, "Literature"])
    message(paste("CpGs in external studies: ", n_literature,
      " (", round(100 * n_literature / n_cpgs, 1), "%)",
      sep = ""
    ))

    # Report breakdown by individual study
    for (study in names(literature_studies)) {
      n_study <- sum(lit_indicators$by_study[, study])
      message(paste("  -", study, ":", n_study, "CpGs"))
    }
  } else {
    message("No literature studies provided - criterion not used")
    criteria_matrix[, "Literature"] <- 0
    n_literature <- 0
    lit_indicators <- NULL
  }

  # ============================================
  # COMBINE CRITERIA
  # Select CpGs meeting at least min_criteria of the three filters
  # ============================================
  message(paste(
    "\n--- Combining criteria (at least", min_criteria, "of 3) ---"
  ))

  # Count how many criteria each CpG meets
  criteria_counts <- rowSums(criteria_matrix)

  # Select CpGs meeting minimum threshold
  selected <- which(criteria_counts >= min_criteria)
  selected_cpgs <- all_cpgs[selected]

  # Create breakdown of which criteria combinations CpGs meet
  combination_breakdown <- table(apply(criteria_matrix, 1, function(x) {
    paste(c("Var", "Assoc", "Lit")[as.logical(x)], collapse = "+")
  }))

  # Report selection statistics
  message("\nCpG selection breakdown:")
  message(paste("  - All 3 criteria:", sum(criteria_counts == 3)))
  message(paste("  - Exactly 2 criteria:", sum(criteria_counts == 2)))
  if (min_criteria <= 1) {
    message(paste("  - Exactly 1 criterion:", sum(criteria_counts == 1)))
  }
  message(paste("  - No criteria:", sum(criteria_counts == 0)))

  # Report specific two-way and three-way combinations
  message("\nSpecific combinations (2+ criteria):")
  for (combo in names(combination_breakdown)) {
    n_combo <- combination_breakdown[combo]
    if (grepl("\\+", combo) && n_combo > 0) {
      message(paste("  -", combo, ":", n_combo, "CpGs"))
    }
  }

  message(paste(
    "\nTotal CpGs selected:", length(selected_cpgs),
    "(", round(100 * length(selected_cpgs) / n_cpgs, 1), "%)"
  ))

  # Compute pairwise and three-way overlap statistics
  overlap_var_assoc <- sum(
    criteria_matrix[, "Variance"] & criteria_matrix[, "Association"]
  )
  overlap_var_lit <- sum(
    criteria_matrix[, "Variance"] & criteria_matrix[, "Literature"]
  )
  overlap_assoc_lit <- sum(
    criteria_matrix[, "Association"] & criteria_matrix[, "Literature"]
  )
  overlap_all <- sum(criteria_counts == 3)

  list(
    filtered_beta = beta_matrix[selected, ],
    selected_cpgs = selected_cpgs,
    n_selected = length(selected_cpgs),
    criteria_matrix = criteria_matrix,
    criteria_counts = criteria_counts,
    p_values = pvals,
    variance_values = cpg_variance,
    literature_info = lit_indicators,
    filters_summary = list(
      n_variance_only = sum(
        criteria_counts == 1 & criteria_matrix[, "Variance"]
      ),
      n_association_only = sum(
        criteria_counts == 1 & criteria_matrix[, "Association"]
      ),
      n_literature_only = sum(
        criteria_counts == 1 & criteria_matrix[, "Literature"]
      ),
      n_var_assoc = overlap_var_assoc,
      n_var_lit = overlap_var_lit,
      n_assoc_lit = overlap_assoc_lit,
      n_all_three = overlap_all,
      n_at_least_two = sum(criteria_counts >= 2),
      n_total = length(selected_cpgs)
    ),
    combination_breakdown = combination_breakdown
  )
}

#' Generate detailed report of CpG filtering results
#'
#' @param filter_results Output from filter_cpg_sites()
#' @param beta_matrix Optional original beta matrix for variance statistics
#' @export
analyze_filtering_results <- function(filter_results, beta_matrix = NULL) {
  cat("\n========== FILTERING ANALYSIS REPORT ==========\n\n")

  # Overall selection statistics
  cat("OVERALL SELECTION:\n")
  cat(sprintf(
    "  Total CpGs selected: %d (%.1f%%)\n",
    filter_results$n_selected,
    100 * filter_results$n_selected / nrow(filter_results$criteria_matrix)
  ))

  # Individual criteria performance
  cat("\nINDIVIDUAL CRITERIA:\n")
  for (crit in colnames(filter_results$criteria_matrix)) {
    n_pass <- sum(filter_results$criteria_matrix[, crit])
    pct_pass <- 100 * n_pass / nrow(filter_results$criteria_matrix)
    cat(sprintf("  %s: %d CpGs (%.1f%%)\n", crit, n_pass, pct_pass))
  }

  # Pairwise and three-way overlap statistics
  cat("\nOVERLAP ANALYSIS:\n")
  cat(sprintf(
    "  Variance \u2229 Association: %d\n",
    filter_results$filters_summary$n_var_assoc
  ))
  cat(sprintf(
    "  Variance \u2229 Literature: %d\n",
    filter_results$filters_summary$n_var_lit
  ))
  cat(sprintf(
    "  Association \u2229 Literature: %d\n",
    filter_results$filters_summary$n_assoc_lit
  ))
  cat(sprintf(
    "  All three criteria: %d\n",
    filter_results$filters_summary$n_all_three
  ))

  # Distribution of criteria counts across all CpGs
  cat("\nCRITERIA COUNT DISTRIBUTION:\n")
  counts_table <- table(filter_results$criteria_counts)
  for (i in names(counts_table)) {
    cat(sprintf("  %s criteria: %d CpGs\n", i, counts_table[i]))
  }

  # Variance statistics for selected CpGs (if beta matrix provided)
  if (!is.null(beta_matrix)) {
    selected_var <-
      filter_results$variance_values[filter_results$selected_cpgs]
    if (length(selected_var) > 0) {
      cat("\nSELECTED CPGs VARIANCE STATISTICS:\n")
      cat(sprintf("  Min variance: %.4f\n", min(selected_var)))
      cat(sprintf("  Median variance: %.4f\n", median(selected_var)))
      cat(sprintf("  Mean variance: %.4f\n", mean(selected_var)))
      cat(sprintf("  Max variance: %.4f\n", max(selected_var)))
    } else {
      cat("\nSELECTED CPGs VARIANCE STATISTICS:\n")
      cat("  No CpGs selected; variance statistics not available.\n")
    }
  }

  cat("\n================================================\n")
}

#' Export filtering results to files for downstream analysis
#'
#' @param filter_results Output from filter_cpg_sites()
#' @param output_dir Directory to save results (created if doesn't exist)
#' @export
export_filtering_results <- function(
  filter_results, output_dir = "filtering_results"
) {
  dir.create(output_dir, showWarnings = FALSE)

  # Create comprehensive data frame with all filtering information
  export_df <- data.frame(
    CpG = rownames(
      filter_results$criteria_matrix
    ), filter_results$criteria_matrix,
    Criteria_Count = filter_results$criteria_counts,
    Selected = as.integer(rownames(filter_results$criteria_matrix) %in%
      filter_results$selected_cpgs),
    P_Value = filter_results$p_values,
    Variance = filter_results$variance_values
  )

  # Add per-study indicators if literature data available
  if (!is.null(filter_results$literature_info)) {
    export_df <- cbind(export_df, filter_results$literature_info$by_study)
    export_df$N_Studies <- filter_results$literature_info$n_studies
  }

  # Save detailed per-CpG results
  write.csv(export_df,
    file.path(output_dir, "cpg_filtering_criteria.csv"),
    row.names = FALSE
  )

  # Save list of selected CpG names only
  write.table(filter_results$selected_cpgs,
    file.path(output_dir, "selected_cpgs.txt"),
    col.names = FALSE, row.names = FALSE, quote = FALSE
  )

  # Save summary statistics
  summary_df <- data.frame(
    Metric = names(unlist(filter_results$filters_summary)),
    Value = unlist(filter_results$filters_summary)
  )
  write.csv(summary_df,
    file.path(output_dir, "filtering_summary.csv"),
    row.names = FALSE
  )

  message(paste("\nResults exported to:", output_dir))
}
