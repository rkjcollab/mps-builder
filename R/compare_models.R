# Step 6: comprehensive model comparison.
# Moved from Aim1_Code/methylation_pipeline_v4_shared.R.

#' Compare performance across all models
#'
#' @param results_continuous List with elastic_net and random_forest results
#' for continuous outcome
#' @param results_dichotomous List with elastic_net and random_forest results
#' for dichotomous outcome
#' @return List with comparison table, best models, and recommendations
#' @export
compare_all_models <- function(results_continuous, results_dichotomous) {
  # Combine performance metrics into single table
  comparison <- rbind(
    cbind(Outcome = "Continuous", results_continuous$elastic_net$performance),
    cbind(
      Outcome = "Continuous", results_continuous$random_forest$performance
    ),
    cbind(
      Outcome = "Dichotomous", results_dichotomous$elastic_net$performance
    ),
    cbind(
      Outcome = "Dichotomous", results_dichotomous$random_forest$performance
    )
  )

  # Calculate composite scores that balance performance and stability
  # Higher R2 with lower SD is preferred for continuous
  if ("R2_Mean" %in% names(comparison)) {
    comparison$Score_Continuous <- comparison$R2_Mean -
      0.5 * comparison$R2_SD
  }

  # Higher AUC with lower SD is preferred for dichotomous
  if ("AUC_Mean" %in% names(comparison)) {
    comparison$Score_Dichotomous <- comparison$AUC_Mean -
      0.5 * comparison$AUC_SD
  }

  # Identify best model for each outcome type
  cont_models <- comparison[comparison$Outcome == "Continuous", ]
  dicho_models <- comparison[comparison$Outcome == "Dichotomous", ]

  best_continuous <- cont_models[which.max(cont_models$Score_Continuous), ]
  best_dichotomous <- dicho_models[which.max(dicho_models$Score_Dichotomous), ]

  list(
    comparison_table = comparison,
    best_continuous_model = best_continuous,
    best_dichotomous_model = best_dichotomous,
    recommendations = paste(
      "For continuous PM2.5:", best_continuous$Model,
      "\nFor dichotomous PM2.5:", best_dichotomous$Model
    )
  )
}

#' Compare CpG filtering results between continuous and dichotomous outcomes
#'
#' @param filter_cont Filtering results for continuous outcome
#' @param filter_dicho Filtering results for dichotomous outcome
#' @return List with overlap statistics and shared CpGs
#' @export
compare_filtering_results <- function(filter_cont, filter_dicho) {
  # Calculate overlap between the two filtering approaches
  overlap_cpgs <- intersect(
    filter_cont$selected_cpgs, filter_dicho$selected_cpgs
  )
  only_cont <- setdiff(filter_cont$selected_cpgs, filter_dicho$selected_cpgs)
  only_dicho <- setdiff(filter_dicho$selected_cpgs, filter_cont$selected_cpgs)

  # Calculate overlap percentage relative to smaller set
  min_sel <- min(filter_cont$n_selected, filter_dicho$n_selected)
  pct_overlap <- if (min_sel > 0) {
    100 * length(overlap_cpgs) / min_sel
  } else {
    NA_real_
  }

  comparison <- list(
    n_continuous = filter_cont$n_selected,
    n_dichotomous = filter_dicho$n_selected,
    n_overlap = length(overlap_cpgs),
    pct_overlap = pct_overlap,
    only_continuous = length(only_cont),
    only_dichotomous = length(only_dicho),
    overlap_cpgs = overlap_cpgs
  )

  message("\n--- Filtering Comparison (Continuous vs Dichotomous) ---")
  message(paste("CpGs selected for continuous:", comparison$n_continuous))
  message(paste("CpGs selected for dichotomous:", comparison$n_dichotomous))
  message(paste(
    "Overlapping CpGs:", comparison$n_overlap,
    "(", ifelse(is.na(comparison$pct_overlap),
      "NA", round(comparison$pct_overlap, 1)
    ), "% overlap)"
  ))

  comparison
}
