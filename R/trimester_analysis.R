# Trimester analysis functions.
# Moved from Aim1_Code/methylation_pipeline_v4_shared.R.

#' Run complete pipeline across all pregnancy trimesters
#'
#' @param beta_matrix Numeric matrix of beta values
#' @param pm25_data Data frame with trimester-specific PM2.5 columns
#' @param trimesters Vector of trimester names to analyze
#' @param k_outer Number of outer CV folds
#' @param k_inner Number of inner CV folds
#' @param var_percentile Variance filter threshold
#' @param pval_threshold P-value threshold
#' @param seed Random seed base (incremented for each trimester)
#' @return List with per-trimester results and cross-trimester comparison
#' @export
run_trimester_analysis <- function(
  beta_matrix, pm25_data,
  trimesters = c("trim1", "trim2", "trim3", "overall"),
  k_outer = 5, k_inner = 5,
  var_percentile = 80,
  pval_threshold = 0.05,
  seed = 123
) {
  trimester_results <- list()

  for (trim in trimesters) {
    message(paste("\n\n===== ANALYZING", toupper(trim), "====="))

    # Extract continuous PM2.5 for this trimester and create dichotomous
    # version
    pm25_cont <- pm25_data[[paste0(trim, "_avg")]]
    pm25_median <- median(pm25_cont, na.rm = TRUE)
    pm25_dicho <- factor(ifelse(pm25_cont >= pm25_median, "High", "Low"))

    # Run full pipeline with different seed for each trimester
    trimester_results[[trim]] <- run_complete_pipeline(
      beta_matrix = beta_matrix,
      pm25_continuous = pm25_cont,
      pm25_dichotomous = pm25_dicho,
      k_outer = k_outer,
      k_inner = k_inner,
      var_percentile = var_percentile,
      pval_threshold = pval_threshold,
      seed = seed + which(trimesters == trim)
    )
  }

  # Compare results across trimesters
  trimester_comparison <- compare_trimester_results(trimester_results)

  list(
    by_trimester = trimester_results,
    trimester_comparison = trimester_comparison
  )
}

#' Compare model performance and CpG selection across trimesters
#'
#' @param trimester_results Output from run_trimester_analysis()
#' @return List with performance table, best models, and consistent CpGs
#' @export
compare_trimester_results <- function(trimester_results) {
  all_metrics <- list()

  # Extract performance metrics for each trimester and model
  for (trim in names(trimester_results)) {
    trim_data <- trimester_results[[trim]]

    # Continuous outcome models
    cont_en <- trim_data$continuous$elastic_net$performance
    cont_en$Trimester <- trim
    cont_en$Outcome <- "Continuous"

    cont_rf <- trim_data$continuous$random_forest$performance
    cont_rf$Trimester <- trim
    cont_rf$Outcome <- "Continuous"

    # Dichotomous outcome models
    dicho_en <- trim_data$dichotomous$elastic_net$performance
    dicho_en$Trimester <- trim
    dicho_en$Outcome <- "Dichotomous"

    dicho_rf <- trim_data$dichotomous$random_forest$performance
    dicho_rf$Trimester <- trim
    dicho_rf$Outcome <- "Dichotomous"

    all_metrics[[trim]] <- rbind(cont_en, cont_rf, dicho_en, dicho_rf)
  }

  comparison_table <- do.call(rbind, all_metrics)

  # Identify best overall model for each outcome type
  best_continuous <- NULL
  best_dichotomous <- NULL

  if ("R2_Mean" %in% names(comparison_table)) {
    cont_subset <- comparison_table[comparison_table$Outcome == "Continuous", ]
    best_continuous <- cont_subset[which.max(cont_subset$R2_Mean), ]
  }

  if ("AUC_Mean" %in% names(comparison_table)) {
    dicho_subset <- comparison_table[
      comparison_table$Outcome == "Dichotomous",
    ]
    best_dichotomous <- dicho_subset[which.max(dicho_subset$AUC_Mean), ]
  }

  # Find CpGs consistently selected across all trimesters
  all_cpgs_en <- list()
  all_cpgs_rf <- list()

  for (trim in names(trimester_results)) {
    all_cpgs_en[[trim]] <-
      trimester_results[[trim]]$continuous$elastic_net$stable_cpgs
    # Take top 50 RF CpGs for comparison
    all_cpgs_rf[[trim]] <-
      trimester_results[[trim]]$continuous$random_forest$top_cpgs[1:50]
  }

  consistent_cpgs_en <- Reduce(intersect, all_cpgs_en)
  consistent_cpgs_rf <- Reduce(intersect, all_cpgs_rf)

  list(
    performance_table = comparison_table,
    best_continuous = best_continuous,
    best_dichotomous = best_dichotomous,
    consistent_cpgs_elastic_net = consistent_cpgs_en,
    consistent_cpgs_random_forest = consistent_cpgs_rf,
    n_consistent_en = length(consistent_cpgs_en),
    n_consistent_rf = length(consistent_cpgs_rf)
  )
}
