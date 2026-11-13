# Main pipeline wrapper.
# Moved from Aim1_Code/methylation_pipeline_v4_shared.R.

#' Run complete methylation analysis pipeline
#'
#' Executes the full analysis workflow: CpG filtering, nested CV with elastic
#' net and random forest, risk score creation, and model comparison for both
#' continuous and dichotomous PM2.5 outcomes.
#'
#' @param beta_matrix Numeric matrix of beta values (CpGs in rows, samples in
#' columns)
#' @param pm25_continuous Numeric vector of continuous PM2.5 exposure values
#' @param pm25_dichotomous Factor with two levels for dichotomous exposure
#' @param covariates Optional data frame of covariates (samples in rows)
#' @param literature_studies Optional named list of CpG vectors from published
#' studies
#' @param k_outer Number of outer CV folds (default 5)
#' @param k_inner Number of inner CV folds (default 5)
#' @param var_percentile Variance filter threshold (default 80 = top 20%)
#' @param pval_threshold P-value threshold for association filter (default
#' 0.05)
#' @param min_criteria Minimum criteria for CpG selection (default 2)
#' @param alpha_range Alpha values to search for elastic net (default 0 to 1
#' by 0.05)
#' @param seed Random seed for reproducibility
#' @return List containing all results for both outcomes and model comparisons
#' @export
run_complete_pipeline <- function(beta_matrix,
                                  pm25_continuous,
                                  pm25_dichotomous,
                                  covariates = NULL,
                                  literature_studies = NULL,
                                  k_outer = 5,
                                  k_inner = 5,
                                  var_percentile = 80,
                                  pval_threshold = 0.05,
                                  min_criteria = 2,
                                  alpha_range = seq(0, 1, 0.05),
                                  seed = 123) {
  set.seed(seed)
  results <- list()

  # Validate covariate dimensions
  if (!is.null(covariates)) {
    if (nrow(covariates) != ncol(beta_matrix)) {
      stop(paste0(
        "Number of rows in covariates must match number of samples (columns) ",
        "in beta_matrix"
      ))
    }
    message(paste("Including", ncol(covariates), "covariates in models"))
    message(paste("Covariates:", paste(colnames(covariates), collapse = ", ")))
  }

  # Ensure dichotomous outcome is a factor
  if (!is.null(pm25_dichotomous) && !is.factor(pm25_dichotomous)) {
    pm25_dichotomous <- as.factor(pm25_dichotomous)
  }

  # ============================================
  # PROCESS CONTINUOUS OUTCOME
  # ============================================
  message("\n=== PROCESSING CONTINUOUS OUTCOME ===")

  # Filter CpGs using multi-criteria approach
  filtered_contin <- filter_cpg_sites(beta_matrix, pm25_continuous,
    literature_studies = literature_studies,
    var_percentile = var_percentile,
    pval_threshold = pval_threshold,
    min_criteria = min_criteria
  )

  results$continuous_filtering <- filtered_contin

  # Create nested CV fold structure
  nested_folds_contin <- create_nested_folds(pm25_continuous,
    k_outer = k_outer,
    k_inner = k_inner,
    seed = seed
  )

  # Run elastic net with nested CV
  message("\nRunning Elastic Net with nested CV...")
  en_results_contin <- elastic_net_nested_cv(
    x = filtered_contin$filtered_beta,
    y = pm25_continuous,
    covariates = covariates,
    nested_folds = nested_folds_contin,
    alpha_range = alpha_range
  )

  # Run random forest with nested CV
  message("\nRunning Random Forest with nested CV...")
  rf_results_contin <- random_forest_nested_cv_tuned(
    x = filtered_contin$filtered_beta,
    y = pm25_continuous,
    covariates = covariates,
    nested_folds = nested_folds_contin
  )

  # Create methylation risk scores
  en_score_contin <- create_mps(
    beta_matrix,
    en_results_contin,
    method = "elastic_net",
    score_type = "weighted"
  )

  rf_score_contin <- create_mps(
    beta_matrix,
    rf_results_contin,
    method = "random_forest",
    score_type = "weighted"
  )

  results$continuous <- list(
    filtered_cpgs = filtered_contin$selected_cpgs,
    n_cpgs_filtered = filtered_contin$n_selected,
    filtering_details = filtered_contin,
    elastic_net = en_results_contin,
    random_forest = rf_results_contin,
    elastic_net_score = en_score_contin,
    random_forest_score = rf_score_contin
  )

  # ============================================
  # PROCESS DICHOTOMOUS OUTCOME
  # ============================================
  message("\n=== PROCESSING DICHOTOMOUS OUTCOME ===")

  # Filter CpGs using multi-criteria approach
  filtered_dichot <- filter_cpg_sites(beta_matrix, pm25_dichotomous,
    literature_studies = literature_studies,
    var_percentile = var_percentile,
    pval_threshold = pval_threshold,
    min_criteria = min_criteria
  )

  results$dichotomous_filtering <- filtered_dichot

  # Create nested CV fold structure (stratified for classification)
  nested_folds_dichot <- create_nested_folds(pm25_dichotomous,
    k_outer = k_outer,
    k_inner = k_inner,
    seed = seed
  )

  # Run elastic net with nested CV
  message("\nRunning Elastic Net with nested CV...")
  en_results_dichot <- elastic_net_nested_cv(
    x = filtered_dichot$filtered_beta,
    y = pm25_dichotomous,
    covariates = covariates,
    nested_folds = nested_folds_dichot,
    alpha_range = alpha_range
  )

  # Run random forest with nested CV
  message("\nRunning Random Forest with nested CV...")
  rf_results_dichot <- random_forest_nested_cv_tuned(
    x = filtered_dichot$filtered_beta,
    y = pm25_dichotomous,
    covariates = covariates,
    nested_folds = nested_folds_dichot
  )

  # Create methylation risk scores
  en_score_dichot <- create_mps(
    beta_matrix,
    en_results_dichot,
    method = "elastic_net",
    score_type = "weighted"
  )

  rf_score_dichot <- create_mps(
    beta_matrix,
    rf_results_dichot,
    method = "random_forest",
    score_type = "weighted"
  )

  results$dichotomous <- list(
    filtered_cpgs = filtered_dichot$selected_cpgs,
    n_cpgs_filtered = filtered_dichot$n_selected,
    filtering_details = filtered_dichot,
    elastic_net = en_results_dichot,
    random_forest = rf_results_dichot,
    elastic_net_score = en_score_dichot,
    random_forest_score = rf_score_dichot
  )

  # ============================================
  # COMPARE ALL MODELS
  # ============================================
  message("\n=== COMPARING ALL MODELS ===")

  results$comparison <- compare_all_models(
    list(
      elastic_net = results$continuous$elastic_net,
      random_forest = results$continuous$random_forest
    ),
    list(
      elastic_net = results$dichotomous$elastic_net,
      random_forest = results$dichotomous$random_forest
    )
  )

  # Compare filtering between continuous and dichotomous approaches
  results$filtering_comparison <- compare_filtering_results(
    filtered_contin, filtered_dichot
  )

  # Print final summary
  message("\n=== PIPELINE COMPLETE ===")
  message("\nModel Performance Summary:")
  print(results$comparison$comparison_table)
  message(paste("\n", results$comparison$recommendations))

  results
}
